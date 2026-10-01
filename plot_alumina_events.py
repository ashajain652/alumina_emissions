#!/Users/ashajain/Developer/miniforge3/envs/WACCM/bin/python
"""Map the locations of alumina reentry emission events, labelled with their mass.

Each event (all rows sharing an `id`) is one dot at its lat/lon. The dot's AREA
scales with the event's total mass (summed over altitudes), and a label shows
that mass.

Input is either:
  - the user CSV: id,start_date,lat,long,altitude_km,emitted_mass_kg
  - the event netCDF from csv_to_alumina_events.py (row dim; event_id, date,
    datesec, lat, lon, altitude, mass)

Examples:
  ./plot_alumina_events.py alumina_test_emissions.csv
  ./plot_alumina_events.py alumina_test_events.nc -o events.png --color-by time
  ./plot_alumina_events.py big_inventory.csv --label top --top-n 20 --projection robinson
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Chart chrome and ink (light surface).
SURFACE = "#fcfcfb"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRID = "#e1e0d9"
LAND = "#efeee9"
MARK = "#2a78d6"  # Single series: one blue.
# Sequential blue for --color-by time (earliest light, latest dark). Starts at
# step 250 so the earliest dot still clears 2:1 against the surface.
SEQ_BLUE = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b"]

# Dot area in points^2 for the largest event; others scale linearly with mass.
MAX_AREA_PT2 = 900.0
MIN_AREA_PT2 = 64.0  # 8 px floor so small events stay visible.


def load_events(path):
    """Return one row per event: id, time, lat, lon, mass_kg, n_levels."""
    path = Path(path)
    if path.suffix.lower() == ".nc":
        import netCDF4

        with netCDF4.Dataset(path) as ds:
            date = np.asarray(ds.variables["date"][:]).astype(int)
            sec = np.asarray(ds.variables["datesec"][:]).astype(int)
            rows = pd.DataFrame({
                "id": np.asarray(ds.variables["event_id"][:]).astype(int),
                "lat": np.asarray(ds.variables["lat"][:], dtype=float),
                "lon": np.asarray(ds.variables["lon"][:], dtype=float),
                "mass": np.asarray(ds.variables["mass"][:], dtype=float),
            })
        # Keep model-calendar dates as strings; pandas cannot represent year 1.
        rows["time"] = [f"{d // 10000:04d}-{d // 100 % 100:02d}-{d % 100:02d} "
                        f"{s // 3600:02d}:{s // 60 % 60:02d}:{s % 60:02d}"
                        for d, s in zip(date, sec)]
    else:
        raw = pd.read_csv(path)
        rows = pd.DataFrame({
            "id": raw["id"].astype(int),
            "time": raw["start_date"].astype(str),
            "lat": raw["lat"].astype(float),
            "lon": raw["long"].astype(float),
            "mass": raw["emitted_mass_kg"].astype(float),
        })

    events = (rows.groupby("id", sort=False)
              .agg(time=("time", "first"), lat=("lat", "first"), lon=("lon", "first"),
                   mass_kg=("mass", "sum"), n_levels=("mass", "size"))
              .reset_index())
    # "YYYY-MM-DD HH:MM:SS" strings sort chronologically, including year 0001.
    return events.sort_values("time", kind="stable").reset_index(drop=True)


def dot_area(mass, mass_max):
    return np.maximum(MIN_AREA_PT2, MAX_AREA_PT2 * mass / mass_max)


def format_mass(kg):
    if kg >= 1e6:
        return f"{kg / 1e6:.3g} kt"
    if kg >= 1e3:
        return f"{kg / 1e3:.3g} t"
    return f"{kg:.3g} kg"


def nice_legend_masses(mass_max):
    """Three round reference masses spanning the data, for the size legend."""
    top = 10 ** np.floor(np.log10(mass_max))
    for step in (5, 2, 1):
        if step * top <= mass_max:
            top *= step
            break
    return [top / 4, top / 2, top]


def auto_central_lon(lons):
    """Centre the map so its edge falls in the widest longitude gap between events.

    Keeps dots away from the left/right edge, where they would be cut in half.
    The result is rounded to a 30° gridline when that keeps the edge in the gap.
    """
    lons = np.sort(((np.asarray(lons) + 180.0) % 360.0) - 180.0)
    if len(lons) == 0:
        return 0.0
    gaps = np.diff(np.append(lons, lons[0] + 360.0))
    k = int(np.argmax(gaps))
    lo, width = lons[k], gaps[k]
    edge = lo + width / 2.0
    central = ((edge + 180.0) + 180.0) % 360.0 - 180.0
    rounded = 30.0 * np.round(central / 30.0)
    rounded_edge = rounded + 180.0
    if (rounded_edge - lo) % 360.0 < width:  # Rounded edge still inside the gap.
        central = rounded
    return float(central)


def place_labels(fig, ax, events, to_label, areas, data_crs):
    """Place mass labels next to their dots without overlapping dots or each other.

    Tries eight positions around each dot (largest events first). If all
    collide, uses the least-crowded position pushed further out, with a leader line.
    """
    renderer = fig.canvas.get_renderer()
    to_disp = data_crs._as_mpl_transform(ax)
    pt2px = fig.dpi / 72.0

    # Obstacles in display pixels: every dot's bounding box, and the figure edge.
    obstacles = []
    for i in events.index:
        x, y = to_disp.transform((events.at[i, "lon"], events.at[i, "lat"]))
        r = np.sqrt(areas[i] / np.pi) * pt2px + 2 * pt2px  # Include the white ring.
        obstacles.append((x - r, y - r, x + r, y + r))
    fig_w, fig_h = fig.bbox.width, fig.bbox.height

    def overlap(b, boxes):
        total = 0.0
        for (x0, y0, x1, y1) in boxes:
            dx = min(b.x1, x1) - max(b.x0, x0)
            dy = min(b.y1, y1) - max(b.y0, y0)
            if dx > 0 and dy > 0:
                total += dx * dy
        # Anything outside the figure counts as heavily overlapping.
        out = (max(0, -b.x0) + max(0, b.x1 - fig_w)) * b.height + \
              (max(0, -b.y0) + max(0, b.y1 - fig_h)) * b.width
        return total + 10 * out

    # (dx, dy) unit directions and the matching text alignment.
    candidates = [(1, 1, "left", "bottom"), (-1, 1, "right", "bottom"),
                  (1, -1, "left", "top"), (-1, -1, "right", "top"),
                  (1, 0, "left", "center"), (-1, 0, "right", "center"),
                  (0, 1, "center", "bottom"), (0, -1, "center", "top")]
    halo = [pe.withStroke(linewidth=3, foreground=SURFACE)]
    placed = []
    for i in sorted(to_label, key=lambda j: -events.at[j, "mass_kg"]):
        ev = events.loc[i]
        r_pt = np.sqrt(areas[i] / np.pi)
        text = f"#{ev['id']}  {format_mass(ev['mass_kg'])}"
        best = None
        for push in (1.0, 2.2):
            for dx, dy, ha, va in candidates:
                d = (r_pt + 4) * push
                off = (dx * d * (0.75 if dy else 1), dy * d * (0.75 if dx else 1))
                ann = ax.annotate(text, xy=(ev["lon"], ev["lat"]), xycoords=to_disp,
                                  xytext=off, textcoords="offset points", ha=ha, va=va,
                                  fontsize=8.5, color=INK_PRIMARY, path_effects=halo,
                                  zorder=4)
                bb = ann.get_window_extent(renderer).expanded(1.05, 1.15)
                score = overlap(bb, obstacles + placed)
                ann.remove()
                if best is None or score < best[0]:
                    best = (score, off, ha, va, bb, push)
                if score == 0:
                    break
            if best[0] == 0:
                break
        _, off, ha, va, bb, push = best
        leader = dict(arrowstyle="-", color=INK_MUTED, linewidth=0.8,
                      shrinkA=0, shrinkB=r_pt + 2) if push > 1 else None
        ax.annotate(text, xy=(ev["lon"], ev["lat"]), xycoords=to_disp,
                    xytext=off, textcoords="offset points", ha=ha, va=va,
                    fontsize=8.5, color=INK_PRIMARY, path_effects=halo, zorder=4,
                    arrowprops=leader)
        placed.append((bb.x0, bb.y0, bb.x1, bb.y1))


def plot(events, out, projection, label_mode, top_n, color_by, title, central_lon):
    if central_lon is None:
        central_lon = auto_central_lon(events["lon"])
    proj = {"platecarree": ccrs.PlateCarree(central_longitude=central_lon),
            "robinson": ccrs.Robinson(central_longitude=central_lon)}[projection]
    data_crs = ccrs.PlateCarree()

    fig = plt.figure(figsize=(12, 6.6), facecolor=SURFACE)
    ax = fig.add_axes([0.03, 0.13, 0.94, 0.75], projection=proj)
    ax.set_global()
    ax.set_facecolor(SURFACE)
    ax.add_feature(cfeature.LAND, facecolor=LAND, edgecolor="none", zorder=0)
    ax.add_feature(cfeature.COASTLINE, edgecolor=INK_MUTED, linewidth=0.5, zorder=1)
    ax.spines["geo"].set_edgecolor(GRID)
    gl = ax.gridlines(draw_labels=True, color=GRID, linewidth=0.6,
                      xlocs=range(-180, 181, 60), ylocs=range(-90, 91, 30), zorder=1)
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {"color": INK_MUTED, "size": 8}

    mass_max = events["mass_kg"].max()
    areas = dot_area(events["mass_kg"].to_numpy(), mass_max)

    if color_by == "time":
        cmap = LinearSegmentedColormap.from_list("seq_blue", SEQ_BLUE)
        order = np.arange(len(events))  # Already sorted chronologically.
        colors = cmap(order / max(len(events) - 1, 1))
    else:
        colors = MARK

    # Draw big dots first so small ones stay on top; the surface-coloured ring
    # keeps overlapping dots distinguishable.
    draw = np.argsort(-areas)
    ax.scatter(events["lon"].to_numpy()[draw], events["lat"].to_numpy()[draw],
               s=areas[draw],
               c=colors if isinstance(colors, str) else colors[draw],
               edgecolors=SURFACE, linewidths=2, alpha=0.9,
               transform=data_crs, zorder=3)

    if label_mode == "all":
        to_label = events.index
    elif label_mode == "top":
        to_label = events["mass_kg"].nlargest(top_n).index
    else:
        to_label = []

    if len(to_label):
        fig.canvas.draw()  # Fix the map layout so label positions can be measured.
        place_labels(fig, ax, events, to_label, areas, data_crs)

    # Title and summary.
    total = events["mass_kg"].sum()
    fig.text(0.03, 0.95, title, fontsize=14, color=INK_PRIMARY, weight="bold")
    fig.text(0.03, 0.915,
             f"{len(events)} events  ·  total {format_mass(total)}  ·  "
             f"{events['time'].iloc[0]} to {events['time'].iloc[-1]}  ·  "
             f"dot area ∝ event mass (summed over altitudes)",
             fontsize=9, color=INK_SECONDARY)

    # Size legend (magnitude is the only encoding when colouring by a single hue).
    ref = nice_legend_masses(mass_max)
    handles = [plt.scatter([], [], s=dot_area(m, mass_max), color=MARK,
                           edgecolors=SURFACE, linewidths=2) for m in ref]
    leg = fig.legend(handles, [format_mass(m) for m in ref], title="Event mass",
                     loc="lower left", bbox_to_anchor=(0.03, 0.0), ncol=3, frameon=False,
                     fontsize=8.5, title_fontsize=8.5, labelcolor=INK_SECONDARY,
                     handletextpad=0.8, columnspacing=1.8, borderpad=0.8)
    leg.get_title().set_color(INK_SECONDARY)

    if color_by == "time":
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
        cax = fig.add_axes([0.62, 0.075, 0.3, 0.018])
        cb = fig.colorbar(sm, cax=cax, orientation="horizontal", ticks=[0, 1])
        cb.ax.set_xticklabels([events["time"].iloc[0], events["time"].iloc[-1]],
                              fontsize=8, color=INK_SECONDARY)
        cb.outline.set_visible(False)
        cb.set_label("Event time (earliest → latest)", fontsize=8.5, color=INK_SECONDARY)

    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", help="event CSV or event netCDF")
    p.add_argument("-o", "--output", help="output image (default: <input stem>_map.png)")
    p.add_argument("--projection", choices=["platecarree", "robinson"], default="platecarree")
    p.add_argument("--label", choices=["auto", "all", "top", "none"], default="auto",
                   help="which events get a mass label; auto = all if <= 30 events, else top N")
    p.add_argument("--top-n", type=int, default=15, help="events to label with --label top/auto")
    p.add_argument("--color-by", choices=["none", "time"], default="none",
                   help="colour dots by event time (sequential blue)")
    p.add_argument("--title", default="Alumina reentry emission events")
    p.add_argument("--central-lon", type=float, default=None,
                   help="map centre longitude (default: auto, keeps events off the map edge)")
    args = p.parse_args()

    events = load_events(args.input)
    label = args.label
    if label == "auto":
        label = "all" if len(events) <= 30 else "top"

    out = args.output or str(Path(args.input).with_name(Path(args.input).stem + "_map.png"))
    plot(events, out, args.projection, label, args.top_n, args.color_by, args.title,
         args.central_lon)

    print(f"{len(events)} events, total {events['mass_kg'].sum():.4f} kg -> {out}")
    print(events.to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
