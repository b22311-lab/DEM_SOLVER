#!/usr/bin/env python3
"""Render particle VTK time series as an animation plus sampled snapshots.

Examples
--------
python visualize_particle_time_series.py --input-dir part19_outputs/dt_0p01
python visualize_particle_time_series.py --input-dir parallel_validation_outputs/part10_n1000/omp2 --sample-count 8
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import numpy as np
from PIL import Image

try:
    import pyvista as pv
except ImportError:  # pragma: no cover - handled at runtime
    pv = None


TIME_PATTERN = re.compile(r"time\s+([\d.eE+-]+)")
FILE_PATTERN = re.compile(r"_(\d+)\.vtk$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a particle animation and sampled snapshots from output_*.vtk files."
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        type=Path,
        help="Directory containing output_*.vtk files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Directory for rendered outputs. Defaults to a sibling '<input>_viz' folder.",
    )
    parser.add_argument(
        "--scalar",
        default="speed",
        help="Scalar used to color particle spheres. Supports speed, vx, vy, vz, or an existing VTK scalar array.",
    )
    parser.add_argument(
        "--sample-count",
        type=int,
        default=6,
        help="Number of snapshot PNGs to export across the full time series.",
    )
    parser.add_argument(
        "--fps",
        type=float,
        default=8.0,
        help="Animation frame rate for the GIF output.",
    )
    parser.add_argument(
        "--window-size",
        nargs=2,
        type=int,
        default=(1200, 900),
        metavar=("WIDTH", "HEIGHT"),
        help="Render window size in pixels.",
    )
    parser.add_argument(
        "--arrow-fraction",
        type=float,
        default=0.08,
        help="Approximate fraction of the domain width occupied by the largest velocity arrow.",
    )
    parser.add_argument(
        "--sphere-resolution",
        type=int,
        default=24,
        help="Sphere theta/phi resolution used for particle glyphs.",
    )
    return parser.parse_args()


def sorted_vtk_files(folder: Path) -> list[Path]:
    files = list(folder.glob("output_*.vtk"))
    files.sort(key=lambda p: int(FILE_PATTERN.search(p.name).group(1)))
    return files


def read_time_from_header(vtk_file: Path) -> float:
    with vtk_file.open("r", encoding="utf-8") as handle:
        handle.readline()
        header = handle.readline().strip()
    match = TIME_PATTERN.search(header)
    if not match:
        raise RuntimeError(f"Could not parse time from VTK header: {vtk_file}")
    return float(match.group(1))


def ensure_scalar(mesh, scalar_name: str) -> str:
    if scalar_name in mesh.array_names:
        return scalar_name

    if "velocity" not in mesh.array_names:
        raise RuntimeError(
            f"Scalar '{scalar_name}' was requested, but the mesh does not contain a 'velocity' vector field."
        )

    velocity = np.asarray(mesh["velocity"], dtype=float)
    if scalar_name == "speed":
        mesh["speed"] = np.linalg.norm(velocity, axis=1)
        return "speed"
    if scalar_name == "vx":
        mesh["vx"] = velocity[:, 0]
        return "vx"
    if scalar_name == "vy":
        mesh["vy"] = velocity[:, 1]
        return "vy"
    if scalar_name == "vz":
        mesh["vz"] = velocity[:, 2]
        return "vz"

    raise RuntimeError(
        f"Scalar '{scalar_name}' is not present in the mesh and is not one of speed/vx/vy/vz."
    )


def scalar_title(scalar_name: str) -> str:
    titles = {
        "speed": "Speed (m/s)",
        "vx": "Vx (m/s)",
        "vy": "Vy (m/s)",
        "vz": "Vz (m/s)",
        "height": "Height (m)",
        "radius": "Radius (m)",
        "kinetic_energy": "Kinetic Energy",
    }
    return titles.get(scalar_name, scalar_name)


def pad_limits(vmin: float, vmax: float) -> tuple[float, float]:
    if math.isclose(vmin, vmax, rel_tol=1.0e-12, abs_tol=1.0e-12):
        delta = max(abs(vmin) * 0.1, 1.0e-6)
        return vmin - delta, vmax + delta
    return vmin, vmax


def combine_bounds(
    current: tuple[float, float, float, float, float, float] | None,
    new_bounds,
) -> tuple[float, float, float, float, float, float]:
    if current is None:
        return tuple(float(x) for x in new_bounds)
    return (
        min(current[0], new_bounds[0]),
        max(current[1], new_bounds[1]),
        min(current[2], new_bounds[2]),
        max(current[3], new_bounds[3]),
        min(current[4], new_bounds[4]),
        max(current[5], new_bounds[5]),
    )


def normalized_bounds(bounds, max_radius: float) -> tuple[float, float, float, float, float, float]:
    xmin, xmax, ymin, ymax, zmin, zmax = [float(x) for x in bounds]
    spans = [xmax - xmin, ymax - ymin, zmax - zmin]
    longest = max(max(spans), max_radius * 2.0, 1.0)
    axis_pad = max(max_radius * 1.25, longest * 0.02, 1.0e-3)

    if spans[0] < 1.0e-12:
        center = 0.5 * (xmin + xmax)
        xmin, xmax = center - axis_pad, center + axis_pad
    if spans[1] < 1.0e-12:
        center = 0.5 * (ymin + ymax)
        ymin, ymax = center - axis_pad, center + axis_pad
    if spans[2] < 1.0e-12:
        center = 0.5 * (zmin + zmax)
        zmin, zmax = center - axis_pad, center + axis_pad

    return xmin, xmax, ymin, ymax, zmin, zmax


def scan_series(vtk_files: list[Path], scalar_name: str) -> dict[str, object]:
    bounds = None
    max_radius = 0.0
    max_speed = 0.0
    scalar_min = math.inf
    scalar_max = -math.inf
    times: list[float] = []

    for vtk_file in vtk_files:
        mesh = pv.read(vtk_file)
        times.append(read_time_from_header(vtk_file))

        if "radius" not in mesh.array_names:
            raise RuntimeError(f"Missing required scalar 'radius' in {vtk_file}")

        radius = np.asarray(mesh["radius"], dtype=float)
        file_pad = float(np.nanmax(radius)) if radius.size else 0.0
        max_radius = max(max_radius, file_pad)

        mesh_scalar = ensure_scalar(mesh, scalar_name)
        scalar_values = np.asarray(mesh[mesh_scalar], dtype=float)
        if scalar_values.size:
            scalar_min = min(scalar_min, float(np.nanmin(scalar_values)))
            scalar_max = max(scalar_max, float(np.nanmax(scalar_values)))

        if "velocity" in mesh.array_names:
            velocity = np.asarray(mesh["velocity"], dtype=float)
            if velocity.size:
                max_speed = max(max_speed, float(np.nanmax(np.linalg.norm(velocity, axis=1))))

        expanded = (
            mesh.bounds[0] - file_pad,
            mesh.bounds[1] + file_pad,
            mesh.bounds[2] - file_pad,
            mesh.bounds[3] + file_pad,
            mesh.bounds[4] - file_pad,
            mesh.bounds[5] + file_pad,
        )
        bounds = combine_bounds(bounds, expanded)

    if bounds is None:
        raise RuntimeError("Could not determine bounds from the input VTK series.")

    clim = pad_limits(scalar_min, scalar_max)
    domain_bounds = normalized_bounds(bounds, max_radius)
    xspan = domain_bounds[1] - domain_bounds[0]
    yspan = domain_bounds[3] - domain_bounds[2]
    zspan = domain_bounds[5] - domain_bounds[4]
    longest_span = max(xspan, yspan, zspan)
    planar = zspan <= max(1.0e-9, 0.02 * max(xspan, yspan))
    arrow_factor = 0.0
    if max_speed > 1.0e-14:
        arrow_factor = longest_span * 0.08 / max_speed

    return {
        "times": times,
        "bounds": domain_bounds,
        "clim": clim,
        "max_radius": max_radius,
        "max_speed": max_speed,
        "planar": planar,
        "longest_span": longest_span,
        "arrow_factor": arrow_factor,
    }


def apply_camera(plotter, bounds, planar: bool) -> None:
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    cz = 0.5 * (zmin + zmax)
    xspan = xmax - xmin
    yspan = ymax - ymin
    zspan = zmax - zmin
    longest = max(xspan, yspan, zspan)

    if planar:
        plotter.enable_parallel_projection()
        plotter.camera_position = [
            (cx, cy, zmax + 3.0 * longest),
            (cx, cy, cz),
            (0.0, 1.0, 0.0),
        ]
        plotter.camera.parallel_scale = 0.58 * max(xspan, yspan)
    else:
        plotter.camera_position = [
            (cx + 1.7 * longest, cy - 1.5 * longest, cz + 1.2 * longest),
            (cx, cy, cz),
            (0.0, 0.0, 1.0),
        ]

    plotter.reset_camera_clipping_range()


def create_scene(
    vtk_file: Path,
    frame_label: str,
    scalar_name: str,
    scalar_label: str,
    clim: tuple[float, float],
    domain_box,
    bounds,
    arrow_factor: float,
    planar: bool,
    window_size: tuple[int, int],
    sphere_resolution: int,
):
    mesh = pv.read(vtk_file)
    mesh_scalar = ensure_scalar(mesh, scalar_name)

    sphere = pv.Sphere(radius=1.0, theta_resolution=sphere_resolution, phi_resolution=sphere_resolution)
    parts = mesh.glyph(scale="radius", geom=sphere, orient=False, factor=1.0)

    arrows = None
    if "velocity" in mesh.array_names and arrow_factor > 0.0:
        ensure_scalar(mesh, "speed")
        arrow = pv.Arrow(tip_length=0.35, tip_radius=0.08, shaft_radius=0.03)
        arrows = mesh.glyph(scale="speed", orient="velocity", geom=arrow, factor=arrow_factor)

    plotter = pv.Plotter(off_screen=True, window_size=window_size)
    plotter.set_background("white")
    plotter.add_mesh(domain_box, style="wireframe", color="black", line_width=2)
    plotter.add_mesh(
        parts,
        scalars=mesh_scalar,
        cmap="viridis",
        clim=clim,
        smooth_shading=True,
        scalar_bar_args={
            "title": scalar_label,
            "position_x": 0.05,
            "position_y": 0.04,
            "width": 0.42,
            "height": 0.08,
            "title_font_size": 14,
            "label_font_size": 11,
        },
    )

    if arrows is not None and arrows.n_points > 0:
        plotter.add_mesh(arrows, color="red", show_scalar_bar=False)

    plotter.add_text(frame_label, position="upper_left", color="black", font_size=14)
    apply_camera(plotter, bounds, planar)
    plotter.render()
    return plotter


def sample_indices(nitems: int, nsamples: int) -> list[int]:
    if nsamples <= 0:
        return []
    if nsamples >= nitems:
        return list(range(nitems))

    indices = np.linspace(0, nitems - 1, num=nsamples)
    unique = []
    seen = set()
    for value in indices:
        idx = int(round(float(value)))
        if idx not in seen:
            seen.add(idx)
            unique.append(idx)
    return unique


def save_gif(frame_paths: list[Path], gif_path: Path, fps: float) -> None:
    if not frame_paths:
        raise RuntimeError("No frames were rendered, so the animation GIF could not be created.")

    duration_ms = max(20, int(round(1000.0 / max(fps, 1.0e-6))))
    images = [Image.open(frame_path) for frame_path in frame_paths]
    try:
        images[0].save(
            gif_path,
            save_all=True,
            append_images=images[1:],
            duration=duration_ms,
            loop=0,
        )
    finally:
        for image in images:
            image.close()


def main() -> None:
    args = parse_args()
    if pv is None:
        raise SystemExit(
            "PyVista is required to run this script. Install it in the active environment, then rerun."
        )

    input_dir = args.input_dir.resolve()
    if not input_dir.is_dir():
        raise SystemExit(f"Input directory does not exist: {input_dir}")

    vtk_files = sorted_vtk_files(input_dir)
    if not vtk_files:
        raise SystemExit(f"No files matching 'output_*.vtk' were found in {input_dir}")

    output_dir = args.output_dir.resolve() if args.output_dir else input_dir.parent / f"{input_dir.name}_viz"
    frames_dir = output_dir / "frames"
    snapshots_dir = output_dir / "snapshots"
    frames_dir.mkdir(parents=True, exist_ok=True)
    snapshots_dir.mkdir(parents=True, exist_ok=True)

    scalar_name = args.scalar
    metadata = scan_series(vtk_files, scalar_name)
    domain_box = pv.Box(bounds=metadata["bounds"])
    scalar_label = scalar_title(scalar_name)
    window_size = tuple(args.window_size)

    arrow_factor = metadata["arrow_factor"]
    if metadata["max_speed"] > 1.0e-14:
        arrow_factor = metadata["longest_span"] * args.arrow_fraction / metadata["max_speed"]

    frame_paths: list[Path] = []
    for iframe, vtk_file in enumerate(vtk_files):
        time_value = metadata["times"][iframe]
        frame_label = f"t = {time_value:.5f} s   frame {iframe + 1}/{len(vtk_files)}"
        plotter = create_scene(
            vtk_file=vtk_file,
            frame_label=frame_label,
            scalar_name=scalar_name,
            scalar_label=scalar_label,
            clim=metadata["clim"],
            domain_box=domain_box,
            bounds=metadata["bounds"],
            arrow_factor=arrow_factor,
            planar=metadata["planar"],
            window_size=window_size,
            sphere_resolution=args.sphere_resolution,
        )
        frame_path = frames_dir / f"frame_{iframe:04d}.png"
        plotter.screenshot(str(frame_path))
        plotter.close()
        frame_paths.append(frame_path)

    gif_path = output_dir / "particle_animation.gif"
    save_gif(frame_paths, gif_path, args.fps)

    snapshot_ids = sample_indices(len(vtk_files), args.sample_count)
    for snapshot_no, idx in enumerate(snapshot_ids, start=1):
        vtk_file = vtk_files[idx]
        time_value = metadata["times"][idx]
        frame_label = f"sample {snapshot_no}/{len(snapshot_ids)}   t = {time_value:.5f} s"
        plotter = create_scene(
            vtk_file=vtk_file,
            frame_label=frame_label,
            scalar_name=scalar_name,
            scalar_label=scalar_label,
            clim=metadata["clim"],
            domain_box=domain_box,
            bounds=metadata["bounds"],
            arrow_factor=arrow_factor,
            planar=metadata["planar"],
            window_size=window_size,
            sphere_resolution=args.sphere_resolution,
        )
        snap_path = snapshots_dir / f"snapshot_{snapshot_no:02d}_t_{time_value:0.5f}.png"
        plotter.screenshot(str(snap_path))
        plotter.close()

    print(f"Rendered {len(frame_paths)} animation frames to {frames_dir}")
    print(f"Saved animated GIF: {gif_path}")
    print(f"Saved {len(snapshot_ids)} sampled snapshots to {snapshots_dir}")


if __name__ == "__main__":
    main()
