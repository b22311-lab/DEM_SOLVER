#!/usr/bin/env python3
"""Compare serial and OpenMP VTK outputs for numerical consistency."""

import argparse
import csv
import glob
import math
import os
import re
import sys

import numpy as np
import pyvista as pv


def sorted_vtk_files(folder):
    files = glob.glob(os.path.join(folder, "output_*.vtk"))
    files.sort(key=lambda x: int(re.search(r"_(\d+)\.vtk$", os.path.basename(x)).group(1)))
    return files


def file_key(path):
    return os.path.basename(path)


def safe_get_array(mesh, name):
    if name in mesh.array_names:
        return np.asarray(mesh[name])
    return None


def compare_dirs(serial_dir, parallel_dir):
    serial_files = sorted_vtk_files(serial_dir)
    parallel_files = sorted_vtk_files(parallel_dir)

    serial_map = {file_key(f): f for f in serial_files}
    parallel_map = {file_key(f): f for f in parallel_files}

    common_keys = sorted(set(serial_map.keys()) & set(parallel_map.keys()), key=lambda x: int(re.search(r"_(\d+)\.vtk$", x).group(1)))

    if not common_keys:
        raise RuntimeError("No matching VTK files between serial and parallel directories")

    max_pos = 0.0
    max_vel = 0.0
    max_ke = 0.0

    sumsq_pos = 0.0
    sumsq_vel = 0.0
    sumsq_ke = 0.0

    n_pos = 0
    n_vel = 0
    n_ke = 0

    for key in common_keys:
        m_s = pv.read(serial_map[key])
        m_p = pv.read(parallel_map[key])

        pos_s = np.asarray(m_s.points)
        pos_p = np.asarray(m_p.points)
        if pos_s.shape != pos_p.shape:
            raise RuntimeError(f"Point shape mismatch in {key}: {pos_s.shape} vs {pos_p.shape}")

        d_pos = pos_s - pos_p
        max_pos = max(max_pos, float(np.max(np.abs(d_pos))))
        sumsq_pos += float(np.sum(d_pos ** 2))
        n_pos += d_pos.size

        vel_s = safe_get_array(m_s, "velocity")
        vel_p = safe_get_array(m_p, "velocity")
        if vel_s is not None and vel_p is not None:
            if vel_s.shape != vel_p.shape:
                raise RuntimeError(f"Velocity shape mismatch in {key}: {vel_s.shape} vs {vel_p.shape}")
            d_vel = vel_s - vel_p
            max_vel = max(max_vel, float(np.max(np.abs(d_vel))))
            sumsq_vel += float(np.sum(d_vel ** 2))
            n_vel += d_vel.size

        ke_s = safe_get_array(m_s, "kinetic_energy")
        ke_p = safe_get_array(m_p, "kinetic_energy")
        if ke_s is not None and ke_p is not None:
            if ke_s.shape != ke_p.shape:
                raise RuntimeError(f"KE shape mismatch in {key}: {ke_s.shape} vs {ke_p.shape}")
            d_ke = ke_s - ke_p
            max_ke = max(max_ke, float(np.max(np.abs(d_ke))))
            sumsq_ke += float(np.sum(d_ke ** 2))
            n_ke += d_ke.size

    rms_pos = math.sqrt(sumsq_pos / n_pos) if n_pos > 0 else float("nan")
    rms_vel = math.sqrt(sumsq_vel / n_vel) if n_vel > 0 else float("nan")
    rms_ke = math.sqrt(sumsq_ke / n_ke) if n_ke > 0 else float("nan")

    return {
        "file_count": len(common_keys),
        "max_abs_position": max_pos,
        "max_abs_velocity": max_vel,
        "max_abs_ke": max_ke,
        "rms_position": rms_pos,
        "rms_velocity": rms_vel,
        "rms_ke": rms_ke,
    }


def main():
    parser = argparse.ArgumentParser(description="Compare serial and parallel VTK outputs")
    parser.add_argument("--serial-dir", required=True)
    parser.add_argument("--parallel-dir", required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--tol-pos", type=float, default=1.0e-5)
    parser.add_argument("--tol-vel", type=float, default=1.0e-5)
    parser.add_argument("--tol-ke", type=float, default=1.0e-7)
    parser.add_argument("--output", default="parallel_validation_report.csv")
    args = parser.parse_args()

    stats = compare_dirs(args.serial_dir, args.parallel_dir)

    status = "PASS"
    if stats["max_abs_position"] > args.tol_pos:
        status = "FAIL"
    if stats["max_abs_velocity"] > args.tol_vel:
        status = "FAIL"
    if stats["max_abs_ke"] > args.tol_ke:
        status = "FAIL"

    row = {
        "case": args.case,
        "threads": args.threads,
        **stats,
        "tol_pos": args.tol_pos,
        "tol_vel": args.tol_vel,
        "tol_ke": args.tol_ke,
        "status": status,
    }

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "case",
                "threads",
                "file_count",
                "max_abs_position",
                "max_abs_velocity",
                "max_abs_ke",
                "rms_position",
                "rms_velocity",
                "rms_ke",
                "tol_pos",
                "tol_vel",
                "tol_ke",
                "status",
            ],
        )
        writer.writeheader()
        writer.writerow(row)

    print(f"Compared {stats['file_count']} files")
    print(f"max_abs_position={stats['max_abs_position']:.6e}, rms_position={stats['rms_position']:.6e}")
    print(f"max_abs_velocity={stats['max_abs_velocity']:.6e}, rms_velocity={stats['rms_velocity']:.6e}")
    print(f"max_abs_ke={stats['max_abs_ke']:.6e}, rms_ke={stats['rms_ke']:.6e}")
    print(f"status={status}")

    if status != "PASS":
        sys.exit(2)


if __name__ == "__main__":
    main()
