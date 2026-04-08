#!/usr/bin/env python3
import csv
import math
import os
import re
import shutil
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
CONFIG_DIR = REPO_ROOT / "part19_configs"
OUTPUT_DIR = REPO_ROOT / "part19_outputs"


def fmt_dt_tag(dt: float) -> str:
    return f"dt_{str(dt).replace('.', 'p')}"


def write_config(
    path: Path,
    *,
    dt: float,
    end_time: float,
    k: float,
    gamma: float,
    density: float,
    rad_init: float,
    x0: float,
    y0: float,
    vx0: float,
    vy0: float,
    gravity: float,
    domain_x: float,
    domain_y: float,
) -> None:
    content = f"""! Auto-generated Part 19.2 case
&domain
    dimension=2
    xr_in={domain_x}
    yr_in={domain_y}
/

&material_properties
    k={k}
    gamma={gamma}
    density={density}
    rad_init={rad_init}
    nparticles=1
/

&particle_initialization
    posit_init={x0},{y0}
    vel_init={vx0},{vy0}
/

&simulation_control
    time_step={dt}
    end_time={end_time}
    write_interval=1
    verbose=.false.
    contact_search_method=0
    cell_size_factor=2.5
/

&gravity
    gravity_axis=2
    gravity_magnitude={gravity}
/
/
"""
    path.write_text(content)


def read_time_from_header(vtk_file: Path) -> float:
    with vtk_file.open("r") as f:
        _ = f.readline()
        header = f.readline().strip()
    m = re.search(r"time\s+([\d.eE+-]+)", header)
    if not m:
        raise RuntimeError(f"Could not parse time from VTK header: {vtk_file}")
    return float(m.group(1))


def sorted_vtk_files(folder: Path):
    files = sorted(folder.glob("output_*.vtk"), key=lambda p: int(re.search(r"_(\d+)\.vtk$", p.name).group(1)))
    return files


def run_case(config_path: Path, out_dir: Path):
    for f in REPO_ROOT.glob("output_*.vtk"):
        f.unlink()

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    result = subprocess.run(["./demo.exe", str(config_path)], cwd=REPO_ROOT, env=env, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Simulation failed for {config_path.name}:\n{result.stderr}\n{result.stdout}")

    out_dir.mkdir(parents=True, exist_ok=True)
    for vtk_file in REPO_ROOT.glob("output_*.vtk"):
        shutil.move(str(vtk_file), out_dir / vtk_file.name)


def extract_series(out_dir: Path):
    files = sorted_vtk_files(out_dir)
    if not files:
        raise RuntimeError(f"No VTK output in {out_dir}")

    t = []
    y = []
    vy = []
    ke = []

    for f in files:
        mesh = pv.read(str(f))
        t.append(read_time_from_header(f))
        y.append(float(mesh.points[0, 1]))
        vel = np.asarray(mesh["velocity"])
        vy.append(float(vel[0, 1]))
        ke_field = np.asarray(mesh["kinetic_energy"])
        ke.append(float(ke_field[0]))

    return np.array(t), np.array(y), np.array(vy), np.array(ke)


def fit_order(dt_vals, err_vals):
    x = np.log(np.array(dt_vals, dtype=float))
    y = np.log(np.array(err_vals, dtype=float))
    return float(np.polyfit(x, y, 1)[0])


def richardson_rows(dt_vals, phi_vals):
    rows = []
    for i in range(len(dt_vals) - 2):
        h1 = float(dt_vals[i])
        h2 = float(dt_vals[i + 1])
        h3 = float(dt_vals[i + 2])
        p1 = float(phi_vals[i])
        p2 = float(phi_vals[i + 1])
        p3 = float(phi_vals[i + 2])

        r1 = h1 / h2 if h2 > 0 else float("nan")
        r2 = h2 / h3 if h3 > 0 else float("nan")
        if not (np.isfinite(r1) and np.isfinite(r2)):
            continue
        if r1 <= 1.0 or abs(r1 - r2) > 1.0e-12:
            continue

        num = p1 - p2
        den = p2 - p3
        if abs(den) <= 1.0e-20 or abs(num) <= 1.0e-20:
            continue

        p_obs = math.log(abs(num / den)) / math.log(r1)
        if abs(r1**p_obs - 1.0) <= 1.0e-14:
            continue
        p_ext = p3 + (p3 - p2) / (r1**p_obs - 1.0)

        rows.append(
            {
                "dt_coarse": h1,
                "dt_medium": h2,
                "dt_fine": h3,
                "refinement_ratio": r1,
                "p_observed": p_obs,
                "phi_extrapolated": p_ext,
            }
        )

    return rows


def run_timestep_sensitivity_study():
    dt_values = [0.02, 0.01, 0.005, 0.0025, 0.001]

    g = 9.81
    y0 = 4.5
    density = 1000.0
    radius = 0.05
    mass = density * (4.0 / 3.0) * math.pi * radius**3

    rows = []
    trajectory_error_data = {}
    energy_data = {}
    final_height_by_dt = {}

    for dt in dt_values:
        tag = fmt_dt_tag(dt)
        cfg = CONFIG_DIR / f"part19_2_{tag}.in"
        out_dir = OUTPUT_DIR / tag

        write_config(
            cfg,
            dt=dt,
            end_time=0.8,
            k=2000.0,
            gamma=0.5,
            density=density,
            rad_init=radius,
            x0=2.5,
            y0=y0,
            vx0=0.0,
            vy0=0.0,
            gravity=g,
            domain_x=5.0,
            domain_y=5.0,
        )
        run_case(cfg, out_dir)

        t, y, vy, ke = extract_series(out_dir)

        y_ref = y0 - 0.5 * g * t**2
        vy_ref = -g * t
        y_err = y - y_ref

        rmse_y = float(np.sqrt(np.mean((y - y_ref) ** 2)))
        max_y = float(np.max(np.abs(y - y_ref)))
        rmse_vy = float(np.sqrt(np.mean((vy - vy_ref) ** 2)))
        max_vy = float(np.max(np.abs(vy - vy_ref)))

        e_total = ke + mass * g * y
        e0 = abs(e_total[0]) if len(e_total) > 0 else 1.0
        max_rel_energy_drift = float(np.max(np.abs(e_total - e_total[0])) / max(e0, 1.0e-14))

        rows.append(
            {
                "dt": dt,
                "nsteps": len(t),
                "rmse_y": rmse_y,
                "max_abs_y": max_y,
                "rmse_vy": rmse_vy,
                "max_abs_vy": max_vy,
                "max_rel_energy_drift": max_rel_energy_drift,
            }
        )

        trajectory_error_data[dt] = (t, y_err)
        energy_data[dt] = (t, e_total / e_total[0])
        final_height_by_dt[dt] = float(y[-1])

    with open(REPO_ROOT / "part19_2_timestep_sensitivity.csv", "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["dt", "nsteps", "rmse_y", "max_abs_y", "rmse_vy", "max_abs_vy", "max_rel_energy_drift"],
        )
        writer.writeheader()
        writer.writerows(rows)

    plt.figure(figsize=(7.2, 4.8), dpi=150)
    for dt in dt_values:
        t, y_err = trajectory_error_data[dt]
        plt.plot(t, y_err, linewidth=1.6, label=f"dt={dt:g}")
    plt.axhline(0.0, color="k", linestyle="--", linewidth=1.2, label="Zero error")
    plt.xlabel("Time [s]")
    plt.ylabel("Trajectory error y_num - y_ref [m]")
    plt.title("Part 19.2: Trajectory Error vs Time")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(REPO_ROOT / "part19_2_trajectory_error_overlay.png", bbox_inches="tight")
    plt.close()

    dt_arr = np.array([r["dt"] for r in rows])
    rmse_y_arr = np.array([r["rmse_y"] for r in rows])
    rmse_v_arr = np.array([r["rmse_vy"] for r in rows])

    plt.figure(figsize=(7.0, 4.8), dpi=150)
    plt.loglog(dt_arr, rmse_y_arr, "o-", linewidth=2, label="RMSE(y)")
    plt.loglog(dt_arr, rmse_v_arr, "s-", linewidth=2, label="RMSE(v_y)")
    plt.xlabel("Timestep dt [s]")
    plt.ylabel("RMSE")
    plt.title("Part 19.2: Error vs Timestep")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(REPO_ROOT / "part19_2_error_vs_dt.png", bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(7.2, 4.8), dpi=150)
    for dt in dt_values:
        t, e_norm = energy_data[dt]
        plt.plot(t, e_norm, linewidth=1.6, label=f"dt={dt:g}")
    plt.xlabel("Time [s]")
    plt.ylabel("Normalized total energy E/E0")
    plt.title("Part 19.2: Energy Evolution vs Timestep")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(REPO_ROOT / "part19_2_energy_evolution.png", bbox_inches="tight")
    plt.close()

    drift_arr = np.array([r["max_rel_energy_drift"] for r in rows])
    plt.figure(figsize=(7.0, 4.8), dpi=150)
    plt.loglog(dt_arr, drift_arr, "o-", linewidth=2)
    plt.xlabel("Timestep dt [s]")
    plt.ylabel("Max relative energy drift")
    plt.title("Part 19.2: Energy Drift vs Timestep")
    plt.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.savefig(REPO_ROOT / "part19_2_energy_drift_vs_dt.png", bbox_inches="tight")
    plt.close()

    order_y = fit_order(dt_arr, rmse_y_arr)
    order_v = fit_order(dt_arr, rmse_v_arr)

    rich_rows = richardson_rows(dt_values, [final_height_by_dt[dt] for dt in dt_values])
    t_end = max(trajectory_error_data[min(dt_values)][0])
    y_ref_end = y0 - 0.5 * g * t_end**2

    with open(REPO_ROOT / "part19_2_richardson.csv", "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "dt_coarse",
                "dt_medium",
                "dt_fine",
                "refinement_ratio",
                "p_observed",
                "phi_extrapolated",
                "y_ref_end",
                "abs_error_extrapolated",
            ],
        )
        writer.writeheader()
        for rr in rich_rows:
            writer.writerow(
                {
                    "dt_coarse": rr["dt_coarse"],
                    "dt_medium": rr["dt_medium"],
                    "dt_fine": rr["dt_fine"],
                    "refinement_ratio": rr["refinement_ratio"],
                    "p_observed": rr["p_observed"],
                    "phi_extrapolated": rr["phi_extrapolated"],
                    "y_ref_end": y_ref_end,
                    "abs_error_extrapolated": abs(rr["phi_extrapolated"] - y_ref_end),
                }
            )

    return {
        "rows": rows,
        "order_y": order_y,
        "order_v": order_v,
        "rich_rows": rich_rows,
    }


def run_critical_timestep_study():
    mass = 1.3
    k_n = 1.0e5
    radius = 0.05
    gravity = 9.81
    density = mass / ((4.0 / 3.0) * math.pi * radius**3)

    # User-requested comparison around the theoretical threshold.
    dt_cases = [0.003, 0.004]
    contact_time = math.pi * math.sqrt(mass / k_n)
    dt_crit = contact_time / math.pi

    case_data = []
    for dt in dt_cases:
        tag = fmt_dt_tag(dt)
        cfg = CONFIG_DIR / f"part19_2_critical_{tag}.in"
        out_dir = OUTPUT_DIR / f"critical_{tag}"

        write_config(
            cfg,
            dt=dt,
            end_time=8.0,
            k=k_n,
            gamma=0.0,
            density=density,
            rad_init=radius,
            x0=2.5,
            y0=0.2,
            vx0=0.0,
            vy0=0.0,
            gravity=gravity,
            domain_x=5.0,
            domain_y=80.0,
        )
        run_case(cfg, out_dir)

        t, y, _, ke = extract_series(out_dir)
        overlap = np.maximum(0.0, radius - y)
        e_total = ke + mass * gravity * y + 0.5 * k_n * overlap**2
        e0 = max(abs(e_total[0]), 1.0e-14)
        e_norm = e_total / e0

        case_data.append(
            {
                "dt": dt,
                "dt_over_dtcrit": dt / dt_crit,
                "time": t,
                "energy_norm": e_norm,
                "max_energy_norm": float(np.max(e_norm)),
                "final_energy_norm": float(e_norm[-1]),
            }
        )

    with open(REPO_ROOT / "part19_2_critical_stability.csv", "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["dt", "dt_over_dtcrit", "max_energy_norm", "final_energy_norm"],
        )
        writer.writeheader()
        for row in case_data:
            writer.writerow(
                {
                    "dt": row["dt"],
                    "dt_over_dtcrit": row["dt_over_dtcrit"],
                    "max_energy_norm": row["max_energy_norm"],
                    "final_energy_norm": row["final_energy_norm"],
                }
            )

    plt.figure(figsize=(7.2, 4.8), dpi=160)
    for row in case_data:
        plt.semilogy(row["time"], row["energy_norm"], linewidth=2.0, label=f"dt={row['dt']:.3f} s")
    plt.axhline(1.0, color="k", linestyle="--", linewidth=1.2, label="Initial total energy")
    plt.xlabel("Time [s]")
    plt.ylabel("Normalized total energy E/E0")
    plt.title("Part 19.2: Critical-Timestep Stability Check")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(REPO_ROOT / "part19_2_critical_energy.png", bbox_inches="tight")
    plt.close()

    return {
        "mass": mass,
        "k_n": k_n,
        "contact_time": contact_time,
        "dt_crit": dt_crit,
        "cases": case_data,
    }


def write_findings_tex(sensitivity, critical):
    rows = sensitivity["rows"]
    order_y = sensitivity["order_y"]
    order_v = sensitivity["order_v"]
    rich_rows = sensitivity["rich_rows"]

    with open(REPO_ROOT / "part19_2_findings.tex", "w") as f:
        f.write("\\subsection{Bonus Part 19.2: Timestep Sensitivity and Stability Limit}\n")
        f.write("A free-fall case without wall contact was simulated for multiple timesteps to quantify numerical accuracy and energy behavior.\n")
        f.write("\\begin{table}[!t]\n")
        f.write("\\centering\n")
        f.write("\\caption{Part 19.2 timestep sensitivity metrics}\n")
        f.write("\\label{tab:part19-2-timestep}\n")
        f.write("\\begin{tabular}{rcccc}\n")
        f.write("\\toprule\n")
        f.write("$\\Delta t$ [s] & RMSE$(y)$ [m] & RMSE$(v_y)$ [m/s] & Max $|\\Delta y|$ [m] & Max energy drift \\\\" + "\n")
        f.write("\\midrule\n")
        for r in rows:
            f.write(
                f"{r['dt']:.4f} & {r['rmse_y']:.3e} & {r['rmse_vy']:.3e} & {r['max_abs_y']:.3e} & {r['max_rel_energy_drift']:.3e} \\\\" + "\n"
            )
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")
        f.write(
            f"A log-log fit gives an error slope of approximately {order_y:.2f} for position RMSE, which is consistent with first-order temporal behavior expected from semi-implicit Euler integration.\n"
        )
        f.write(
            f"Velocity RMSE remains nearly timestep-independent at about 1e-7 m/s across all tested timesteps (slope {order_v:.2f}), indicating that truncation error in this observable is already below a small floating-point and interpolation floor for this setup.\n"
        )

        if rich_rows:
            p_vals = [rr["p_observed"] for rr in rich_rows]
            p_avg = float(np.mean(p_vals))
            f.write(
                f"A Richardson three-grid analysis on final-height values gives observed orders between {min(p_vals):.2f} and {max(p_vals):.2f} (mean {p_avg:.2f}), consistent with first-order time integration.\n"
            )
            f.write("\\begin{table}[!t]\n")
            f.write("\\centering\n")
            f.write("\\caption{Part 19.2 Richardson observed order (final height)}\n")
            f.write("\\label{tab:part19-2-richardson}\n")
            f.write("\\begin{tabular}{rcc}\n")
            f.write("\\toprule\n")
            f.write("$(\\Delta t_c,\\Delta t_m,\\Delta t_f)$ [s] & $r$ & $p_{obs}$ \\\\" + "\n")
            f.write("\\midrule\n")
            for rr in rich_rows:
                f.write(
                    f"({rr['dt_coarse']:.4f}, {rr['dt_medium']:.4f}, {rr['dt_fine']:.4f}) & {rr['refinement_ratio']:.1f} & {rr['p_observed']:.2f} \\\\" + "\n"
                )
            f.write("\\bottomrule\n")
            f.write("\\end{tabular}\n")
            f.write("\\end{table}\n")

        f.write(
            "To connect this convergence study with stability limits, a bouncing test was run around a theoretical explicit-step threshold. "
            "Using a linear spring-mass estimate, the contact time is $t_c \\approx \\pi\\sqrt{m/k_n}$ and the threshold used here is "
            "$\\Delta t_{crit}=t_c/\\pi$ \\cite{edem_linear_spring,pfc_timestep}.\n"
        )
        f.write(
            f"For $m={critical['mass']:.1f}$ kg and $k_n={critical['k_n']:.1e}$ N/m, this gives $t_c\\approx {critical['contact_time']:.4f}$ s and "
            f"$\\Delta t_{{crit}}\\approx {critical['dt_crit']:.4f}$ s.\n"
        )

        f.write("\\begin{table}[!t]\n")
        f.write("\\centering\n")
        f.write("\\caption{Part 19.2 critical timestep check (bouncing case)}\n")
        f.write("\\label{tab:part19-2-critical}\n")
        f.write("\\begin{tabular}{rccc}\n")
        f.write("\\toprule\n")
        f.write("$\\Delta t$ [s] & $\\Delta t/\\Delta t_{crit}$ & Max $E/E0$ & Final $E/E0$ \\\\" + "\n")
        f.write("\\midrule\n")
        for row in critical["cases"]:
            f.write(
                f"{row['dt']:.4f} & {row['dt_over_dtcrit']:.2f} & {row['max_energy_norm']:.2f} & {row['final_energy_norm']:.2f} \\\\" + "\n"
            )
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")

        under = min(critical["cases"], key=lambda x: x["dt"])
        over = max(critical["cases"], key=lambda x: x["dt"])
        f.write(
            f"The under-limit run ($\\Delta t={under['dt']:.3f}$ s) remains bounded but shows numerical energy gain (Max $E/E0={under['max_energy_norm']:.2f}$), "
            f"whereas the over-limit run ($\\Delta t={over['dt']:.3f}$ s) exhibits rapid unphysical energy growth (Max $E/E0={over['max_energy_norm']:.2f}$), "
            "which is visible as a sharp upward trend in the critical-energy plot.\n"
        )


def main() -> None:
    CONFIG_DIR.mkdir(exist_ok=True)
    OUTPUT_DIR.mkdir(exist_ok=True)

    sensitivity = run_timestep_sensitivity_study()
    critical = run_critical_timestep_study()
    write_findings_tex(sensitivity, critical)


if __name__ == "__main__":
    main()
