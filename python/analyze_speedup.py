#!/usr/bin/env python3
import csv
import math
import sys
from collections import defaultdict

import matplotlib.pyplot as plt


def to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def amdahl_speedup(threads: int, serial_fraction: float) -> float:
    if threads <= 0:
        return float("nan")
    return 1.0 / (serial_fraction + (1.0 - serial_fraction) / threads)


def estimate_serial_fraction(case_rows):
    """Estimate Amdahl serial fraction f from measured speedups."""
    estimates = []
    for row in case_rows:
        p = row["threads"]
        s = row.get("speedup", float("nan"))
        if p <= 1 or not math.isfinite(s) or s <= 0.0:
            continue
        denom = 1.0 - 1.0 / p
        if abs(denom) < 1.0e-14:
            continue
        f = (1.0 / s - 1.0 / p) / denom
        if math.isfinite(f):
            estimates.append(min(1.0, max(0.0, f)))

    if not estimates:
        return float("nan")

    return sum(estimates) / len(estimates)


def main() -> None:
    input_csv = sys.argv[1] if len(sys.argv) > 1 else "performance_raw.csv"

    rows = []
    with open(input_csv, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            row["threads"] = int(row["threads"])
            row["total_runtime_s"] = to_float(row["total_runtime_s"])
            rows.append(row)

    by_case = defaultdict(list)
    for row in rows:
        by_case[row["case"]].append(row)

    output_rows = []
    by_case_output = defaultdict(list)
    for case, case_rows in by_case.items():
        case_rows.sort(key=lambda x: x["threads"])
        base = next((r["total_runtime_s"] for r in case_rows if r["threads"] == 1), None)
        if base is None or base <= 0.0:
            continue

        for row in case_rows:
            runtime = row["total_runtime_s"]
            threads = row["threads"]
            speedup = base / runtime if runtime > 0.0 else float("nan")
            efficiency = speedup / threads if threads > 0 else float("nan")

            output_rows.append(
                {
                    "case": case,
                    "threads": threads,
                    "runtime_s": runtime,
                    "speedup": speedup,
                    "efficiency": efficiency,
                }
            )

            by_case_output[case].append(
                {
                    "case": case,
                    "threads": threads,
                    "runtime_s": runtime,
                    "speedup": speedup,
                    "efficiency": efficiency,
                }
            )

    with open("performance_speedup.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["case", "threads", "runtime_s", "speedup", "efficiency"])
        writer.writeheader()
        writer.writerows(output_rows)

    amdahl_rows = []
    for case in sorted(by_case_output.keys()):
        f_serial = estimate_serial_fraction(by_case_output[case])
        amdahl_rows.append({"case": case, "serial_fraction": f_serial})

    with open("performance_amdahl_fit.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["case", "serial_fraction"])
        writer.writeheader()
        writer.writerows(amdahl_rows)

    # Plot speedup
    plt.figure(figsize=(8, 5), dpi=150)
    for case in sorted({r["case"] for r in output_rows}):
        case_data = [r for r in output_rows if r["case"] == case]
        case_data.sort(key=lambda x: x["threads"])
        x = [r["threads"] for r in case_data]
        y = [r["speedup"] for r in case_data]
        plt.plot(x, y, marker="o", linewidth=2, label=case)

    all_threads = sorted({r["threads"] for r in output_rows})
    if all_threads:
        plt.plot(all_threads, all_threads, "k--", linewidth=1.5, label="Ideal")

    plt.xlabel("Threads")
    plt.ylabel("Speedup")
    plt.title("OpenMP Speedup")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("speedup_plot.png", bbox_inches="tight")
    plt.close()

    # Plot efficiency
    plt.figure(figsize=(8, 5), dpi=150)
    for case in sorted({r["case"] for r in output_rows}):
        case_data = [r for r in output_rows if r["case"] == case]
        case_data.sort(key=lambda x: x["threads"])
        x = [r["threads"] for r in case_data]
        y = [100.0 * r["efficiency"] for r in case_data]
        plt.plot(x, y, marker="o", linewidth=2, label=case)

    if all_threads:
        plt.plot(all_threads, [100.0 for _ in all_threads], "k--", linewidth=1.5, label="Ideal (100%)")

    plt.xlabel("Threads")
    plt.ylabel("Efficiency [%]")
    plt.title("OpenMP Parallel Efficiency")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("efficiency_plot.png", bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    main()
