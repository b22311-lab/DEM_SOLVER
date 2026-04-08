#!/usr/bin/env python3
import csv
import math
import sys

import matplotlib.pyplot as plt


def to_float(v: str) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def main() -> None:
    input_csv = sys.argv[1] if len(sys.argv) > 1 else "weak_scaling_raw.csv"

    rows = []
    with open(input_csv, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            row["threads"] = int(row["threads"])
            row["nparticles"] = int(row["nparticles"])
            row["total_runtime_s"] = to_float(row["total_runtime_s"])
            row["particle_contacts_s"] = to_float(row["particle_contacts_s"])
            row["contact_candidates_total"] = to_float(row["contact_candidates_total"])
            row["contacts_detected_total"] = to_float(row["contacts_detected_total"])
            rows.append(row)

    rows.sort(key=lambda r: r["threads"])
    if not rows:
        raise RuntimeError("No weak-scaling rows found")

    t1 = next((r["total_runtime_s"] for r in rows if r["threads"] == 1), float("nan"))
    c1 = next((r["particle_contacts_s"] for r in rows if r["threads"] == 1), float("nan"))

    summary = []
    for r in rows:
        tp = r["total_runtime_s"]
        cp = r["particle_contacts_s"]
        weak_eff_total = t1 / tp if tp > 0 and math.isfinite(t1) else float("nan")
        weak_eff_contact = c1 / cp if cp > 0 and math.isfinite(c1) else float("nan")

        summary.append(
            {
                "threads": r["threads"],
                "nparticles": r["nparticles"],
                "runtime_s": tp,
                "normalized_runtime": tp / t1 if t1 > 0 and math.isfinite(t1) else float("nan"),
                "weak_efficiency_total": weak_eff_total,
                "contact_runtime_s": cp,
                "weak_efficiency_contact": weak_eff_contact,
                "candidates_total": r["contact_candidates_total"],
                "contacts_total": r["contacts_detected_total"],
            }
        )

    with open("weak_scaling_summary.csv", "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "threads",
                "nparticles",
                "runtime_s",
                "normalized_runtime",
                "weak_efficiency_total",
                "contact_runtime_s",
                "weak_efficiency_contact",
                "candidates_total",
                "contacts_total",
            ],
        )
        writer.writeheader()
        writer.writerows(summary)

    x = [r["threads"] for r in summary]

    plt.figure(figsize=(7.2, 4.8), dpi=150)
    plt.plot(x, [r["runtime_s"] for r in summary], "o-", linewidth=2, label="Total runtime")
    plt.plot(x, [r["contact_runtime_s"] for r in summary], "s-", linewidth=2, label="Contact runtime")
    plt.plot(x, [summary[0]["runtime_s"] for _ in x], "k--", linewidth=1.2, label="Ideal weak scaling")
    plt.xlabel("Threads")
    plt.ylabel("Runtime [s]")
    plt.title("Weak Scaling Runtime")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("weak_scaling_runtime.png", bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(7.2, 4.8), dpi=150)
    plt.plot(x, [100.0 * r["weak_efficiency_total"] for r in summary], "o-", linewidth=2, label="Total runtime")
    plt.plot(x, [100.0 * r["weak_efficiency_contact"] for r in summary], "s-", linewidth=2, label="Contact runtime")
    plt.plot(x, [100.0 for _ in x], "k--", linewidth=1.2, label="Ideal")
    plt.xlabel("Threads")
    plt.ylabel("Weak-scaling efficiency [%]")
    plt.title("Weak-Scaling Efficiency")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("weak_scaling_efficiency.png", bbox_inches="tight")
    plt.close()

    with open("weak_scaling_findings.tex", "w") as f:
        f.write("\\subsection{Weak-Scaling Study}\n")
        f.write("Weak scaling was evaluated by keeping particles per thread constant ($N/p=500$) while increasing total particles and domain size with thread count.\n")
        f.write("\\begin{table}[!t]\n")
        f.write("\\centering\n")
        f.write("\\caption{Weak-scaling results ($N \\propto p$)}\n")
        f.write("\\label{tab:weak-scaling}\n")
        f.write("\\begin{tabular}{rcccc}\n")
        f.write("\\toprule\n")
        f.write("Threads & $N$ & Runtime [s] & Norm. runtime & Weak eff. [\\%] \\\\" + "\n")
        f.write("\\midrule\n")
        for r in summary:
            f.write(
                f"{r['threads']} & {r['nparticles']} & {r['runtime_s']:.4f} & {r['normalized_runtime']:.3f} & {100.0*r['weak_efficiency_total']:.1f} \\\\" + "\n"
            )
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")


if __name__ == "__main__":
    main()
