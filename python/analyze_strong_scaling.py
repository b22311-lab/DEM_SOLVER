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
    input_csv = sys.argv[1] if len(sys.argv) > 1 else "strong_scaling_raw.csv"

    rows = []
    with open(input_csv, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            row["threads"] = int(row["threads"])
            row["nparticles"] = int(row["nparticles"])
            row["total_runtime_s"] = to_float(row["total_runtime_s"])
            row["particle_contacts_s"] = to_float(row["particle_contacts_s"])
            rows.append(row)

    rows.sort(key=lambda r: r["threads"])
    if not rows:
        raise RuntimeError("No strong-scaling rows found")

    t1 = next((r["total_runtime_s"] for r in rows if r["threads"] == 1), float("nan"))
    c1 = next((r["particle_contacts_s"] for r in rows if r["threads"] == 1), float("nan"))

    summary = []
    for r in rows:
        p = r["threads"]
        total = r["total_runtime_s"]
        contact = r["particle_contacts_s"]
        s_total = t1 / total if total > 0 and math.isfinite(t1) else float("nan")
        e_total = s_total / p if p > 0 and math.isfinite(s_total) else float("nan")
        s_contact = c1 / contact if contact > 0 and math.isfinite(c1) else float("nan")
        e_contact = s_contact / p if p > 0 and math.isfinite(s_contact) else float("nan")

        summary.append(
            {
                "case": r["case"],
                "nparticles": r["nparticles"],
                "threads": p,
                "runtime_s": total,
                "speedup_total": s_total,
                "efficiency_total": e_total,
                "contact_runtime_s": contact,
                "speedup_contact": s_contact,
                "efficiency_contact": e_contact,
                "candidates_total": r["contact_candidates_total"],
                "contacts_total": r["contacts_detected_total"],
            }
        )

    with open("strong_scaling_summary.csv", "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "case",
                "nparticles",
                "threads",
                "runtime_s",
                "speedup_total",
                "efficiency_total",
                "contact_runtime_s",
                "speedup_contact",
                "efficiency_contact",
                "candidates_total",
                "contacts_total",
            ],
        )
        writer.writeheader()
        writer.writerows(summary)

    x = [r["threads"] for r in summary]
    y_total = [r["speedup_total"] for r in summary]
    y_contact = [r["speedup_contact"] for r in summary]

    plt.figure(figsize=(7.2, 4.8), dpi=150)
    plt.plot(x, y_total, "o-", linewidth=2, label="Total runtime")
    plt.plot(x, y_contact, "s-", linewidth=2, label="Contact stage")
    plt.plot(x, x, "k--", linewidth=1.2, label="Ideal")
    plt.xlabel("Threads")
    plt.ylabel("Speedup")
    plt.title("Strong Scaling")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("strong_scaling_speedup.png", bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(7.2, 4.8), dpi=150)
    plt.plot(x, [100.0 * r["efficiency_total"] for r in summary], "o-", linewidth=2, label="Total runtime")
    plt.plot(x, [100.0 * r["efficiency_contact"] for r in summary], "s-", linewidth=2, label="Contact stage")
    plt.plot(x, [100.0 for _ in x], "k--", linewidth=1.2, label="Ideal")
    plt.xlabel("Threads")
    plt.ylabel("Efficiency [%]")
    plt.title("Strong-Scaling Efficiency")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("strong_scaling_efficiency.png", bbox_inches="tight")
    plt.close()

    case = summary[0]["case"]
    n = summary[0]["nparticles"]
    case_tex = case.replace("_", "\\_")
    with open("strong_scaling_findings.tex", "w") as f:
        f.write("\\subsection{Strong-Scaling Study}\n")
        f.write(f"Strong scaling was evaluated by fixing the problem size (case \\texttt{{{case_tex}}}, $N={n}$) and varying OpenMP threads.\n")
        f.write("\\begin{table}[!t]\n")
        f.write("\\centering\n")
        f.write("\\caption{Strong-scaling results for fixed particle count}\n")
        f.write("\\label{tab:strong-scaling}\n")
        f.write("\\begin{tabular}{rcccc}\n")
        f.write("\\toprule\n")
        f.write("Threads & Runtime [s] & $S_p$ & $E_p$ [\\%] & Contact speedup \\\\" + "\n")
        f.write("\\midrule\n")
        for r in summary:
            f.write(
                f"{r['threads']} & {r['runtime_s']:.4f} & {r['speedup_total']:.3f} & {100.0*r['efficiency_total']:.1f} & {r['speedup_contact']:.3f} \\\\" + "\n"
            )
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")


if __name__ == "__main__":
    main()
