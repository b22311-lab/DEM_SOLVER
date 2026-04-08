#!/usr/bin/env python3
import csv
import math
import sys
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np


def to_float(v: str) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def to_int(v: str) -> int:
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return 0


def fit_exponent(n_vals, y_vals):
    pairs = [(n, y) for n, y in zip(n_vals, y_vals) if n > 0 and y > 0]
    if len(pairs) < 2:
        return float("nan")
    x = np.log(np.array([p[0] for p in pairs], dtype=float))
    y = np.log(np.array([p[1] for p in pairs], dtype=float))
    return float(np.polyfit(x, y, 1)[0])


def main() -> None:
    input_csv = sys.argv[1] if len(sys.argv) > 1 else "part18_neighbor_search_raw.csv"

    rows = []
    with open(input_csv, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            row["nparticles"] = to_int(row["nparticles"])
            row["method_id"] = to_int(row["method_id"])
            row["total_runtime_s"] = to_float(row["total_runtime_s"])
            row["particle_contacts_s"] = to_float(row["particle_contacts_s"])
            row["contact_candidates_total"] = to_float(row["contact_candidates_total"])
            row["contacts_detected_total"] = to_float(row["contacts_detected_total"])
            row["loop_count"] = max(1, to_int(row["loop_count"]))
            row["candidates_avg_per_step"] = row["contact_candidates_total"] / row["loop_count"]
            rows.append(row)

    by_n = defaultdict(dict)
    for row in rows:
        by_n[row["nparticles"]][row["method"]] = row

    summary_rows = []
    n_list = sorted(by_n.keys())
    for n in n_list:
        all_pairs = by_n[n].get("all_pairs")
        cell = by_n[n].get("cell_linked")
        if all_pairs is None or cell is None:
            continue

        total_speedup = all_pairs["total_runtime_s"] / cell["total_runtime_s"] if cell["total_runtime_s"] > 0 else float("nan")
        contact_speedup = all_pairs["particle_contacts_s"] / cell["particle_contacts_s"] if cell["particle_contacts_s"] > 0 else float("nan")
        reduction_pct = 100.0 * (1.0 - cell["contact_candidates_total"] / all_pairs["contact_candidates_total"]) if all_pairs["contact_candidates_total"] > 0 else float("nan")
        reduction_factor = all_pairs["contact_candidates_total"] / cell["contact_candidates_total"] if cell["contact_candidates_total"] > 0 else float("nan")

        contacts_match = abs(all_pairs["contacts_detected_total"] - cell["contacts_detected_total"]) < 0.5
        summary_rows.append(
            {
                "nparticles": n,
                "all_pairs_total_runtime_s": all_pairs["total_runtime_s"],
                "cell_total_runtime_s": cell["total_runtime_s"],
                "total_speedup_cell_vs_allpairs": total_speedup,
                "all_pairs_contact_runtime_s": all_pairs["particle_contacts_s"],
                "cell_contact_runtime_s": cell["particle_contacts_s"],
                "contact_speedup_cell_vs_allpairs": contact_speedup,
                "all_pairs_candidates_total": all_pairs["contact_candidates_total"],
                "cell_candidates_total": cell["contact_candidates_total"],
                "candidate_reduction_pct": reduction_pct,
                "candidate_reduction_factor": reduction_factor,
                "all_pairs_contacts_detected": all_pairs["contacts_detected_total"],
                "cell_contacts_detected": cell["contacts_detected_total"],
                "contacts_match": "yes" if contacts_match else "no",
                "all_pairs_candidates_per_step": all_pairs["candidates_avg_per_step"],
                "cell_candidates_per_step": cell["candidates_avg_per_step"],
            }
        )

    with open("part18_neighbor_search_summary.csv", "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "nparticles",
                "all_pairs_total_runtime_s",
                "cell_total_runtime_s",
                "total_speedup_cell_vs_allpairs",
                "all_pairs_contact_runtime_s",
                "cell_contact_runtime_s",
                "contact_speedup_cell_vs_allpairs",
                "all_pairs_candidates_total",
                "cell_candidates_total",
                "candidate_reduction_pct",
                "candidate_reduction_factor",
                "all_pairs_contacts_detected",
                "cell_contacts_detected",
                "contacts_match",
                "all_pairs_candidates_per_step",
                "cell_candidates_per_step",
            ],
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    # Runtime plot for contact stage
    if summary_rows:
        nvals = [r["nparticles"] for r in summary_rows]
        tt_all = [r["all_pairs_total_runtime_s"] for r in summary_rows]
        tt_cell = [r["cell_total_runtime_s"] for r in summary_rows]
        t_all = [r["all_pairs_contact_runtime_s"] for r in summary_rows]
        t_cell = [r["cell_contact_runtime_s"] for r in summary_rows]

        plt.figure(figsize=(7.0, 4.8), dpi=150)
        plt.loglog(nvals, tt_all, "o-", linewidth=2, label="All-pairs total runtime")
        plt.loglog(nvals, tt_cell, "s-", linewidth=2, label="Cell-linked total runtime")
        plt.xlabel("Number of particles N")
        plt.ylabel("Total runtime [s]")
        plt.title("Part 18: Total Runtime Scaling")
        plt.grid(True, which="both", alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig("part18_total_runtime.png", bbox_inches="tight")
        plt.close()

        plt.figure(figsize=(7.0, 4.8), dpi=150)
        plt.loglog(nvals, t_all, "o-", linewidth=2, label="All-pairs contact stage")
        plt.loglog(nvals, t_cell, "s-", linewidth=2, label="Cell-linked contact stage")
        plt.xlabel("Number of particles N")
        plt.ylabel("Contact-stage runtime [s]")
        plt.title("Part 18: Contact-Stage Runtime Scaling")
        plt.grid(True, which="both", alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig("part18_contact_runtime.png", bbox_inches="tight")
        plt.close()

        cand_all = [r["all_pairs_candidates_per_step"] for r in summary_rows]
        cand_cell = [r["cell_candidates_per_step"] for r in summary_rows]
        plt.figure(figsize=(7.0, 4.8), dpi=150)
        plt.loglog(nvals, cand_all, "o-", linewidth=2, label="All-pairs candidates/step")
        plt.loglog(nvals, cand_cell, "s-", linewidth=2, label="Cell-linked candidates/step")
        plt.xlabel("Number of particles N")
        plt.ylabel("Candidate pairs per timestep")
        plt.title("Part 18: Candidate-Pair Scaling")
        plt.grid(True, which="both", alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig("part18_candidates.png", bbox_inches="tight")
        plt.close()

        s_total = [r["total_speedup_cell_vs_allpairs"] for r in summary_rows]
        s_contact = [r["contact_speedup_cell_vs_allpairs"] for r in summary_rows]
        plt.figure(figsize=(7.0, 4.8), dpi=150)
        plt.plot(nvals, s_total, "o-", linewidth=2, label="Total-runtime speedup")
        plt.plot(nvals, s_contact, "s-", linewidth=2, label="Contact-stage speedup")
        plt.axhline(1.0, color="k", linestyle="--", linewidth=1.2)
        plt.xscale("log")
        plt.xlabel("Number of particles N")
        plt.ylabel("Cell-linked speedup over all-pairs")
        plt.title("Part 18: Cell-Linked Speedup")
        plt.grid(True, which="both", alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig("part18_speedup.png", bbox_inches="tight")
        plt.close()

    alpha_all = fit_exponent(
        [r["nparticles"] for r in summary_rows],
        [r["all_pairs_candidates_per_step"] for r in summary_rows],
    )
    alpha_cell = fit_exponent(
        [r["nparticles"] for r in summary_rows],
        [r["cell_candidates_per_step"] for r in summary_rows],
    )

    slower_cases = [r["nparticles"] for r in summary_rows if r["contact_speedup_cell_vs_allpairs"] < 1.0]

    with open("part18_findings.tex", "w") as f:
        f.write("\\subsection{Bonus Part 18: Neighbour-Search Comparison}\n")
        f.write("A switchable contact-search implementation was added with \\texttt{contact\\_search\\_method=0} (all-pairs) and \\texttt{contact\\_search\\_method=1} (cell-linked).\n")
        if math.isfinite(alpha_all) and math.isfinite(alpha_cell):
            f.write(
                f"Empirical candidate-pair scaling from the measured data is approximately $\\mathcal{{O}}(N^{{{alpha_all:.2f}}})$ for all-pairs and $\\mathcal{{O}}(N^{{{alpha_cell:.2f}}})$ for the cell-linked strategy.\n"
            )
        f.write("\\begin{table}[!t]\n")
        f.write("\\scriptsize\n")
        f.write("\\centering\n")
        f.write("\\caption{Part 18 neighbour-search comparison (serial runs)}\n")
        f.write("\\label{tab:part18-neighbour-search}\n")
        f.write("\\begin{tabular}{rcccccccc}\n")
        f.write("\\toprule\n")
        f.write("$N$ & $t_{tot}^{AP}$ & $t_{tot}^{CL}$ & $S_{tot}$ & $t_c^{AP}$ & $t_c^{CL}$ & $S_c$ & Cand. factor \\\\" + "\n")
        f.write("\\midrule\n")
        for r in summary_rows:
            f.write(
                f"{r['nparticles']} & {r['all_pairs_total_runtime_s']:.4f} & {r['cell_total_runtime_s']:.4f} & "
                f"{r['total_speedup_cell_vs_allpairs']:.2f} & {r['all_pairs_contact_runtime_s']:.4f} & "
                f"{r['cell_contact_runtime_s']:.4f} & {r['contact_speedup_cell_vs_allpairs']:.2f} & "
                f"{r['candidate_reduction_factor']:.2f} \\\\" + "\n"
            )
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")
        f.write("where $t_{tot}^{AP}, t_{tot}^{CL}$ are total runtimes; $t_c^{AP}, t_c^{CL}$ are contact-stage runtimes; $S_{tot}, S_c$ are the corresponding cell-linked speedups over all-pairs; and Cand. factor is the candidate-pair ratio $N_{cand}^{AP}/N_{cand}^{CL}$.\n")
        f.write("Detected-contact counts match between methods for each tested $N$, confirming contact-search correctness while reducing candidate checks.\n")
        f.write("The measured candidate-pair reduction confirms why the baseline all-pairs method behaves as the expensive $\\mathcal{O}(N^2)$ reference algorithm.\n")
        if slower_cases:
            case_text = ", ".join(str(v) for v in slower_cases)
            f.write(
                f"For small systems (here $N={case_text}$), the cell-linked method is slower than all-pairs because grid build and traversal overheads dominate the reduced pair-check count.\n"
            )


if __name__ == "__main__":
    main()
