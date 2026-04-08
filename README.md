# DEM Solver (Fortran + Python)

This repository contains a modular Discrete Element Method (DEM) solver for spherical particles with:

- gravity, particle-particle contact, and wall contact
- semi-implicit Euler time integration
- OpenMP parallelization for selected kernels
- reproducible verification, profiling, scaling, and report-generation scripts

## 1) Portable environment setup

Use the provided Conda specification to reproduce the Python and build-tool dependencies.

```bash
conda env create -f environment.yml
conda activate dem-solver
```

If your Conda environment name differs, set it when running scripts:

```bash
CONDA_ENV_NAME=my-env scripts/run_part19_2_study.sh
```

## 2) Build the solver

```bash
make -j
```

This generates `demo.exe`.
The Makefile applies OpenMP via `OMPFLAGS` so the build remains valid even if Conda predefines `FFLAGS`.

Cleanup commands:

```bash
make clean-build      # removes *.o, *.mod, demo.exe
make clean-artifacts  # removes generated png/csv/log/vtk and output folders
make clean            # runs both clean-build and clean-artifacts
```

## 3) One-command full artifact generation

To generate figures, CSVs, LaTeX snippets, and compile the paper PDF in one run:

```bash
./run_all_artifacts.sh
```

The script runs:

1. `scripts/run_all_tests.sh`
2. `scripts/run_part10_14.sh`
3. `scripts/run_part12_profile.sh`
4. `scripts/run_scaling_modes.sh`
5. `scripts/run_part18_study.sh`
6. `scripts/run_part19_2_study.sh`
7. `scripts/compile_latex_all.sh` (if available)

## 4) Key outputs

- Main paper: `ieee_dem_whitepaper.pdf`
- Verification plots: `test1_freefall_verification.png`, `test2_constant_velocity_verification.png`, `test3_particle_bounce_verification.png`
- Performance: `performance_raw.csv`, `performance_speedup.csv`, `speedup_plot.png`, `efficiency_plot.png`
- Scaling: `strong_scaling_*.png/csv`, `weak_scaling_*.png/csv`
- Part 18: `part18_neighbor_search_summary.csv`, `part18_*.png`, `part18_findings.tex`
- Part 19.2: `part19_2_timestep_sensitivity.csv`, `part19_2_richardson.csv`, `part19_2_critical_stability.csv`, `part19_2_*.png`, `part19_2_findings.tex`

## 5) Notes

- LaTeX compilation requires a local `pdflatex` installation with IEEEtran support.
- Scripts are path-portable and can be launched from any working directory.
- Most analysis scripts auto-select Conda (`dem-solver`) when available, otherwise they fall back to `python3`.

## 6) Benchmarked OpenMP Decisions

These implementation choices were benchmarked and are now documented here directly.

### 6.1 Why `if(nparticles > 32)` is used

Forcing OpenMP on very small particle counts was tested against the guarded implementation (`if(nparticles > 32)`) using `OMP_NUM_THREADS=8` and repeated runs.

- For `N=8,16,24,32`, forced parallel was slower by about `1.7x`, `3.4x`, `2.4x`, and `1.9x` (median runtime ratios).
- For larger `N` (`48+`), parallel begins to amortize overhead and becomes beneficial.

This supports keeping the `> 32` crossover guard for small-loop overhead control.

### 6.2 Why `force_private(:,:,tid)` is used

The all-pairs kernel was benchmarked against a direct OpenMP array-reduction variant (`reduction(+:force)`) at `OMP_NUM_THREADS=8`.

- At `N=1000`, array-reduction median total runtime was about `1.51x` slower.
- At `N=5000`, array-reduction median total runtime was about `1.10x` slower.
- At `N=10000`, array-reduction median total runtime was about `1.15x` slower.

In this setup, `force_private(:,:,tid)` provided better runtime behavior for contact-heavy runs, so it remains the default strategy.

## 7) Practical Benchmark Defaults

The benchmark generator scripts now apply the following defaults:

- Use `1, 6, 8` threads by default (avoid high thread counts like `12` unless intentionally testing them).
- Run repeated measurements and aggregate with **median**:
	- `REPEATS=5` by default (set `REPEATS=3` for faster checks).
- Enable affinity pinning by default for cleaner runs:
	- `OMP_PROC_BIND=true`
	- `OMP_PLACES=cores`
- Reduce timing contamination from file output by setting a very large `write_interval` in temporary timing configs.

This makes runtime comparisons target solver compute behavior rather than VTK I/O overhead.

Example override:

```bash
THREADS_LIST="1 6 8" REPEATS=5 OMP_PROC_BIND=true OMP_PLACES=cores scripts/run_part10_14.sh
```

For best reproducibility, benchmark when the machine is otherwise quiet.


# DEM_SOLVER
