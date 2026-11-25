# DME-TIR-Interval

Reproducible analysis code for **within-eye sensitivity of first-year visual-acuity time in range to assessment schedules** using DRCR Protocol T and AC public trial datasets.

The formal analysis compares 4-, 8- and 12-week assessment schedules imposed on existing measurements. The locked strict samples are Protocol T **401 eyes** (primary two-arm treatment model: **274 eyes**) and Protocol AC **157 eyes / 139 participants**. The unique primary estimand is the mean absolute within-eye difference between the 8- and 4-week schedules. Actual examination days, exact piecewise-linear integration at VA ≥69 letters, a shared 0–364-day window and 5,000 participant bootstrap resamples are retained. No participant-level data are distributed.

## Setup and synthetic checks

Use Python 3.12 or later (the release is checked with the original project interpreter). NumPy is pinned to the original analysis version.

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python reproduce.py --synthetic-only
```

Synthetic fixtures are invented and labeled SYNTHETIC ONLY. Tests cover threshold equality, hidden intermediate dips and peaks, actual irregular dates, annual clipping/tail restrictions, calendar alignment after missing visits, and participant-level bilateral resampling.

## Full reproduction with separately obtained data

Read [data access, source attribution and archive hashes](docs/DATA_ACCESS.md), then place both official ZIP files in the documented `data/raw/` directories. Do not commit them.

```bash
.venv/bin/python reproduce.py --out work/run_001
```

The output directory must be new and inside this checkout. The runner copies code/configuration into an isolated local run, links the two verified source packages, regenerates qualification/masks and executes formal analysis, independent numerical checks, frozen secondary/sensitivity analyses, and the post hoc clinical/selection/common-support extensions. All 5,000 bootstrap replicates are retained. It then generates two summary figures as PNG/SVG and compares generated TSV digests to frozen final-table fingerprints. Logs and the machine-readable verification report are inside the run. Keep the entire run directory private: it contains individual-level intermediates.

For only the formal analysis, use `--formal-only`. Repeating a run requires a new output path. No manuscript writing, submission or external upload is triggered by the runner.

## Code map

- `phase0/`: frozen rules, eligibility/calendar logic and synthetic checks.
- `phase0_real_2026-11-24/`: pre-access population addendum and baseline implementation note.
- `phase1/`: original formal analysis, treatment models, holdout validation and sensitivities.
- `phase2_ajo/`, `phase2_ajo_revision/`: separately labeled post hoc endpoint, clinical, selection-weighting and common-baseline-support analyses. Historical directory names are retained for imports; repository naming is journal-neutral.
- `docs/`: raw-field mappings, source attribution and release engineering notes.
- `provenance/`: whole-file digests, without individual data.
- `plot_results.py`: figures generated directly from local final tables.

AC selection-weighting results are diagnostic only because the prespecified qualification gate failed. Original randomized strategy remains fixed after switching. All interpretation is limited to thinning records along existing treatment paths.

## License and citation

Original code is released under the [MIT license](LICENSE). Dataset terms remain with the DRCR/Jaeb source (see [data access](docs/DATA_ACCESS.md)). Archive creator/author is the **DME-TIR-Interval** project. Zenodo metadata are controlled by `.zenodo.json`; citation metadata are in `CITATION.cff`.

The release validation record is in `provenance/release_validation.json`; all 60 regenerated final TSVs were byte-identical to the existing frozen results. GitHub Actions checks synthetic cases without downloading trial data.
