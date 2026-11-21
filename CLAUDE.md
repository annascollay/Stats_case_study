# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A group coursework project (Statistical Case Studies 2026/27, University of Edinburgh; authors Anna Scollay, Lewis Seath, Zacharias Onyejiaka) analysing 15-minute river discharge for the **River Don at Parkhill** (SEPA). It is an R / R Markdown project, not a software package: there is no build system, linter or test suite. The deliverable is a knitted PDF report.

- `Project 1 Part 1/` — finished Part 1 (extreme flow peaks, threshold u = 65 m³/s, time-to-peak). Reference only.
- `Project 1 Part 2/` — **current work**: representing and simulating the *shape* of flood hydrographs for flood inundation modelling, compared against SEPA's current approach.
- `02. Rider_RMDs/` folders and root-level `SCS_P1P2_Lewis.R` are individual members' scratch work; the `01. Master_RMD/` folder is the canonical code + report.

## Running / knitting (R 4.6.1, Windows)

Always run from `Project 1 Part 2/01. Master_RMD/` — every path (data file, `cache/`, logo, `.bib`) is relative to that folder.

Knit the report (works from any directory; PDF is written next to the Rmd):

```bash
"/c/Program Files/R/R-4.6.1/bin/x64/Rscript" -e "setwd('C:/Users/ZachO/OneDrive - University of Edinburgh/01. Dissertations/Statistical Case Studies/SCS Projects/Stats_case_study/Project 1 Part 2/01. Master_RMD'); rmarkdown::render('SCS_P1P2_REPORT.Rmd')"
```

- Do not pass `output_format = 'pdf_document'`: it overrides the YAML's `bookdown::pdf_document2` and every `\@ref()` becomes "??". Don't pass `output_dir` either.
- Tables that are cross-referenced must be built with `kable(..., caption = ...)` directly inside their Rmd chunk (not as a pre-built object from the script), otherwise references render as "Table ??". Keep the table *data* in the script and the `kable()` call in the Rmd.

Run the analysis alone: `Rscript master_code_part2.R` from the same folder.

Packages: tidyverse, lubridate, knitr, kableExtra, patchwork, wavelets, splines, rlang, bookdown.

Data: `river_don_parkhill_flow.rds` (gitignored) is loaded if present in the working directory; otherwise the script downloads 2001–2026 from the SEPA KiWIS API and saves it there (slow).

## Architecture

**Split between script and report.** `master_code_part2.R` does *all* computation and builds every table, plot and caption object. `SCS_P1P2_REPORT.Rmd` just `source()`s it in the setup chunk and then prints those objects in chunks (`echo=FALSE`) interleaved with prose. To change a figure/table, edit the R script; to change wording/captions/layout, edit the Rmd. Inline `r ...` values in the prose reference variables from the script, so renaming a variable can silently break the report.

**Pipeline in `master_code_part2.R`** (numbered `# N.` section headers):
1–2. Load data; keep quality codes 50/100/200/254, set aside Suspect (150) and keep `river_data_raw` for checking it.
3. Independent peaks (Robson & Reed 1999 criteria, u = 65) — same logic as Part 1.
5. Checks of Unchecked/Suspect readings against linear interpolation inside event windows.
7. Event start/end trough-to-trough (`rise_start_row`, `fall_end_row`), one row per event in `events`.
8. Standardise each event → matrix **`Y`** (events × `n_grid` = 64 points): time rescaled to [0,1], peak realigned to `u_star` (mean peak position), flow divided by peak flow.
9. `reconstruct_hydrograph()` inverts the standardisation back to m³/s and hours.
10. Shared evaluation framework used by every method so they compare fairly: fitted curves are reconstructed to real units and scored by `flood_errors()` (RMSE, peak error, timing error, % volume error); `loocv()` for leave-one-out; fixed train/test split (`set.seed(2026)`); `plot_choice`, `train_test_choice`, `insample_loo_choice`, `cmp_row`.
11. SEPA baseline: one average (unaligned) shape scaled by peak flow.
12–16. Shape methods: polynomial, PCA, wavelets, gamma-type hydrograph function, B-splines. Per-event methods use train/test; PCA and wavelets (fitted across events) use leave-one-out.
18. Validation plots and the comparison of all methods (accuracy, parsimony, paired tests, synthesis).
19+. Uncertainty: distribution/dependence of shape parameters, simulation via multivariate Normal vs smoothed bootstrap, held-out bands, final simulated floods, sampling uncertainty via refits.

New methods should plug into the section 10 helpers rather than computing errors their own way.

## Caching

Expensive steps are wrapped in `cached(name, expr, key)`, which writes `cache/<name>_<hash>.rds`, the hash being of `Y` and `key`. The cache refreshes automatically when the data or `key` settings change, but **not when function code changes** — after editing a function used inside a `cached()` call, delete the matching `cache/<name>_*.rds` files (or set `use_cache <- FALSE`). The top-level `cache/` at the repo root is a stale copy from running with the wrong working directory; the live one is `Project 1 Part 2/01. Master_RMD/cache/`.

## Report conventions

- Methods and results are described in plain language aimed at SEPA (e.g. "extreme flow event", "standardised time"), consistent with existing captions.
- Cross-references use bookdown labels (`{#method-pca}`, `\@ref(...)`); citations come from `Project_2_references.bib` (APA CSL).
- Figures float with `H` placement; figure sizes are set per chunk.
