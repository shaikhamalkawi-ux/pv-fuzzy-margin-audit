# PV Fuzzy Margin Audit

Public development repository for a measurement-anchored fuzzy audit of weak-grid stability margins in reduced-order photovoltaic models.

## Current scope

This repository contains **reproducibility material only** during the blind-review stage. The manuscript PDF, author metadata, and citation metadata are intentionally not mirrored here and the repository is **not cited in the submitted paper at this stage**.

The analysis asks whether a reduced-order PV model that appears conservative at a crisp minute-mean operating point also preserves the uncertainty envelope of the model-computed critical short-circuit-ratio (SCR) boundary when measured irradiance varies within that minute.

### Frozen headline results

- 12 outcome-blind one-minute OSMG windows.
- 720 simultaneous measured ten-station irradiance profiles.
- Maximum one-sided correction for the 12 crisp minute means: **0.869%**.
- Maximum correction across the measured within-window support: **5.572%**.
- **9/720** measured profiles require more than the frozen **2.5%** comparator guard.
- All nine exceedances occur in window **W11**.

These are model-computed stability boundaries conditioned on measured irradiance inputs; they are **not field measurements of a stability boundary**, and 9/720 is **not a probability-of-failure estimate**.

## Repository layout

- `code/` — analysis and verification scripts.
- `protocol/` — frozen window-selection and analysis protocol.
- `results/` — compact derived result tables used for verification.
- `THIRD_PARTY_DATA.md` — source and redistribution boundary for OSMG data.
- `BLIND_REVIEW_NOTICE.md` — current review-stage publication boundary.

## Data source

Oahu Solar Measurement Grid (OSMG), 1-second solar irradiance measurements, DOI: **10.7799/1052451**.

Raw OSMG measurements are third-party data and are not redistributed in this repository. The scripts operate on a locally acquired source archive.

## Reproducibility boundary

The public repository intentionally omits the manuscript and author-facing submission files during blind review. A post-acceptance release can add the archival paper citation, CITATION.cff, DOI-tagged release, and the remaining public-safe upstream model package.

## Important interpretation

Fuzzy membership and possibility are used as graded representations of measured within-window irradiance uncertainty. They are **not probability distributions**. The frozen 2.5% guard is evaluated as a comparator and is not re-tuned from these fuzzy results.

