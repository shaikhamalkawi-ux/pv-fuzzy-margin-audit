# PV Fuzzy Margin R1 — Protocol Lock

## Project separation
R1 is a new fuzzy/possibilistic paper line. It does not modify the frozen upstream ROM bases, critical-mode anchors, or the 2.5% comparator guard.

## Primary data and frozen selection
- Data: Oahu Solar Measurement Grid, July 2010, 1-s GHI.
- Source DOI: `10.7799/1052451`.
- Ten stations: DH1, DH9, DH8, DH2, DH4, DH3, AP1, AP3, AP4, AP6.
- Twelve 60-s windows are selected outcome-blind: four mean-GHI bins crossed with three empirical spatial-CV tertiles.
- All 720 observed one-second profiles in the frozen windows are used.

## Primary fuzzy definition
For station n inside a frozen 60-s window:
`(a,b,c,d) = (min, Q10, Q90, max)`.

Joint profile membership is the minimum of the ten station memberships, evaluated only at actually observed simultaneous profiles.

Primary alpha levels: `0, 0.25, 0.50, 0.75, 1.00`.

The alpha=0 set is the full finite observed universe of 60 simultaneous profiles in that minute.

## Physical calculation
Each measured ten-position vector is used directly and repeated across the 40 strings. No power matching, reshaping, retraining, re-anchoring, or guard tuning is allowed.

Models:
- exact FOM coherent/transverse oracle for repeated profiles;
- Aggregate15;
- SPPK35;
- generic SPPK45;
- anchored 45-state ROM.

The predeclared Xg scan is 0.15–1.8 pu. Only unique stable-to-unstable crossings are admitted to boundary-error summaries; nonunique cases are HOLD.

## Primary outputs
- fuzzy critical-SCR alpha-cut intervals;
- normalized alpha-cut Hausdorff error;
- interval overlap;
- paired false-stable/false-unstable profile counts and false-stable union widths;
- one-sided required guard `g_req=max(0,SCR_F/SCR_R-1)`;
- possibility/necessity that the frozen 2.5% comparator is exceeded.

## Claim boundary
Allowed:
- measurement-anchored fuzzy/possibilistic input uncertainty;
- fuzzy FOM/ROM physical-boundary fidelity;
- graded decision-envelope and guard-scope diagnostics;
- comparison with crisp, support-only, and empirical-percentile summaries.

Not allowed:
- fuzzy-controller claims;
- field/hardware validation of benchmark critical SCR;
- probability interpretation of membership or possibility;
- universal correctness of one membership function;
- post-hoc recalibration of the inherited 2.5% guard.
