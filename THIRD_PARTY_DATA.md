# Third-party data boundary

## Oahu Solar Measurement Grid (OSMG)

The analysis uses measured irradiance from the Oahu Solar Measurement Grid (OSMG).

- Dataset: **Oahu Solar Measurement Grid (1-Year Archive): 1-Second Solar Irradiance; Oahu, Hawaii**
- DOI: **10.7799/1052451**
- Source organization: NREL/NLR data archive

Raw OSMG measurements are **not redistributed** in this repository.

The repository may contain:
- frozen station/window mappings;
- derived profile identifiers;
- fuzzy-input parameters;
- model-computed boundary summaries;
- verification scripts.

To reproduce the analysis, acquire the OSMG source data from the official archive and run the preprocessing script locally.

The measured data provide irradiance inputs only. The reported critical-SCR/stability boundaries remain model-computed quantities and are not field-measured stability labels.
