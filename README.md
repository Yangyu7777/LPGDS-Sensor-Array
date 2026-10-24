# L-PGDS–DPH Sensor Array

Computationally engineered L-PGDS receptor arrays for fluorescence fingerprinting of fatty acids (FAs) and per- and polyfluoroalkyl substances (PFAS).

## Repository contents

| Folder | Contents |
| --- | --- |
| [Fatty_acids/](Fatty_acids/) | FA chain-length, unsaturation, concentration, and binary-mixture analyses; fluorescence profiles, LDA and held-out classification |
| [PFAS/](PFAS/) | PFOA, PFNA, PFDA, PFHxS and PFOS analyses; concentration series, PFOA/PFOS mixtures, LDA and confusion matrices |
| [DPH_LPGDS_WT_MD/](DPH_LPGDS_WT_MD/) | WT L-PGDS–DPH simulation, trajectory archives, structural analysis, figures and SI methods |
| [DPH_FA_LPGDS_WT_MD/](DPH_FA_LPGDS_WT_MD/) | WT L-PGDS–DPH–C18:2 simulation, restored trajectories, shielding/hydration/contact analyses, figures and SI methods |

## Data analysis

All sensing analyses use six channels (WT, LA, LB, LC, LD and LE). Each analyte condition has ten measurements, split into six training and four held-out samples. Standardization and linear discriminant analysis (LDA) are fitted using training data only. The 95% ellipses summarize training-score distributions; single-split test accuracy is not external validation.

### Reproduce the analysis

```bash
pip install -r Fatty_acids/requirements.txt
python Fatty_acids/analysis_pipeline.py --input fatty_acids.xlsx --output-dir Fatty_acids/analysis_output --seed 20260921

pip install -r PFAS/requirements.txt
python PFAS/analyze_pfas.py --input pfas.xlsx --output PFAS/results --seed 20261017
```

Supply the cleaned experimental workbooks separately; they are not included in this repository. See the README in each analysis folder for required sheets and outputs. Check experimental summary sheets to establish whether plotted error bars denote SD or SEM.

## Molecular dynamics

**WT–DPH:** The [binary-system module](DPH_LPGDS_WT_MD/) contains one simulation trajectory, analysis files and [English SI methods](DPH_LPGDS_WT_MD/methods/SI_methods_EN.md). Restore compressed trajectory parts using `python3 DPH_LPGDS_WT_MD/scripts/restore_data.py` and verify checksums using `python3 DPH_LPGDS_WT_MD/scripts/verify_data.py`.

**WT–DPH–C18:2:** The [ternary-system module](DPH_FA_LPGDS_WT_MD/) uses **9Z,12Z-linoleate (CAS 60-33-3)**, replacing the earlier E/E isomer (CAS 506-21-8). It contains the selected 5 ns trajectory, available same-batch repeat summaries, contact/hydration/shielding analyses, figures and [English SI methods](DPH_FA_LPGDS_WT_MD/methods/SI_methods_EN.md). Earlier versions remain accessible through Git history.

These short MD trajectories provide descriptive structural observations; they do not independently establish binding affinity, cooperative binding, or the origin of fluorescence changes.
