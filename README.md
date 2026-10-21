# LPGDS sensor array

**Engineered L-PGDS–DPH fluorescence sensor arrays for fatty acids and PFAS**

This repository contains reproducible Python workflows for six engineered L-PGDS receptor channels: WT, LA, LB, LC, LD, and LE.

## Research modules

| Directory | Analytes | Workflow |
| --- | --- | --- |
| [Fatty_acids](Fatty_acids/) | Fatty-acid chain lengths, unsaturation, concentration gradients and mixtures | LDA, response profiles, training and held-out predictions |
| [PFAS](PFAS/) | PFOA, PFNA, PFDA, PFHxS, PFOS; PFOA/PFOS mixtures and concentration gradients | LDA, response profiles, held-out confusion matrices |
| [DPH_LPGDS_WT_MD](DPH_LPGDS_WT_MD/) | DPH and wild-type L-PGDS | WT–DPH MD simulation data, figures, restorable raw trajectories, and English SI methods |
| [DPH_FA_LPGDS_WT_MD](DPH_FA_LPGDS_WT_MD/) | C18:2 added near WT–DPH | FA_R2 MD data, English figures, restorable raw trajectories, and English SI methods |

## Reproduce

Fatty acids:

```bash
pip install -r Fatty_acids/requirements.txt
python Fatty_acids/analysis_pipeline.py --input fatty_acids.xlsx --output-dir Fatty_acids/analysis_output --seed 20260921
```

PFAS:

```bash
pip install -r PFAS/requirements.txt
python PFAS/analyze_pfas.py --input pfas.xlsx --output PFAS/results --seed 20261017
```

Use the actual cleaned workbooks as inputs; these are not part of this code repository. For sheet names and workflow details, consult each module's README.

## Shared analysis principles

- Six fluorescence-response channels: WT, LA, LB, LC, LD and LE.
- Ten replicates per analyte class, with six randomly assigned to training and four held out for prediction.
- Standardization and LDA are fitted on training observations only.
- The 95% ellipses are descriptive covariance regions for the training scores.
- Classification accuracy from a single split does not constitute independent external validation.
- Bar-chart errors come from the supplied experimental summary sheets; confirm whether their uncertainties are SD or SEM.

## Scope

The repository documents **both fatty-acid and PFAS fingerprint discrimination** on the L-PGDS–DPH platform. The two datasets have separate scripts and dependencies so either study can be reproduced independently.

## DPH–L-PGDS (WT) molecular dynamics

The [MD module](DPH_LPGDS_WT_MD/) documents one WT–DPH molecular dynamics trajectory, including simulation inputs, analysis, figures, and [English SI methods](DPH_LPGDS_WT_MD/methods/SI_methods_EN.md). Original large trajectories and XML files are losslessly stored as gzip parts. Run `python3 DPH_LPGDS_WT_MD/scripts/restore_data.py` to restore the archived files and `python3 DPH_LPGDS_WT_MD/scripts/verify_data.py` to check the deposited SHA-256 hashes. This single trajectory provides descriptive structural observations, not an independently replicated binding assessment or an affinity measurement.


## C18:2 near WT–DPH molecular dynamics

The [FA_R2 MD module](DPH_FA_LPGDS_WT_MD/) documents C18:2 behavior after addition near the selected WT–DPH endpoint, including complete restorable trajectories, analysis, English figures and [SI methodology](DPH_FA_LPGDS_WT_MD/methods/SI_methods_EN.md). This selected 5 ns trajectory describes local contacts and geometric shielding; it does not establish cooperative affinity, enclosed-pocket residence or a fluorescence response.
