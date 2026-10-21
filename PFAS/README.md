# PFAS fluorescence fingerprint analysis

Reproducible LDA and bar-chart analysis for the six-channel L-PGDS–DPH sensor array (WT, LA, LB, LC, LD, LE).

## Experiments
- Five PFAS: PFOA, PFNA, PFDA, PFHxS, PFOS
- PFOA/PFOS mixtures: 100:0, 75:25, 50:50, 25:75, 0:100
- PFOA concentration: 5, 7.5, 10, 15, 20 µM
- PFOS concentration: 5, 7.5, 10, 15, 20 µM

## Reproduce
Install Python 3.10+ and run:

```bash
pip install -r requirements.txt
python analyze_pfas.py --input "20261017-ALL LDA and Bar chart.xlsx" --output results --seed 20261017
```

The Excel input workbook is **not committed** here. Use the cleaned Excel data supplied with the associated research project and place it next to the script, or set `--input` to its path.

The script creates 12 figures in PNG/PDF/SVG formats (4 LDA, 4 grouped response profiles, 4 confusion matrices), one combined multipage PDF, and an LDA results workbook.

## Methods
Each class contains ten observations, randomly split into six training and four held-out prediction observations using a fixed random seed. Standardization and the LDA classifier are fitted **only** on training data. LDA plots show filled training markers, open test markers, and descriptive 95% chi-square covariance ellipses from training scores. Held-out accuracy is based on a **single** stratified random split, not independently collected validation data. Error bars are taken directly from the supplied summary sheets and are not assumed to be either SD or SEM.

The plotting style mirrors the fatty-acid array plots: fixed six-channel colors, white background, and all legends outside the data axes.
