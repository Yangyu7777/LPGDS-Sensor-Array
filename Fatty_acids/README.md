# Fatty acid analysis — L-PGDS–DPH sensor array

This directory contains the **original fatty-acid pipeline**, preserved when the repository expanded to include PFAS.

## Studies
- Unsaturation discrimination (C18:0, C18:1, C18:2)
- Chain-length discrimination (C10, C12, C16, C18)
- Fatty-acid mixture-ratio discrimination
- C18:0 concentration series
- C18:2 concentration series

## Files
- `analysis_pipeline.py` — data loading, train-only scaling, LDA, 95% training-score ellipses, grouped response plots, statistics consistency checks, and exports
- `config.example.json` — worksheet-name and figure configuration
- `VALIDATION.md` — historical results and noted differences between raw and summary statistics
- `requirements.txt` — Python dependencies

## Run

From the repository root:

```bash
pip install -r Fatty_acids/requirements.txt
python Fatty_acids/analysis_pipeline.py --input "20260921-ALL区分汇总LDA和柱状图.xlsx" --output-dir Fatty_acids/analysis_output --seed 20260921
```

The input Excel workbook is not committed. Default sheet names and the expected columns are specified by `config.example.json` and `analysis_pipeline.py`.

Each analyte class is split into training (n=6) and held-out prediction (n=4) replicates. All scaling and model fitting occur on training data only. The historical results in `VALIDATION.md` reflect one split and are not repeated independent validation.
