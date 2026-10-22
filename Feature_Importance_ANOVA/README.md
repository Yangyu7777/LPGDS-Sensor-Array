# L-PGDS feature importance: fatty acids and PFAS

This folder contains a reproducible, six-channel **SelectKBest (one-way ANOVA F-test)** analysis covering five fatty-acid discrimination tasks and four PFAS discrimination tasks.

## Statistical interpretation
For each task, each channel (WT, LA, LB, LC, LD, LE) is scored independently using `sklearn.feature_selection.f_classif`. The normalized score is `100 × F_i / sum(F_WT, ..., F_LE)`. Scores sum to 100% *within each task* and represent relative univariate discriminative ability, **not** the channel's contribution to LDA classification accuracy. This analysis uses all replicates descriptively; it is separate from LDA training/testing.

## Files
- `analyze_feature_importance.py`: source code for scoring and plotting.
- `requirements.txt`: dependencies.
- `results/feature_importance_detailed.csv`: raw ANOVA scores and normalized percentages.
- `results/feature_importance_matrix.csv`: nine tasks by six channel matrix.
- `results/FAs_PFAS_combined_heatmap.svg`: editable vector heatmap.
- `results/FAs_PFAS_combined_heatmap.png`: preview image.

## Inputs (not duplicated in this folder)
- Fatty acids: `20260922-ALL区分汇总LDA和柱状图.xlsx`
- PFAS: `20261017-ALL LDA and Bar chart(1).xlsx`

Input sheet names and columns are configured in the script. The data matrix assumes category labels in column A and channel responses WT, LA, LB, LC, LD, LE in columns B–G.

## Usage
```
python -m pip install -r requirements.txt
python analyze_feature_importance.py --fas '20260922-ALL区分汇总LDA和柱状图.xlsx' --pfas '20261017-ALL LDA and Bar chart(1).xlsx' --out results
```
Use `--annotate` to put scores inside the heatmap cells. All plots use a shared color scale; darker blue means greater within-task normalized F-score.

**Verification:** calculations should be regenerated from the source Excel workbooks before manuscript submission; manually transcribed older heatmap figures may contain mismatched numbers.
