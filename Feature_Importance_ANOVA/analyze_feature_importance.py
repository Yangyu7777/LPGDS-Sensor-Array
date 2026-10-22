#!/usr/bin/env python3
"""Reproduce six-channel ANOVA feature importance for fatty acids and PFAS."""
import argparse
import csv
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from artifact_tool import Blob, SpreadsheetFile
from sklearn.feature_selection import SelectKBest, f_classif

CHANNELS = ("WT", "LA", "LB", "LC", "LD", "LE")
SHEETS = {
    "FAs": [
        ("不同链长LDA", "Chain length"),
        ("不饱和度LDA", "Unsaturation"),
        ("不同比例LDA", "Mixture ratio"),
        ("C18浓度LDA", "C18:0 concentration"),
        ("C18-2浓度LDA", "C18:2 concentration"),
    ],
    "PFAS": [
        ("5 PFAS LDA", "PFAS identity"),
        ("PFOA PFOS ratio LDA", "PFOA/PFOS mixture"),
        ("PFOA concentration LDA", "PFOA concentration"),
        ("PFOS concentration LDA", "PFOS concentration"),
    ],
}


def read_task(book, sheet_name):
    data = book.worksheets.get_item(sheet_name).get_range("A1:G200").values
    valid = []
    for row in data:
        if len(row) < 7 or not isinstance(row[0], str) or not row[0].strip():
            continue
        values = row[1:7]
        if all(isinstance(v, (float, int)) and not isinstance(v, bool)
               and np.isfinite(v) for v in values):
            valid.append((row[0].strip(), values))
    if not valid:
        raise ValueError(f"No samples in worksheet: {sheet_name}")
    labels = np.asarray([v[0] for v in valid])
    x = np.asarray([v[1] for v in valid], dtype=float)
    counts = Counter(labels)
    if len(counts) < 2 or min(counts.values()) < 2:
        raise ValueError(f"Insufficient class observations: {sheet_name}: {counts}")
    return x, labels, counts


def analyze(path, group):
    book = SpreadsheetFile.import_xlsx(Blob.load(str(path)))
    outputs = []
    for sheet, task in SHEETS[group]:
        x, y, counts = read_task(book, sheet)
        scores = np.asarray(SelectKBest(f_classif, k="all").fit(x, y).scores_, dtype=float)
        if not np.all(np.isfinite(scores)) or np.any(scores < 0) or scores.sum() <= 0:
            raise ValueError(f"Invalid ANOVA scores: {sheet}")
        outputs.append(dict(group=group, task=task, sheet=sheet, n=len(y),
                            counts=counts, scores=scores,
                            normalized=100 * scores / scores.sum()))
    return outputs


def save_tables(records, out):
    with (out / "feature_importance_detailed.csv").open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["Group", "Task", "Worksheet", "Channel", "N samples",
                    "Class counts", "ANOVA F-score", "Normalized F-score (%)"])
        for r in records:
            for channel, score, percent in zip(CHANNELS, r["scores"], r["normalized"]):
                w.writerow([r["group"], r["task"], r["sheet"], channel, r["n"],
                            "; ".join(f"{k}:{v}" for k,v in sorted(r["counts"].items())),
                            f"{score:.9g}", f"{percent:.6f}"])
    with (out / "feature_importance_matrix.csv").open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["Group", "Task", *CHANNELS])
        for r in records:
            w.writerow([r["group"], r["task"], *(f"{v:.6f}" for v in r["normalized"])])


def plot(rows, out, annotate=False, vmax=70):
    fig, ax = plt.subplots(figsize=(9.2, 1.1 + len(rows) * 0.75))
    fig.patch.set_alpha(0)
    ax.patch.set_alpha(0)
    z = np.stack([r["normalized"] for r in rows])
    im = ax.imshow(z, cmap="Blues", vmin=0, vmax=vmax, aspect="auto", interpolation="nearest")
    ax.set_xticks(range(6), CHANNELS)
    ax.set_yticks(range(len(rows)), [r["task"] for r in rows])
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    if annotate:
        for i, row in enumerate(z):
            for j, value in enumerate(row):
                ax.text(j, i, f"{value:.2f}", ha="center", va="center",
                        color="white" if value >= vmax * .55 else "black", fontsize=9)
    cb = fig.colorbar(im, ax=ax, fraction=.04, pad=.03)
    cb.set_label("Normalized ANOVA F-score (%)")
    cb.outline.set_visible(False)
    fig.tight_layout()
    for suffix in ("png", "svg", "pdf"):
        fig.savefig(out.with_suffix("." + suffix), dpi=350, transparent=True)
    plt.close(fig)


def plot_combined(records, out, annotate=False, vmax=70):
    fig, axes = plt.subplots(2, 1, figsize=(10, 9),
                              gridspec_kw={"height_ratios": (5, 4)})
    fig.patch.set_alpha(0)
    im = None
    for ax, group, label in zip(axes, ("FAs", "PFAS"), ("A  Fatty acids", "B  PFAS")):
        rows = [r for r in records if r["group"] == group]
        z = np.stack([r["normalized"] for r in rows])
        im = ax.imshow(z, cmap="Blues", vmin=0, vmax=vmax, aspect="auto",
                       interpolation="nearest")
        ax.set_xticks(range(6), CHANNELS)
        ax.set_yticks(range(len(rows)), [r["task"] for r in rows])
        ax.set_title(label, loc="left", fontweight="bold")
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        if annotate:
            for i, row in enumerate(z):
                for j, value in enumerate(row):
                    ax.text(j, i, f"{value:.2f}", ha="center", va="center",
                            color="white" if value >= vmax * .55 else "black", fontsize=9)
    fig.subplots_adjust(left=.27, right=.82, top=.95, bottom=.08, hspace=.43)
    cb = fig.colorbar(im, ax=axes, fraction=.03, pad=.035)
    cb.set_label("Normalized ANOVA F-score (%)")
    cb.outline.set_visible(False)
    for suffix in ("png", "svg", "pdf"):
        fig.savefig(out.with_suffix("." + suffix), dpi=350, transparent=True)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fas", type=Path, required=True)
    p.add_argument("--pfas", type=Path, required=True)
    p.add_argument("--out", type=Path, default=Path("results"))
    p.add_argument("--annotate", action="store_true")
    p.add_argument("--vmax", type=float, default=70.0)
    a = p.parse_args()
    if a.vmax <= 0:
        p.error("--vmax must be positive")
    a.out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
    records = analyze(a.fas, "FAs") + analyze(a.pfas, "PFAS")
    save_tables(records, a.out)
    for group in ("FAs", "PFAS"):
        plot([r for r in records if r["group"] == group], a.out / (group + "_heatmap"),
             a.annotate, a.vmax)
    plot_combined(records, a.out / "FAs_PFAS_combined_heatmap", a.annotate, a.vmax)
    for r in records:
        print(r["group"], r["task"], dict(zip(CHANNELS, np.round(r["normalized"], 2))))
    print("Output:", a.out.resolve())


if __name__ == "__main__":
    main()
