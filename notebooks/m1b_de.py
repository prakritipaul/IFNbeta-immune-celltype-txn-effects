# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # M1b: pseudobulk differential expression
#
# Contract: `audit/contracts/M1b_de.md`. Python: `src/m1b_de.py`; DESeq2 + qvalue: `R/m1b_deseq2.R` (tested in `tests/test_m1b.py`).
# This notebook only reads their output files and shows them. **DE** = q < 0.05 and |log2FC| ≥ 1 (unshrunk log2FC, treatment vs control).
# No verdicts here: robustness and negative-control checks (R1–R3, N1) are computed in M2.

# %%
import sys
from pathlib import Path

ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from IPython.display import display

from src.m0_qc import ribosomal_mask

pd.set_option("display.width", 160, "display.max_columns", 20)
PROC = ROOT / "data/processed"
FILES = ["m1b_de.csv", "m1b_runs.csv", "m1b_fits.csv"]
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#e6e6e3", "grid.linewidth": 0.6, "axes.axisbelow": True, "figure.dpi": 110})

de = pd.read_csv(PROC / "m1b_de.csv")
runs = pd.read_csv(PROC / "m1b_runs.csv", usecols=["run", "cell_type", "gene", "de", "pvalue", "cooks_flag"])
fits = pd.read_csv(PROC / "m1b_fits.csv")
CELL_TYPES = sorted(de["cell_type"].unique())

# %% [markdown]
# ## 1. Every output file: shape, head and tail

# %%
for f, df in [("m1b_de.csv", de), ("m1b_fits.csv", fits)]:
    print(f"{f}: {df.shape[0]} rows × {df.shape[1]} columns")
    display(df.head()); display(df.tail())

# m1b_runs.csv is ~0.8 GB: shape from the column subset already loaded, head/tail of all columns read cheaply
n_cols = len(pd.read_csv(PROC / "m1b_runs.csv", nrows=0).columns)
print(f"m1b_runs.csv: {len(runs)} rows × {n_cols} columns")
display(pd.read_csv(PROC / "m1b_runs.csv", nrows=5))
display(pd.read_csv(PROC / "m1b_runs.csv", skiprows=range(1, len(runs) - 4)))

# %% [markdown]
# ## 2. DE genes per cell type (main run), next to Kang 2018 Supp Table 1
# Paper counts = rows per sheet of `data/kang2018_supp/MOESM4.xlsx`, split by the sign of its log2FoldChange. Counts only; gene-level overlap (V9) is in 2b.

# %%
SHEETS = {"cd19": "B cells", "cd14": "CD14+ Monocytes", "cd4": "CD4 T cells", "cd8": "CD8 T cells",
          "dc": "Dendritic cells", "cd16 mono": "FCGR3A+ Monocytes", "megakaryocyte": "Megakaryocytes", "nk": "NK cells"}
paper = pd.read_excel(ROOT / "data/kang2018_supp/MOESM4.xlsx", sheet_name=None)
paper_counts = pd.DataFrame({SHEETS[s]: {"paper_up": int((d["log2FoldChange"] > 0).sum()),
                                          "paper_down": int((d["log2FoldChange"] < 0).sum())} for s, d in paper.items()}).T
ours = de.groupby("cell_type")["de"].value_counts().unstack(fill_value=0).reindex(columns=["up", "down"], fill_value=0)
ours.columns = ["ours_up", "ours_down"]
counts_table = ours.join(paper_counts)
counts_table["ours_total"] = counts_table["ours_up"] + counts_table["ours_down"]
counts_table["paper_total"] = counts_table["paper_up"] + counts_table["paper_down"]
counts_table["genes_tested"] = de.groupby("cell_type")["tested"].sum()
counts_table["significant"] = de.groupby("cell_type")["significant"].sum()
display(counts_table)

# %% [markdown]
# ## 2b. Gene-level overlap with Supp Table 1 (V9)
# `in_both` = DE in both lists (any direction); `opposite_sign` = of those, how many disagree in direction.
# `% paper recovered` is V9 (✓ if > 70%). Paper symbols mapped to our names: 5 Excel date-mangled (by Ensembl ID) and 4 make-unique `.1` names (as PSMA2, B4).

# %%
# Excel turned these symbols into dates; mapped by the sheet's Ensembl ID
DATE_FIX = {"ENSG00000184640": "SEPT9", "ENSG00000099785": "MARCH2", "ENSG00000125354": "SEPT6",
            "ENSG00000108387": "SEPT4", "ENSG00000145416": "MARCH1"}
SUFFIX_FIX = {"CLN3": "CLN3.1", "MRPS17": "MRPS17.1", "NDUFB8": "NDUFB8.1", "PSMA2": "PSMA2.1", "TAP2": "TAP2.1"}
rows = []
for s, d in paper.items():
    ct = SHEETS[s]
    sym = d["Unnamed: 0"].map(DATE_FIX).fillna(d["featureData.symbol"]).replace(SUFFIX_FIX)
    p_sign = pd.Series(np.sign(d["log2FoldChange"].to_numpy()), index=sym)
    assert p_sign.index.is_unique and p_sign.index.isin(de["gene"]).all(), ct
    o = de[(de["cell_type"] == ct) & (de["de"] != "none")].set_index("gene")["log2FC"].pipe(np.sign)
    both = o.index.intersection(p_sign.index)
    rows.append({"cell_type": ct, "ours_only": len(o) - len(both), "in_both": len(both),
                 "paper_only": len(p_sign) - len(both), "opposite_sign": int((o[both] != p_sign[both]).sum()),
                 "% paper recovered": round(100 * len(both) / len(p_sign), 1),
                 "% ours in paper": round(100 * len(both) / len(o), 1)})
overlap = pd.DataFrame(rows).set_index("cell_type").loc[counts_table.index]
overlap["V9 (>70%)"] = np.where(overlap["% paper recovered"] > 70, "✓", "✗")
display(counts_table.join(overlap))

venn = overlap[["ours_only", "in_both", "paper_only"]]
fig, axes = plt.subplots(1, 2, figsize=(14, 4), sharey=True)
for ax, data, xlabel in [(axes[0], venn, "DE genes"),
                         (axes[1], 100 * venn.div(venn.sum(axis=1), axis=0), "% of genes DE in either list")]:
    left = np.zeros(len(data))
    for col, c in zip(data.columns, ["#1b74ac", "#7a5aa6", "#dd1b26"]):
        bars = ax.barh(data.index, data[col], left=left, color=c, label=col.replace("_", " "))
        if ax is axes[1]:  # label segments wide enough to hold the text
            ax.bar_label(bars, labels=[f"{v:.0f}%" if v >= 4 else "" for v in data[col]],
                         label_type="center", color="white", fontsize=8)
        left += data[col].to_numpy()
    ax.set_xlabel(xlabel)
axes[0].invert_yaxis(); axes[1].legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.12))
fig.suptitle("Our DE genes vs Kang 2018 Supp Table 1 (Venn counts as bars)")
plt.tight_layout()
fig.savefig(ROOT / "figures/m1b_concordance.png", dpi=200, bbox_inches="tight"); plt.show()

# %% [markdown]
# ## 3. Fits: dispersion trend used, samples, genes tested (every run × cell type)

# %%
display(fits[fits["run"] == "main"])
print("fit_type_used over all runs × cell types:")
display(fits.groupby(["fit_type_used"]).size().rename("fits"))
display(fits[fits["fit_type_used"] != "parametric"])

# %% [markdown]
# ## 4. NA p-values and Cook's flags
# `cooks_flag` = DESeq2 replaced an outlier count, or results() set the p-value NA by its Cook's filter. Every other NA p-value must be an all-zero gene (the R script stops otherwise).
# DESeq2 computes Cook's distances only for samples whose design cell (unique row of the model matrix) holds ≥ 3 samples; in `~ donor + condition` every sample is alone in its cell, so paired runs cannot flag anything.

# %%
na = runs.assign(na_p=runs["pvalue"].isna())
na_table = na.groupby(["run", "cell_type"])[["na_p", "cooks_flag"]].sum()
display(na_table.loc[["main", "no_outliers", "unpaired"]].unstack(0))
print("Cook's-flagged genes, all runs:", int(runs["cooks_flag"].sum()))
print("NA p-value genes with any non-zero count, main run:",
      int((de["pvalue"].isna() & (de["baseMean"] > 0)).sum()))

# %% [markdown]
# ## 5. Volcano plots (main run, tested genes)
# Red: up; blue: down; grey: not DE. Dashed: |log2FC| = 1. y = −log10 q.

# %%
fig, axes = plt.subplots(2, 4, figsize=(16, 7.5), sharex=True)
COLORS = {"up": "#dd1b26", "down": "#1b74ac", "none": "#b0b0ad"}
for ax, ct in zip(axes.flat, CELL_TYPES):
    d = de[(de["cell_type"] == ct) & de["tested"] & de["qvalue"].notna()]
    for call in ["none", "down", "up"]:
        s = d[d["de"] == call]
        ax.scatter(s["log2FC"], -np.log10(s["qvalue"].clip(lower=1e-300)), s=4, c=COLORS[call], alpha=0.6, lw=0,
                   label=f"{call} ({len(s)})")
    ax.axvline(-1, ls="--", c="k", lw=0.6); ax.axvline(1, ls="--", c="k", lw=0.6)
    ax.set_title(ct); ax.legend(fontsize=7, loc="upper left", markerscale=3)
for ax in axes[1]:
    ax.set_xlabel("log2FC (stim vs ctrl)")
for ax in axes[:, 0]:
    ax.set_ylabel("−log10 q")
plt.tight_layout(); plt.show()

# %% [markdown]
# ## 6. Top 20 DE genes per cell type (main run, by q, then |log2FC|)
# CLAUDE.md sanity check: classic ISGs (ISG15, IFIT*, MX1, OAS*, STAT1 …) should dominate. Last row counts ribosomal genes among the top 50 (`ribosomal_mask`, decision D4).

# %%
top = {}
ribo50 = {}
for ct in CELL_TYPES:
    d = de[(de["cell_type"] == ct) & (de["de"] != "none")].copy()
    d["abs"] = d["log2FC"].abs()
    d = d.sort_values(["qvalue", "abs"], ascending=[True, False])
    top[ct] = [f"{g} ({'+' if l > 0 else '−'}{abs(l):.1f})" for g, l in zip(d["gene"].head(20), d["log2FC"].head(20))]
    top[ct] += [""] * (20 - len(top[ct]))
    ribo50[ct] = int(ribosomal_mask(d["gene"].head(50)).sum())
top_table = pd.DataFrame(top, index=range(1, 21))
top_table.loc["RPL/RPS in top 50"] = ribo50
display(top_table)

# %% [markdown]
# ## 7. DE counts per run (up + down), per cell type
# LODO and shuffle columns give min–max over their runs. Table only; interpretation (R1–R3, N1) is M2's.

# %%
per_run = runs[runs["de"] != "none"].groupby(["run", "cell_type"]).size()
per_run = per_run.reindex(pd.MultiIndex.from_frame(fits[["run", "cell_type"]]), fill_value=0)
per_run = per_run.rename("n_de").reset_index()
per_run["family"] = per_run["run"].str.replace(r"_.*", "", regex=True).where(
    per_run["run"].str.match(r"^(lodo|shuffle)_"), per_run["run"])


def _range(s):
    return f"{s.min()}–{s.max()} (n={len(s)})"


summary = pd.DataFrame({
    "main": per_run[per_run["family"] == "main"].set_index("cell_type")["n_de"],
    "no_outliers": per_run[per_run["family"] == "no_outliers"].set_index("cell_type")["n_de"],
    "unpaired": per_run[per_run["family"] == "unpaired"].set_index("cell_type")["n_de"],
    "lodo min–max": per_run[per_run["family"] == "lodo"].groupby("cell_type")["n_de"].agg(_range),
    "shuffle min–max": per_run[per_run["family"] == "shuffle"].groupby("cell_type")["n_de"].agg(_range),
})
display(summary)

# %% [markdown]
# ## 8. Do technical covariates threaten the DE calls? Variance fractions of DE vs non-DE genes
#
# M1a's variance partitioning (section 5 of `m1a_pseudobulk`) shows condition explaining a median of ~0 of each gene's
# variance, and the technical covariates (`log n_cells`, `log total_counts`) explaining more. Three things to keep in mind:
#
# 1. **The median gene doesn't respond to IFN-β.** Most genes are not DE, so condition's share is ~0 for them. For the
#    genes we call DE, condition explains the majority of the variance (median 57–83% per cell type), and the
#    technical covariates only ~3–9%.
# 2. **Donor's share is not really 0.** The covariates are per-sample and almost identical to donor (canonical
#    correlation 0.80–0.99, M1a section 6): a big donor contributes many cells and counts to every cell type.
#    With two terms carrying the same information, variance partitioning can't split it, and the covariates took it.
#    Adding them back (donor + covariates, bottom panel) gives donor's real share: median ~15–34% per gene.
# 3. **Why this doesn't threaten the DE calls:** `~ donor + condition` compares each donor's stim with their own ctrl,
#    so donor-level differences (including `n_cells` and `total_counts`) are removed. Within a donor the covariates
#    barely differ between conditions (correlation with condition ≤ 0.19; megakaryocytes 0.31), and size factors
#    already correct for depth. Direct evidence: label shuffles give 0–2 DE genes (section 7).
#
# **The plot also shows:** DE genes have *less* technical-covariate and donor variance than non-DE genes (bottom two
# panels), and a small tail of DE genes where condition explains < 10% (197 of 4,947, 4%; 105 of them in CD14+
# monocytes). Those are either noisy low-count genes (large residual) or genes with big donor differences where the
# paired test still finds a consistent within-donor shift. Leave-one-donor-out (R1, in M2) will show whether they hold.
#
# Remaining caveats: megakaryocytes (covariates track condition a little more, 4 donors), and condition = 10x run,
# which no model can separate.

# %%
import seaborn as sns

vp = pd.read_csv(PROC / "m1a_varpart.csv")
vp = vp.merge(de[["cell_type", "gene", "de"]], on=["cell_type", "gene"], how="left", validate="one_to_one")
assert vp["de"].notna().all(), "every gene in variance partitioning should have a main-run DE call"
vp["DE status"] = np.where(vp["de"] != "none", "DE", "not DE")
vp["technical covariates"] = vp["log_n_cells"] + vp["log_total_counts"]
vp["donor + covariates"] = vp["donor"] + vp["technical covariates"]

panels = [("condition", "Condition"), ("technical covariates", "Technical covariates (log n_cells + log total_counts)"),
          ("donor + covariates", "Donor + technical covariates (donor's real share)")]
fig, axes = plt.subplots(len(panels), 1, figsize=(12, 10), sharex=True)
for ax, (col, title) in zip(axes, panels):
    sns.violinplot(data=vp, x="cell_type", y=col, hue="DE status", hue_order=["DE", "not DE"],
                   palette={"DE": "#dd1b26", "not DE": "#b0b0ac"}, split=True, inner="quart", cut=0,
                   density_norm="width", linewidth=0.8, ax=ax)
    ax.set_title(title, loc="left", fontsize=10)
    ax.set_ylabel("variance fraction"); ax.set_xlabel(""); ax.set_ylim(0, 1)
    ax.legend(loc="upper right", fontsize=8, frameon=False) if ax is axes[0] else ax.get_legend().remove()
axes[-1].tick_params(axis="x", rotation=30)
plt.tight_layout(); plt.show()

table = vp.groupby(["cell_type", "DE status"])[["condition", "technical covariates", "donor + covariates"]].median()
table["genes"] = vp.groupby(["cell_type", "DE status"]).size()
display(table.round(3))
