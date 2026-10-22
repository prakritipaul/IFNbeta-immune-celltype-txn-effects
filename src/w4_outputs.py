"""Wave 4: Fig 3 (IFN signature), Fig 4 (DE heatmap), Fig 5 (programs heatmap) and the final table.

Contract: audit/contracts/W4_figures_table.md (decisions F1–F6, F5 in audit/decisions.md).
Nothing here computes a new result: every plotted or tabulated number is read from an approved module's output,
and each builder asserts that what it shows equals its source file.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch, Rectangle
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

from src.figures import COND_COLORS, CT_COLORS, STACK_ORDER

LINKAGE = "average"  # M2c D2 row
SCORES = {"ifna_score": "A. Hallmark IFN-α response score",
          "ifng_score": "B1. Hallmark IFN-γ response score, full set (200 genes; 73 shared with IFN-α)",
          "ifng_noshared_score": "B2. Hallmark IFN-γ response score, genes shared with IFN-α removed (specificity control)"}
LABEL_COLORS = {"shared": "#7b3294", "shared*": "#c2a5cf", "partial": "#a6dba0", "specific": "#008837",
                "specific*": "#b8e186", "mixed": "#f4a582"}
LABEL_ORDER = ["shared", "shared*", "partial", "mixed", "specific", "specific*"]  # F6 sort order
PROGRAM_LABELS = ["shared", "partial", "specific", "mixed"]
NT_GREY = "#bdbdbd"


# ---------- Fig 3 ----------
def kept_cells(cell_scores, sample_scores):
    """Cells of kept samples only (A1: DCs keep 6 donors, megakaryocytes 4).

    In:  m3_cell_scores (barcode, donor, condition, cell_type, scores); m3_sample_scores (cell_type, donor, condition, ...).
    Out: the rows of cell_scores whose (cell_type, donor, condition) is a kept sample; asserts per-sample cell counts
         and score means equal sample_scores (e.g. B cells / patient_101 / ctrl: 107 cells, mean IFN-α 0.160).
    """
    key = ["cell_type", "donor", "condition"]
    k = cell_scores.merge(sample_scores[key], on=key, how="inner", validate="many_to_one")
    means = k.groupby(key)[list(SCORES)].agg(["mean", "size"])
    s = sample_scores.set_index(key).loc[means.index]
    assert len(means) == len(sample_scores), "a kept sample has no cells"
    assert (means[(list(SCORES)[0], "size")] == s["n_cells"]).all(), "cell counts differ from m3_sample_scores"
    for col in SCORES:
        np.testing.assert_allclose(means[(col, "mean")], s[col], rtol=0, atol=1e-9)
    return k


def plot_fig3(cells, sample_scores):
    """Fig 3: per cell type, ctrl and stim violins of per-cell scores; donor means as dots joined ctrl→stim.

    In:  cells from kept_cells; sample_scores (m3_sample_scores.csv, the donor means plotted).
    Out: Figure, one row per SCORES entry × 8 cell types (STACK_ORDER), y shared within a row.
    """
    fig, axes = plt.subplots(len(SCORES), len(STACK_ORDER), figsize=(17, 3.75 * len(SCORES)), sharey="row")
    for r, (col, title) in enumerate(SCORES.items()):
        for c, ct in enumerate(STACK_ORDER):
            ax = axes[r, c]
            for x, cond in enumerate(COND_COLORS):
                v = ax.violinplot(cells.loc[(cells["cell_type"] == ct) & (cells["condition"] == cond), col],
                                  positions=[x], widths=0.8, showextrema=False, bw_method="scott", points=100)
                for b in v["bodies"]:
                    b.set(facecolor=COND_COLORS[cond], edgecolor="none", alpha=0.45)
            wide = sample_scores[sample_scores["cell_type"] == ct].pivot(index="donor", columns="condition", values=col)
            for _, d in wide.iterrows():
                ax.plot([0, 1], [d["ctrl"], d["stim"]], color="#333", lw=0.8, marker="o", ms=3.5, mfc="white", zorder=3)
            ax.set_xticks([0, 1], list(COND_COLORS))
            ax.set_xlim(-0.6, 1.6)
            ax.grid(False)
            if r == 0:
                ax.set_title(f"{ct}\n(n = {len(wide)} donors)", fontsize=9, color="black")
        axes[r, 0].set_ylabel("score")
        axes[r, 0].text(0, 1.22 if r == 0 else 1.08, title, transform=axes[r, 0].transAxes, fontsize=11, fontweight="bold")
    fig.legend(handles=[Patch(color=COND_COLORS[k], alpha=0.45, label=f"{k} cells") for k in COND_COLORS]
               + [plt.Line2D([], [], color="#333", marker="o", mfc="white", label="donor mean (ctrl→stim)")],
               loc="lower center", ncol=3, frameon=False)
    fig.subplots_adjust(left=0.05, right=0.99, top=0.90, bottom=0.07, wspace=0.18, hspace=0.42)
    return fig


# ---------- shared heatmap helpers ----------
def _column_tree(fig, rect, Z):
    """Dendrogram above heatmap columns; leaves at 5, 15, ... so they line up with column centres."""
    ax = fig.add_axes(rect)
    dendrogram(Z, ax=ax, color_threshold=0, above_threshold_color="#555", no_labels=True)
    ax.set_axis_off()
    return ax


def _heatmap(ax, mat, mask, lim, label):
    """Diverging heatmap of mat (NaN = untested, grey hatched); returns the image."""
    ax.set_facecolor(NT_GREY)
    im = ax.imshow(mat.to_numpy(dtype=float), cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto", interpolation="nearest")
    for y, x in zip(*np.nonzero(mat.isna().to_numpy())):
        ax.add_patch(Rectangle((x - 0.5, y - 0.5), 1, 1, facecolor=NT_GREY, hatch="///", edgecolor="white", lw=0))
    ys, xs = np.nonzero(mask.to_numpy())
    ax.scatter(xs, ys, s=7, c="black", lw=0)
    ax.set_xticks(range(mat.shape[1]), mat.columns, rotation=55, ha="right", fontsize=9)
    ax.tick_params(axis="x", length=0)
    ax.grid(False)
    return im


def _label_bar(ax, labels, colors, title):
    """One-column bar coloured by each row's label."""
    ax.imshow(np.array([[list(colors).index(l)] for l in labels]), aspect="auto", interpolation="nearest",
              cmap=plt.matplotlib.colors.ListedColormap(list(colors.values())), vmin=-0.5, vmax=len(colors) - 0.5)
    ax.set_xticks([0], [title], rotation=55, ha="right", fontsize=9)
    ax.set_yticks([])
    ax.tick_params(axis="x", length=0)


def _group_brackets(ax, groups):
    """Row-group names left of the gene labels: groups = list of (name, first_row, last_row)."""
    for name, a, b in groups:
        ax.annotate("", xy=(-0.30, a - 0.4), xytext=(-0.30, b + 0.4), xycoords=("axes fraction", "data"),
                    arrowprops=dict(arrowstyle="-", color="#555", lw=1))
        ax.text(-0.32, (a + b) / 2, name, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=8)


# ---------- Fig 4 ----------
def jaccard_order(jaccard):
    """Column order and linkage from the up-gene Jaccard matrix (F5): average linkage on 1 − J.

    In:  m2c_jaccard.csv read with index_col=[0, 1] (direction, cell_type).
    Out: (cell types in dendrogram leaf order, linkage matrix).
    """
    m = jaccard.loc["up"]
    m = m.loc[:, list(m.index)]
    assert not m.isna().any().any(), "Jaccard undefined"
    Z = linkage(squareform(1 - m.to_numpy(), checks=True), method=LINKAGE)
    return [m.index[i] for i in dendrogram(Z, no_plot=True)["leaves"]], Z


def fig4_rows(table, calls, col_order, n_shared=25, n_specific=5):
    """Rows of Fig 4 (F3, F4).

    In:  m2_table, m2_calls, column order (specific blocks follow it).
    Out: DataFrame gene, label, group, display (gene with * for shared*/specific*), rank_value.
         Shared: top n_shared of shared + shared* by mean |log2FC| over the 8 cell types (all tested for these genes).
         Specific, per cell type (F3 revised): genes DE only there (specific / specific*) and tested in all 8
         cell types, top n_specific by |log2FC| in that cell type (fewer if fewer exist).
    """
    lfc = calls.pivot(index="gene", columns="cell_type", values="log2FC")
    tested = calls.pivot(index="gene", columns="cell_type", values="tested")
    sh = table.loc[table["shared_response"].isin(["shared", "shared*"]), ["gene", "shared_response"]].copy()
    assert tested.loc[sh["gene"]].all().all(), "a shared gene is untested somewhere"
    sh["rank_value"] = lfc.loc[sh["gene"]].abs().mean(axis=1).to_numpy()
    sh = sh.sort_values("rank_value", ascending=False, kind="stable").head(n_shared).assign(group="shared")
    sp = table.loc[table["shared_response"].isin(["specific", "specific*"])].copy()
    sp = sp[tested.loc[sp["gene"]].all(axis=1).to_numpy()]  # F3 revised: tested in all 8
    sp["group"] = sp[STACK_ORDER].apply(lambda r: [c for c in STACK_ORDER if r[c] != "–"], axis=1)
    assert (sp["group"].str.len() == 1).all()
    sp["group"] = sp["group"].str[0]
    sp["rank_value"] = [abs(lfc.at[g, ct]) for g, ct in zip(sp["gene"], sp["group"])]
    blocks = [sp[sp["group"] == ct].sort_values("rank_value", ascending=False, kind="stable").head(n_specific)
              for ct in col_order]
    rows = pd.concat([sh] + blocks)[["gene", "shared_response", "group", "rank_value"]].reset_index(drop=True)
    rows = rows.rename(columns={"shared_response": "label"})
    rows["display"] = rows["gene"] + np.where(rows["label"].str.endswith("*"), "*", "")
    assert not rows["gene"].duplicated().any()
    return rows


def fig4_matrix(rows, calls, col_order):
    """log2FC (NaN where untested) and DE mask, rows × col_order; asserts every tested cell equals m2_calls.

    Out: (lfc DataFrame, de bool DataFrame), indexed by gene.
    """
    c = calls.set_index(["gene", "cell_type"])
    lfc = calls.pivot(index="gene", columns="cell_type", values="log2FC").loc[rows["gene"], col_order]
    tested = calls.pivot(index="gene", columns="cell_type", values="tested").loc[rows["gene"], col_order].astype(bool)
    de = calls.pivot(index="gene", columns="cell_type", values="de").loc[rows["gene"], col_order] != "none"
    assert not (de & ~tested).any().any(), "DE call on an untested gene"
    lfc = lfc.where(tested)
    for g in lfc.index:  # cell-by-cell check against the source file
        for ct in col_order:
            src = c.loc[(g, ct)]
            assert (np.isnan(lfc.at[g, ct]) and not src["tested"]) or lfc.at[g, ct] == src["log2FC"]
            assert de.at[g, ct] == (src["de"] != "none")
    return lfc, de


def plot_fig4(rows, lfc, de, Z, specific_lim=3):
    """Fig 4: two stacked heatmaps sharing columns (Jaccard dendrogram above).

    Shared block: symmetric limits ±ceil(max |log2FC|); specific block: its own ±specific_lim (values beyond saturate).
    Each block has its own colour bar, `Shared response?` side bar and row-group brackets.
    """
    n, gap = len(rows), 0.5
    H = 0.165 * (n + gap) + 2.6
    fig = plt.figure(figsize=(8.5, H))
    left, width, top_h, bottom = 0.40, 0.40, 0.9 / H, 1.55 / H
    map_h = 1 - bottom - top_h - 0.75 / H
    row_h = map_h / (n + gap)
    is_shared = (rows["group"] == "shared").to_numpy()
    blocks = [(is_shared, np.ceil(np.nanmax(np.abs(lfc.to_numpy()[is_shared])))), (~is_shared, specific_lim)]
    y_top = bottom + map_h
    for k, (sel, lim) in enumerate(blocks):
        r, h = rows[sel].reset_index(drop=True), row_h * sel.sum()
        y0 = y_top - h
        ax = fig.add_axes([left, y0, width, h])
        im = _heatmap(ax, lfc[sel], de[sel], lim, "log2FC")
        ax.set_yticks(range(len(r)), r["display"], fontsize=7, fontstyle="italic")
        if k == 0:
            ax.set_xticks([])
        groups, start = [], 0
        for name, m in r.groupby("group", sort=False).size().items():
            groups.append(("shared (top %d)" % m if name == "shared" else f"specific to\n{name}", start, start + m - 1))
            if start:
                ax.axhline(start - 0.5, color="black", lw=0.8)
            start += m
        _group_brackets(ax, groups)
        bar = fig.add_axes([left + width + 0.01, y0, 0.025, h])
        _label_bar(bar, r["label"], {l: LABEL_COLORS[l] for l in ["shared", "shared*", "specific", "specific*"]},
                   "Shared response?")
        if k == 0:
            bar.set_xticks([])
        cax = fig.add_axes([left + width + 0.07, y0 + h * 0.2, 0.025, h * 0.6])
        fig.colorbar(im, cax=cax, label="log2FC (stim vs ctrl)", extend="neither" if k == 0 else "both")
        y_top = y0 - gap * row_h
    _column_tree(fig, [left, bottom + map_h + 0.05 / H, width, top_h], Z)
    fig.legend(handles=[Patch(color=LABEL_COLORS[k], label=k) for k in ["shared", "shared*", "specific", "specific*"]]
               + [plt.Line2D([], [], marker="o", ls="", color="black", ms=3, label="DE (q < 0.05, |log2FC| ≥ 1)")],
               loc="upper left", bbox_to_anchor=(left + width + 0.045, bottom + map_h + 0.6 / H), frameon=False, fontsize=8)
    fig.suptitle("Top shared and cell-type-specific IFN-β-responsive genes", y=1 - 0.2 / H, fontsize=12)
    return fig


# ---------- Fig 5 ----------
def nes_order(gsea, programs):
    """Column order from NES profiles (F5): average linkage, 1 − Pearson r, over affected programs tested in all 8.

    Out: (cell types in leaf order, linkage, number of programs used).
    """
    nes = gsea.pivot(index="pathway", columns="cell_type", values="NES").loc[programs["pathway"], sorted(STACK_ORDER)]
    complete = nes[nes.notna().all(axis=1)]  # display only (M2c D2): correlation needs complete rows
    Z = linkage(complete.T.to_numpy(), method=LINKAGE, metric="correlation")
    return [complete.columns[i] for i in dendrogram(Z, no_plot=True)["leaves"]], Z, len(complete)


def fig5_matrix(gsea, programs, col_order, padj=0.05):
    """NES matrix and significance mask, rows grouped by label (shared, partial, specific, mixed).

    Within a label: by number of affected cell types, then mean NES (descending).
    Asserts every cell equals m2c_gsea NES, and that arrows rebuilt from NES sign + padj equal m2c_programs.
    Out: (rows DataFrame from programs, NES DataFrame (NaN = n/t), significant bool DataFrame), indexed by pathway.
    """
    nes = gsea.pivot(index="pathway", columns="cell_type", values="NES").reindex(columns=col_order)
    p = gsea.pivot(index="pathway", columns="cell_type", values="padj").reindex(columns=col_order)
    rows = programs.assign(n_aff=programs["n_up"] + programs["n_down"],
                           mean_nes=nes.loc[programs["pathway"]].mean(axis=1).to_numpy(),
                           lab=pd.Categorical(programs["label"], PROGRAM_LABELS, ordered=True))
    assert rows["lab"].notna().all()
    rows = rows.sort_values(["lab", "n_aff", "mean_nes"], ascending=[True, False, False], kind="stable").reset_index(drop=True)
    m, padj_m = nes.loc[rows["pathway"]], p.loc[rows["pathway"]]
    sig = padj_m < padj
    arrows = pd.DataFrame(np.where(m.isna(), "n/t", np.where(sig, np.where(m > 0, "↑", "↓"), "–")),
                          index=m.index, columns=col_order)
    assert (arrows.to_numpy() == rows.set_index("pathway").loc[m.index, col_order].to_numpy()).all(), "arrows differ"
    g = gsea.set_index(["pathway", "cell_type"])["NES"]
    for pw in m.index:
        for ct in col_order:
            assert (np.isnan(m.at[pw, ct]) and (pw, ct) not in g.index) or m.at[pw, ct] == g.loc[(pw, ct)]
    return rows, m, sig


def plot_fig5(rows, m, sig, Z):
    """Fig 5: NES heatmap with column dendrogram (NES profiles), label side bar and row-group brackets."""
    n = len(rows)
    H = 0.16 * n + 2.6
    fig = plt.figure(figsize=(10, H))
    left, width, top_h, bottom = 0.46, 0.34, 0.9 / H, 1.55 / H
    map_h = 1 - bottom - top_h - 0.75 / H
    ax = fig.add_axes([left, bottom, width, map_h])
    _column_tree(fig, [left, bottom + map_h + 0.05 / H, width, top_h], Z)
    lim = np.ceil(np.nanmax(np.abs(m.to_numpy())))
    im = _heatmap(ax, m, sig, lim, "NES")
    for y, x in zip(*np.nonzero(m.isna().to_numpy())):
        ax.text(x, y, "n/t", ha="center", va="center", fontsize=5)
    ax.set_yticks(range(n), [pw.replace("HALLMARK_", "H: ").replace("KEGG_", "K: ").replace("_", " ") for pw in rows["pathway"]],
                  fontsize=7)
    groups, start = [], 0
    for name, k in rows.groupby("label", sort=False).size().items():
        groups.append((f"{name} ({k})", start, start + k - 1))
        if start:
            ax.axhline(start - 0.5, color="black", lw=0.8)
        start += k
    bar = fig.add_axes([left + width + 0.01, bottom, 0.025, map_h])
    pcol = {"shared": LABEL_COLORS["shared"], "partial": LABEL_COLORS["partial"], "specific": LABEL_COLORS["specific"],
            "mixed": LABEL_COLORS["mixed"]}
    _label_bar(bar, rows["label"], pcol, "label")
    for name, a, b in groups:
        bar.text(1.0, (a + b) / 2, " " + name, ha="left", va="center", fontsize=8, transform=bar.get_yaxis_transform())
    cax = fig.add_axes([left + width + 0.12, bottom + map_h * 0.62, 0.02, map_h * 0.25])
    fig.colorbar(im, cax=cax, label="NES")
    fig.legend(handles=[plt.Line2D([], [], marker="o", ls="", color="black", ms=3, label="BH-adjusted p < 0.05"),
                        Patch(facecolor=NT_GREY, hatch="///", edgecolor="white", label="n/t: not tested\n(set size outside 15–500)")],
               loc="upper left", bbox_to_anchor=(left + width + 0.09, bottom + map_h * 0.55), frameon=False, fontsize=8)
    fig.suptitle("IFN-β-responsive programs across cell types (GSEA, Hallmark + KEGG)", y=1 - 0.2 / H, fontsize=12)
    return fig


# ---------- final table (F6) ----------
def final_table(table, calls):
    """results/final_table.csv (F6).

    In:  m2_table (arrows, labels, flags), m2_calls (log2FC, qvalue).
    Out: gene, 8 arrow columns (m2_table order), "Shared response?", n_cell_types, expressed_in_n_cell_types
         (cell types where tested, baseMean > 1), megakaryocyte_lower_confidence,
         then "<cell type> log2FC" and "<cell type> q" per cell type (verbatim from m2_calls; q empty = untested).
         Sorted shared → shared* → partial → mixed → specific → specific*, then n_cell_types (desc),
         then mean |log2FC| over tested cell types (desc), then gene.
    """
    cts = [c for c in table.columns if c in STACK_ORDER]
    assert sorted(cts) == sorted(STACK_ORDER)
    lfc = calls.pivot(index="gene", columns="cell_type", values="log2FC").loc[table["gene"], cts]
    q = calls.pivot(index="gene", columns="cell_type", values="qvalue").loc[table["gene"], cts]
    tested = calls.pivot(index="gene", columns="cell_type", values="tested").loc[table["gene"], cts].astype(bool)
    out = table[["gene"] + cts].copy()
    out["Shared response?"] = table["shared_response"].to_numpy()
    out["n_cell_types"] = table["n_cell_types"].to_numpy()
    out["expressed_in_n_cell_types"] = tested.sum(axis=1).to_numpy()  # tested = baseMean > 1
    out["megakaryocyte_lower_confidence"] = table["megakaryocyte_lower_confidence"].to_numpy()
    for ct in cts:
        out[f"{ct} log2FC"] = lfc[ct].to_numpy()
        out[f"{ct} q"] = q[ct].to_numpy()
    sort_key = lfc.where(tested).abs().mean(axis=1).to_numpy()
    out = (out.assign(_lab=pd.Categorical(out["Shared response?"], LABEL_ORDER, ordered=True), _m=sort_key)
           .sort_values(["_lab", "n_cell_types", "_m", "gene"], ascending=[True, False, False, True], kind="stable")
           .drop(columns=["_lab", "_m"]).reset_index(drop=True))
    assert out["Shared response?"].notna().all() and set(out["Shared response?"]) <= set(LABEL_ORDER)
    # arrows identical to m2_table, gene by gene
    ref = table.set_index("gene").loc[out["gene"], cts]
    assert (out[cts].to_numpy() == ref.to_numpy()).all(), "arrows differ from m2_table"
    # every arrow traceable: ↑ / ↓ exactly where m2_calls says de up / down
    de = calls.pivot(index="gene", columns="cell_type", values="de").loc[out["gene"], cts]
    sym = de.replace({"up": "↑", "down": "↓", "none": "–"})
    assert (out[cts].to_numpy() == sym.to_numpy()).all(), "arrows differ from m2_calls de"
    return out
