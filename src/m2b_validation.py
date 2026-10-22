"""M2b: validation checks on the M1b DE calls (report tier: numbers and plots, only V9 has a pass mark).

Contract: audit/contracts/M2_shared_specific.md (part M2b). Definitions: audit/validation_plan.md.
Checks whose definition is still open (P4, I6, F3b/module inputs, I5 ranking, P2 q method) take the open
choice as a required argument: there are no defaults here for anything that changes a number.
"""
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

from src.m0_qc import ribosomal_mask
from src.m1b_de import classify
from src.m3_signature import GENE_MAPPING, IFNA, MSIGDB_GMT, MSIGDB_SHA256, load_gene_sets, map_gene_sets

_ROOT = Path(__file__).resolve().parents[1]
CANONICAL_ISGS = ["ISG15", "IFI6", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "RSAD2", "IRF7", "STAT1"]
# Supp Table 1 sheet -> our cell type; symbol fixes approved for V9 (decisions.md, V9 symbol mapping; same as notebooks/m1b_de.py)
SHEETS = {"cd19": "B cells", "cd14": "CD14+ Monocytes", "cd4": "CD4 T cells", "cd8": "CD8 T cells",
          "dc": "Dendritic cells", "cd16 mono": "FCGR3A+ Monocytes", "megakaryocyte": "Megakaryocytes", "nk": "NK cells"}
DATE_FIX = {"ENSG00000184640": "SEPT9", "ENSG00000099785": "MARCH2", "ENSG00000125354": "SEPT6",
            "ENSG00000108387": "SEPT4", "ENSG00000145416": "MARCH1"}
SUFFIX_FIX = {"CLN3": "CLN3.1", "MRPS17": "MRPS17.1", "NDUFB8": "NDUFB8.1", "PSMA2": "PSMA2.1", "TAP2": "TAP2.1"}
KEGG_GMT = "data/msigdb/c2.cp.kegg_legacy.v2026.1.Hs.symbols.gmt"
KEGG_SHA256 = "4a87c8d1260e637ec174f5960885519d46fa29da650e70d3a3641a4d41c4fa5e"  # O5


def ifna_genes(genes):
    """Hallmark IFN-alpha set in our symbols (B2 file, B4 mapping; 100% coverage or stop). Out: set of 97 genes."""
    sets = load_gene_sets(_ROOT / MSIGDB_GMT, MSIGDB_SHA256, names=(IFNA,))
    mapped, _ = map_gene_sets(sets, genes, GENE_MAPPING)
    return set(mapped[IFNA])


def p1_isgs(de):
    """P1: canonical ISG calls per cell type.

    In:  main-run de (cell_type, gene, de). Out: cell_type x 10 ISGs of "up"/"down"/"none", plus n_de (count != none).
    """
    sub = de[de["gene"].isin(CANONICAL_ISGS)]
    wide = sub.pivot(index="cell_type", columns="gene", values="de")[CANONICAL_ISGS]
    assert wide.notna().all().all(), "every canonical ISG must have a call in every cell type"
    wide["n_de"] = (wide != "none").sum(axis=1)
    return wide


def _fisher(table):
    """Odds ratio and one-sided (greater) and two-sided Fisher p for a 2x2 [[a, b], [c, d]]."""
    res_g = stats.fisher_exact(table, alternative="greater")
    res_2 = stats.fisher_exact(table, alternative="two-sided")
    return res_g.statistic, res_g.pvalue, res_2.pvalue


def bh(p):
    """Benjamini-Hochberg adjusted p-values (scipy)."""
    return stats.false_discovery_control(np.asarray(p, dtype=float), method="bh")


def p2_p3_fisher(de, ifna):
    """P2 (up-DEGs vs all tested genes) and P3 (up vs down) Fisher tests on Hallmark IFN-alpha, per cell type.

    In:  main-run de (cell_type, gene, tested, de); ifna = set of IFN-alpha genes.
    Out: one row per cell type: n_tested, n_ifna_tested, n_up, up_ifna, n_down, down_ifna,
         p2_odds_ratio, p2_p_greater, p2_p_two_sided, p3_odds_ratio, p3_p_greater, p3_p_two_sided.
         P2 table: rows up / not up (tested genes), columns IFN / not IFN. P3 table: rows up / down.
    """
    rows = []
    for ct, d in de[de["tested"]].groupby("cell_type"):
        ifn = d["gene"].isin(ifna)
        up, down = d["de"] == "up", d["de"] == "down"
        a, b = int((up & ifn).sum()), int((up & ~ifn).sum())
        c, dd = int((~up & ifn).sum()), int((~up & ~ifn).sum())
        e, f = int((down & ifn).sum()), int((down & ~ifn).sum())
        assert a + b + c + dd == len(d)
        or2, p2g, p22 = _fisher([[a, b], [c, dd]])
        or3, p3g, p32 = _fisher([[a, b], [e, f]])
        rows.append({"cell_type": ct, "n_tested": len(d), "n_ifna_tested": int(ifn.sum()), "n_up": a + b, "up_ifna": a,
                     "n_down": e + f, "down_ifna": e, "p2_odds_ratio": or2, "p2_p_greater": p2g, "p2_p_two_sided": p22,
                     "p3_odds_ratio": or3, "p3_p_greater": p3g, "p3_p_two_sided": p32})
    out = pd.DataFrame(rows)
    out.insert(out.columns.get_loc("p2_p_greater") + 1, "p2_q_bh", bh(out["p2_p_greater"]))  # P2 decision: BH over 8, one-sided
    return out


def rederive(df):
    """Re-derive `significant`/`de` with the current DE rule (m1b_de.classify: q < 0.05 and |log2FC| >= 1).

    Out: (reclassified copy, number of `de` calls that differ from the file's).
    """
    out = classify(df)
    return out, int((out["de"] != df["de"]).sum())


def read_runs(path, families):
    """Rows of m1b_runs.csv for runs whose name starts with one of `families` (8 columns only; file is 0.84 GB),
    with `de` re-derived by the current rule.

    Out: (DataFrame run, cell_type, gene, baseMean, log2FC, qvalue, tested, cooks_flag, de; n calls changed).
    """
    cols = ["run", "cell_type", "gene", "baseMean", "log2FC", "qvalue", "tested", "cooks_flag", "de"]
    parts = [c[c["run"].str.startswith(tuple(families))] for c in pd.read_csv(path, usecols=cols, chunksize=1_000_000)]
    return rederive(pd.concat(parts, ignore_index=True))


def n1_counts(runs, fits, de):
    """N1: DE genes per shuffle run, next to the real (main) count.

    In:  runs (shuffle rows), fits (m1b_fits.csv), main de. Out: run, cell_type, n_de, kind ("shuffle"/"real").
    Asserts the expected 35/10/3 shuffles per cell type (contract assumption).
    """
    shuffle_fits = fits[fits["run"].str.startswith("shuffle_")][["run", "cell_type"]]
    n = runs[runs["de"] != "none"].groupby(["run", "cell_type"]).size()
    n = n.reindex(pd.MultiIndex.from_frame(shuffle_fits), fill_value=0).rename("n_de").reset_index()
    expected = {"Dendritic cells": 10, "Megakaryocytes": 3}
    per_ct = n.groupby("cell_type").size()
    assert all(per_ct[ct] == expected.get(ct, 35) for ct in per_ct.index) and len(per_ct) == 8, per_ct
    real = (de["de"] != "none").groupby(de["cell_type"]).sum().rename("n_de").reset_index().assign(run="main")
    return pd.concat([n.assign(kind="shuffle"), real.assign(kind="real")], ignore_index=True)


def r1_lodo(de, lodo):
    """R1 (O4): per cell type, fraction of main DE calls whose log2FC keeps its sign in every LODO refit.

    In:  main de; lodo = runs rows of lodo_* runs (run, cell_type, gene, log2FC).
    Out: (summary per cell type: n_lodo_runs, n_de, n_stable, pct_stable; per-call table: cell_type, gene, de, lodo_stable).
    Asserts 8/6/4 LODO runs per cell type and a non-NaN log2FC for every DE gene in every refit.
    """
    calls = de[de["de"] != "none"][["cell_type", "gene", "de", "log2FC"]]
    m = calls.merge(lodo, on=["cell_type", "gene"], suffixes=("", "_lodo"), validate="one_to_many")
    assert m["log2FC_lodo"].notna().all(), "a DE gene has NaN log2FC in a LODO refit"
    n_runs = lodo.groupby("cell_type")["run"].nunique()
    assert n_runs.to_dict() == {ct: {"Dendritic cells": 6, "Megakaryocytes": 4}.get(ct, 8) for ct in n_runs.index}
    per_gene = m.groupby(["cell_type", "gene"]).size()
    assert len(per_gene) == len(calls) and (per_gene.to_numpy() == n_runs[per_gene.index.get_level_values(0)].to_numpy()).all(), \
        "every DE gene must appear in every LODO run of its cell type"
    m["same"] = np.sign(m["log2FC_lodo"]) == np.sign(m["log2FC"])
    per = m.groupby(["cell_type", "gene"])["same"].all().rename("lodo_stable").reset_index()
    per = calls[["cell_type", "gene", "de"]].merge(per, on=["cell_type", "gene"], validate="one_to_one")
    s = per.groupby("cell_type")["lodo_stable"].agg(n_de="size", n_stable="sum")
    s.insert(0, "n_lodo_runs", n_runs)
    s["pct_stable"] = 100 * s["n_stable"] / s["n_de"]
    return s.reset_index(), per


def overlap(de, other, name):
    """R2/R3: overlap of main DE calls with another run's DE calls, per cell type (genes, either direction).

    In:  main de; other = that run's rows (cell_type, gene, de). Out per cell type: n_main, n_<name>, n_both,
         jaccard, pct_main_recovered (= n_both / n_main), sign_flips (DE in both, opposite direction).
    """
    rows = []
    for ct in sorted(de["cell_type"].unique()):
        a = de[(de["cell_type"] == ct) & (de["de"] != "none")].set_index("gene")["de"]
        b = other[(other["cell_type"] == ct) & (other["de"] != "none")].set_index("gene")["de"]
        both = a.index.intersection(b.index)
        rows.append({"cell_type": ct, "n_main": len(a), f"n_{name}": len(b), "n_both": len(both),
                     "jaccard": len(both) / len(a.index.union(b.index)), "pct_main_recovered": 100 * len(both) / len(a),
                     "sign_flips": int((a[both] != b[both]).sum())})
    return pd.DataFrame(rows)


def paper_supp1():
    """Kang 2018 Supp Table 1 in long form with our symbols (approved V9 mapping).

    Out: cell_type, gene, ensembl, paper_log2FC (one row per paper DE gene per cell type).
    """
    out = []
    for s, d in pd.read_excel(_ROOT / "data/kang2018_supp/MOESM4.xlsx", sheet_name=None).items():
        sym = d["Unnamed: 0"].map(DATE_FIX).fillna(d["featureData.symbol"]).replace(SUFFIX_FIX)
        out.append(pd.DataFrame({"cell_type": SHEETS[s], "gene": sym.to_numpy(), "ensembl": d["Unnamed: 0"].to_numpy(),
                                 "paper_log2FC": d["log2FoldChange"].to_numpy()}))
    out = pd.concat(out, ignore_index=True)
    assert not out.duplicated(["cell_type", "gene"]).any()
    return out


def v9(de, paper):
    """V9 overlap (as in notebooks/m1b_de.py) and V9/V9b log2FC correlations, per cell type.

    In:  main de (with log2FC, log2FC_shrunk, tested); paper from paper_supp1().
    Out: (overlap table: ours_only, in_both, paper_only, opposite_sign, pct_paper_recovered, pct_ours_in_paper, V9;
          correlation table: cell_type, gene_set, n, pearson/spearman/median |paper - ours| for unshrunk and shrunken;
          gene-level table used for the scatter).
         gene_set "in_both" = DE in both lists; "paper_tested" = every paper DEG we tested (baseMean > 1).
    """
    m = paper.merge(de[["cell_type", "gene", "de", "tested", "log2FC", "log2FC_shrunk"]], on=["cell_type", "gene"],
                    how="left", validate="one_to_one", indicator=True)
    assert (m["_merge"] == "both").all(), "every paper gene must be in our results (approved mapping)"
    m["in_both"] = m["de"] != "none"
    rows, cors = [], []
    for ct, g in m.groupby("cell_type"):
        ours = de[(de["cell_type"] == ct) & (de["de"] != "none")]
        nb = int(g["in_both"].sum())
        rows.append({"cell_type": ct, "ours_only": len(ours) - nb, "in_both": nb, "paper_only": len(g) - nb,
                     "opposite_sign": int((np.sign(g.loc[g["in_both"], "log2FC"]) != np.sign(g.loc[g["in_both"], "paper_log2FC"])).sum()),
                     "pct_paper_recovered": round(100 * nb / len(g), 1), "pct_ours_in_paper": round(100 * nb / len(ours), 1)})
        for name, sel in [("in_both", g["in_both"]), ("paper_tested", g["tested"])]:
            s = g[sel]
            assert s[["log2FC", "log2FC_shrunk"]].notna().all().all(), (ct, name)
            r = {"cell_type": ct, "gene_set": name, "n": len(s)}
            for col, lab in [("log2FC", "unshrunk"), ("log2FC_shrunk", "shrunk")]:
                r[f"pearson_{lab}"] = stats.pearsonr(s["paper_log2FC"], s[col]).statistic
                r[f"spearman_{lab}"] = stats.spearmanr(s["paper_log2FC"], s[col]).statistic
                r[f"median_abs_diff_{lab}"] = (s["paper_log2FC"] - s[col]).abs().median()  # scale check: r ignores shrinkage
            cors.append(r)
    ov = pd.DataFrame(rows)
    ov["V9"] = np.where(ov["pct_paper_recovered"] > 70, "pass", "fail")
    return ov, pd.DataFrame(cors), m.drop(columns="_merge")


def authors_modules(paper):
    """Authors' 7-cluster labels from data/kang2018_code/final.gene.gct.txt, mapped to our symbols via Supp Table 1's Ensembl IDs.

    Out: DataFrame ensembl, gene, module (string label as in the file, "0".."6"), row_order (file order).
    """
    g = pd.read_csv(_ROOT / "data/kang2018_code/final.gene.gct.txt", sep="\t", skiprows=2, dtype=str)
    g = g[g["V1"].str.startswith("ENSG")]
    sym = paper.drop_duplicates("ensembl").set_index("ensembl")["gene"]
    assert len(g) == 3055 and g["V1"].is_unique and g["V1"].isin(sym.index).all()
    assert sym.index.is_unique and sym.is_unique
    return pd.DataFrame({"ensembl": g["V1"].to_numpy(), "gene": sym[g["V1"]].to_numpy(),
                         "module": g["7 Clusters"].to_numpy(), "row_order": np.arange(len(g))})


def f3b_matrix(counts, samples, genes, pseudocount, use_dropped_samples):
    """F3b: log2 CPM per sample (+ all-PBMC group per donor x condition), rows z-scored over samples.

    In:  counts (gene x sample raw pseudobulk, m1a_counts.csv indexed by gene); samples (m1a_samples.csv);
         genes to keep; pseudocount (OPEN); use_dropped_samples (OPEN: include the A1-dropped samples, kept=False?).
    Out: genes x columns DataFrame of z-scores (mean 0, population sd 1 per row; ddof=0).
    Routine: 1) all-PBMC = sum of every sample of a donor x condition (all cells, whatever the choice in 2)
             2) cell-type samples, optionally only kept ones  3) CPM on each column's total over all genes
             4) log2(CPM + pseudocount)  5) keep `genes`  6) z-score rows; constant rows raise.
    """
    assert set(samples["sample"]) == set(counts.columns)
    meta = samples.set_index("sample")
    pbmc = counts.T.groupby([meta.loc[counts.columns, "donor"], meta.loc[counts.columns, "condition"]]).sum().T
    pbmc.columns = [f"All PBMCs|{d}|{c}" for d, c in pbmc.columns]
    cols = meta.index[meta["kept"]] if not use_dropped_samples else meta.index
    mat = pd.concat([counts[cols], pbmc], axis=1)
    lcpm = np.log2(mat / mat.sum(axis=0) * 1e6 + pseudocount).loc[list(genes)]
    sd = lcpm.std(axis=1, ddof=0)
    assert (sd > 0).all(), f"constant rows: {list(sd.index[sd == 0])}"
    return lcpm.sub(lcpm.mean(axis=1), axis=0).div(sd, axis=0)


def compare_modules(z, authors, k=7, random_state=0, n_init=25):
    """Modules (O7): k-means on the z-scored rows vs the authors' labels ("method inferred").

    In:  z (genes x samples), authors (gene, module). Out: (labels Series, crosstab ours x authors, ARI on shared genes).
    """
    labels = pd.Series(KMeans(n_clusters=k, random_state=random_state, n_init=n_init).fit_predict(z.to_numpy()),
                       index=z.index, name="our_module")
    a = authors.set_index("gene")["module"]
    common = labels.index.intersection(a.index)
    return labels, pd.crosstab(labels[common], a[common]), adjusted_rand_score(a[common], labels[common])


def kegg_sets():
    """All KEGG legacy sets (O5 file, sha256 checked), with the 6 B4 renames applied (M2b follow-ups).

    Out: {name: set of genes in our symbols}, e.g. a set holding RIGI -> holds DDX58.
    """
    path = _ROOT / KEGG_GMT
    assert hashlib.sha256(path.read_bytes()).hexdigest() == KEGG_SHA256
    return {f[0]: {GENE_MAPPING.get(g, g) for g in f[2:]} for f in (l.split("\t") for l in path.read_text().splitlines())}


def kegg_enrichment(module_genes, universe, sets):
    """Module biology: one-sided Fisher per KEGG set, module genes vs universe (OPEN: which universe), BH over sets.

    Out: set, overlap, set_size_in_universe, odds_ratio, p, q_bh; sorted by p.
    """
    module_genes, universe = set(module_genes), set(universe)
    assert module_genes <= universe
    rows = []
    for name, s in sets.items():
        s = s & universe
        a = len(module_genes & s)
        res = stats.fisher_exact([[a, len(module_genes) - a], [len(s) - a, len(universe) - len(module_genes) - len(s) + a]],
                                 alternative="greater")
        rows.append({"set": name, "overlap": a, "set_size_in_universe": len(s), "odds_ratio": res.statistic, "p": res.pvalue})
    out = pd.DataFrame(rows)
    out["q_bh"] = bh(out["p"])
    return out.sort_values("p").reset_index(drop=True)


def i5_ribo(table, score, n=50):
    """I5: ribosomal (D4 definition) and MT- genes among the top n shared genes of m2_table.csv.

    In:  table (m2_table.csv: gene, shared_response); score = Series gene -> ranking value, higher = more top (OPEN).
    Out: dict n_shared, n_ranked, n_ribo, n_mt, ribo_genes.
    """
    shared = table.loc[table["shared_response"].isin(["shared", "shared*"]), "gene"]
    top = score.reindex(shared).sort_values(ascending=False).head(n).index
    return {"n_shared": len(shared), "n_ranked": len(top), "n_ribo": int(ribosomal_mask(top).sum()),
            "n_mt": int(pd.Index(top).str.startswith("MT-").sum()), "ribo_genes": list(top[ribosomal_mask(top)])}


def size_factors(counts):
    """DESeq2 median-of-ratios size factors (estimateSizeFactors, type="ratio") for one fit's samples.

    In:  counts genes x samples (raw). Out: Series per sample. Genes with a zero in any sample are left out of the
         reference (DESeq2: geometric mean 0). E.g. columns [1,2],[2,4] -> sf 1/sqrt(2), sqrt(2).
    """
    with np.errstate(divide="ignore"):
        logc = np.log(counts.to_numpy(dtype=float))
    ref = logc.mean(axis=1)
    ok = np.isfinite(ref)
    assert ok.sum() > 0
    return pd.Series(np.exp(np.median(logc[ok] - ref[ok, None], axis=0)), index=counts.columns)


def p4_features(counts, samples, tested_genes, detection):
    """P4 matching features for one cell type.

    In:  counts (genes x samples, raw, all samples); samples (m1a_samples rows of THIS cell type);
         tested_genes; detection = Series gene -> fraction of the cell type's control cells (kept donors) with count > 0.
    Out: DataFrame indexed by tested gene: control_mean (mean DESeq2-normalised count over kept control samples;
         size factors from all kept samples of the cell type and all genes, as in the main fit), detection, and their
         rank-percentiles pct_mean, pct_detection (ties averaged; rank of log(mean) = rank of mean).
    """
    kept = samples[samples["kept"]]
    sf = size_factors(counts[kept["sample"]])
    ctrl = kept.loc[kept["condition"] == "ctrl", "sample"]
    norm = counts.loc[list(tested_genes), ctrl] / sf[ctrl]
    f = pd.DataFrame({"control_mean": norm.mean(axis=1), "detection": detection.reindex(norm.index)})
    assert f.notna().all().all()
    f["pct_mean"] = f["control_mean"].rank(method="average", pct=True)
    f["pct_detection"] = f["detection"].rank(method="average", pct=True)
    return f


def p4_null(features, up_genes, de_genes, ifna, k=30, n_draws=1000, seed=0):
    """P4: expression-matched null for the number of IFN-alpha genes among the up-DEGs.

    In:  features from p4_features (index = all tested genes); up_genes; de_genes (up + down, for diagnostics); ifna set.
    Routine: per up-DEG, all tested genes sorted by Euclidean distance on (pct_mean, pct_detection), focal excluded
             (stable sort: ties keep file order); pool = first k. Per draw (one rng = default_rng(seed) per call):
             visit up-DEGs in a random order, pick one unused gene uniformly from its pool; if every pool gene is
             used, take the nearest unused gene outside the pool. Statistic = IFN-alpha genes in the matched set.
    Out: (dict n_up, observed, null_mean, null_median, null_max, p_empirical = (1 + #null >= observed) / (n_draws + 1),
          pct_null_genes_de, pct_null_genes_ifna, pct_degs_ifna, n_pool_exhausted; array of null counts).
    """
    genes = features.index
    pos = genes.get_indexer(list(up_genes))
    assert (pos >= 0).all() and len(set(pos)) == len(pos)
    xy = features[["pct_mean", "pct_detection"]].to_numpy()
    d = np.sqrt(((xy[pos, None, :] - xy[None, :, :]) ** 2).sum(axis=2))
    d[np.arange(len(pos)), pos] = np.inf  # focal excluded: sorts last
    order = np.argsort(d, axis=1, kind="stable")[:, :-1]
    is_ifn, is_de = genes.isin(ifna), genes.isin(list(de_genes))
    rng = np.random.default_rng(seed)
    null, n_de_null, exhausted = np.empty(n_draws, int), 0, 0
    for t in range(n_draws):
        used = np.zeros(len(genes), bool)
        for i in rng.permutation(len(pos)):
            pool = order[i, :k]
            avail = pool[~used[pool]]
            if len(avail):
                g = rng.choice(avail)
            else:
                exhausted += 1
                rest = order[i, k:]
                g = rest[~used[rest]][0]
            used[g] = True
        null[t] = int(is_ifn[used].sum())
        n_de_null += int(is_de[used].sum())
    observed, total = int(is_ifn[pos].sum()), n_draws * len(pos)
    return ({"n_up": len(pos), "observed": observed, "null_mean": null.mean(), "null_median": float(np.median(null)),
             "null_max": int(null.max()), "p_empirical": (1 + int((null >= observed).sum())) / (n_draws + 1),
             "pct_null_genes_de": 100 * n_de_null / total, "pct_null_genes_ifna": 100 * null.sum() / total,
             "pct_degs_ifna": 100 * observed / len(pos), "n_pool_exhausted": exhausted}, null)


def i6_effects(counts, samples, de):
    """I6 (M2b follow-ups): per donor x cell type (kept samples), the effect size of that cell type's DE genes
    vs log of the mean of ctrl and stim n_cells / total_counts.

    In:  counts (genes x samples raw), samples (m1a_samples), main de. CPM = count / column total x 1e6.
    Out: one row per donor x cell type: n_de, median_log2_ratio_aligned (headline: median over DE genes of
         sign(call) x log2((stim CPM + 1) / (ctrl CPM + 1)), so up and down genes don't cancel), median_log2_ratio
         (same without the sign), log_n_cells, log_total_counts.
         E.g. an up gene 3 -> 15 CPM gives +2; a down gene 15 -> 3 CPM also gives +2 aligned.
    """
    tot = counts.sum(axis=0)
    assert (tot[samples["sample"]].to_numpy() == samples["total_counts"].to_numpy()).all()
    lcpm = np.log2(counts / tot * 1e6 + 1)
    rows = []
    for (ct, donor), s in samples[samples["kept"]].groupby(["cell_type", "donor"]):
        s = s.set_index("condition")
        calls = de[(de["cell_type"] == ct) & (de["de"] != "none")].set_index("gene")["de"]
        r = lcpm.loc[calls.index, s.loc["stim", "sample"]] - lcpm.loc[calls.index, s.loc["ctrl", "sample"]]
        assert np.isfinite(r).all()
        rows.append({"cell_type": ct, "donor": donor, "n_de": len(calls),
                     "median_log2_ratio_aligned": (r * np.where(calls == "up", 1, -1)).median(),
                     "median_log2_ratio": r.median(),
                     "log_n_cells": np.log(s["n_cells"].mean()), "log_total_counts": np.log(s["total_counts"].mean())})
    return pd.DataFrame(rows)
