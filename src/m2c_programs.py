"""M2c: shared programs across cell types, at the level of top genes and of gene sets (GSEA).

Contract: audit/contracts/M2c_programs.md (decisions G1-G6). Tests: tests/test_m2c.py (written blind by Codex).
GSEA itself is R/m2c_fgsea.R; this module validates inputs, ranks genes, maps gene-set symbols (B4),
calls R, and labels programs. Invalid input raises ValueError before R runs or anything is written.
"""
import hashlib
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from src.m3_signature import GENE_MAPPING

_ROOT = Path(__file__).resolve().parents[1]
R_SCRIPT = _ROOT / "R/m2c_fgsea.R"
RSCRIPT = "/usr/local/bin/Rscript"
# collection name -> (GMT file under <in_dir>/msigdb, sha256)
GMTS = {
    "hallmark": ("h.all.v2026.1.Hs.symbols.gmt", "eecaf6dad908334ae885406ec72bdc0646d8917588ed7c219fac92fc5363f596"),
    "kegg_legacy": ("c2.cp.kegg_legacy.v2026.1.Hs.symbols.gmt",
                    "4a87c8d1260e637ec174f5960885519d46fa29da650e70d3a3641a4d41c4fa5e"),
}
MIN_BASEMEAN = 1.0  # ranked = baseMean > 1 (M1b "tested")
SYMBOL = {1: "↑", -1: "↓", 0: "–"}


def top_genes(calls, n=20):
    """Top n up and top n down DE genes per cell type (G5).

    In:  calls with cell_type, gene, log2FC, qvalue, de ("up"/"down"/"none"); one row per cell type x gene.
    Out: cell_type, direction, rank, gene, log2FC, qvalue. Only DE genes; ranked by q ascending, then
         |log2FC| descending, then gene name. Fewer than n rows when fewer genes are DE (no padding),
         e.g. 25 up + 4 down -> ranks 1-20 up, 1-4 down.
    """
    _check_unique(calls)
    de = calls.loc[calls["de"].isin(["up", "down"]),
                   ["cell_type", "de", "gene", "log2FC", "qvalue"]].rename(columns={"de": "direction"})
    de = de.assign(abs_lfc=de["log2FC"].abs()).sort_values(
        ["cell_type", "direction", "qvalue", "abs_lfc", "gene"], ascending=[True, True, True, False, True])
    top = de.groupby(["cell_type", "direction"], sort=False).head(n).copy()
    top.insert(2, "rank", top.groupby(["cell_type", "direction"]).cumcount() + 1)
    return top[["cell_type", "direction", "rank", "gene", "log2FC", "qvalue"]].reset_index(drop=True)


def top_sharing(top):
    """Gene x cell type table of top-list membership: ↑ / ↓ / –, plus n_up, n_down (cell types).

    E.g. A in top up of cell_1 and cell_2, top down of cell_3 -> A: cell_1 ↑, cell_2 ↑, cell_3 ↓, n_up 2, n_down 1.
    """
    sign = top.assign(s=top["direction"].map({"up": 1, "down": -1}))
    wide = sign.pivot(index="gene", columns="cell_type", values="s")
    n_up, n_down = (wide == 1).sum(axis=1), (wide == -1).sum(axis=1)
    table = wide.apply(lambda col: col.map(lambda v: SYMBOL[0] if pd.isna(v) else SYMBOL[int(v)]))
    table["n_up"], table["n_down"] = n_up, n_down
    return table.sort_values(["n_up", "n_down"], ascending=False).reset_index()


def jaccard_matrix(calls, direction):
    """Jaccard index |A ∩ B| / |A ∪ B| between the full DE sets (one direction) of every pair of cell types.

    In:  calls (cell_type, gene, de); direction "up" or "down".
    Out: square DataFrame, cell types (all in calls, sorted) on both axes; diagonal 1, or NaN for a cell type with no DE genes in that direction (0/0).
         E.g. {a,b,c} vs {b,c,d} -> 0.5. Two empty sets -> NaN (0/0 is undefined; none in the real data).
    """
    if direction not in ("up", "down"):
        raise ValueError(f"direction must be 'up' or 'down', got {direction!r}")
    cell_types = sorted(calls["cell_type"].unique())
    sets = {ct: set(calls.loc[(calls["cell_type"] == ct) & (calls["de"] == direction), "gene"]) for ct in cell_types}
    return pd.DataFrame([[len(sets[a] & sets[b]) / len(sets[a] | sets[b]) if sets[a] | sets[b] else np.nan
                          for b in cell_types] for a in cell_types],
                        index=pd.Index(cell_types, name="cell_type"), columns=cell_types)


def label_programs(gsea, cell_types, padj_threshold=0.05):
    """One row per program affected (padj < threshold) in at least one cell type (G4).

    In:  gsea with cell_type, collection, pathway, NES, padj (one row per tested cell type x set);
         cell_types = the full list (all 8); a cell type with no gsea rows is n/t everywhere, never dropped.
    Out: collection, pathway, one column per cell type in cell_types (↑ / ↓ affected; – tested, not affected;
         n/t not tested there, i.e. set size outside 15-500 (D2); n/t counts as not affected),
         n_up, n_down, label: shared (affected in every cell type, same sign) / partial (>= 2, same sign) /
         specific (1) / mixed (affected up in some and down in others). Signs of unaffected cells are ignored.
    """
    unknown = sorted(set(gsea["cell_type"]) - set(cell_types))
    if unknown:
        raise ValueError(f"gsea cell types not in cell_types: {unknown}")
    cell_types = sorted(cell_types)
    hit = gsea.loc[gsea["padj"] < padj_threshold]
    sign = hit.assign(s=(hit["NES"] > 0).astype(int) * 2 - 1)
    wide = sign.pivot(index=["collection", "pathway"], columns="cell_type", values="s").reindex(columns=cell_types)
    n_up, n_down = (wide == 1).sum(axis=1), (wide == -1).sum(axis=1)
    n_aff = n_up + n_down
    label = pd.Series("partial", index=wide.index)
    label[n_aff == 1] = "specific"
    label[(n_aff == len(cell_types)) & ((n_up == 0) | (n_down == 0))] = "shared"
    label[(n_up > 0) & (n_down > 0)] = "mixed"
    out = wide.apply(lambda col: col.map(lambda v: SYMBOL[0] if pd.isna(v) else SYMBOL[int(v)]))
    tested = gsea.assign(t=True).pivot(index=["collection", "pathway"], columns="cell_type", values="t").reindex(
        index=wide.index, columns=cell_types)
    out = out.mask(tested.isna(), "n/t")  # D2: set size outside 15-500 in that cell type
    out["n_up"], out["n_down"], out["label"] = n_up, n_down, label
    out.columns.name = None
    return out.reset_index()


def rank_genes(calls):
    """Ranking input for GSEA (G2): tested genes (baseMean > 1), stat = log2FC / lfcSE (DESeq2 Wald statistic).

    Out: cell_type, gene, stat; sorted by cell type, stat descending, then gene (fixes the order of ties).
    Raises ValueError naming cell type and gene for a NaN / infinite stat or zero lfcSE.
    """
    _check_unique(calls)
    ranked = calls.loc[calls["baseMean"] > MIN_BASEMEAN, ["cell_type", "gene", "log2FC", "lfcSE"]]
    bad = ranked.loc[~np.isfinite(ranked["log2FC"]) | ~np.isfinite(ranked["lfcSE"]) | (ranked["lfcSE"] <= 0)]
    if len(bad):
        raise ValueError(f"invalid Wald stat (NaN log2FC/lfcSE or lfcSE <= 0) for cell_type x gene: "
                         f"{list(zip(bad['cell_type'], bad['gene']))}")
    ranked = ranked.assign(stat=ranked["log2FC"] / ranked["lfcSE"])
    return ranked.sort_values(["cell_type", "stat", "gene"], ascending=[True, False, True])[
        ["cell_type", "gene", "stat"]].reset_index(drop=True)


def load_gmt(path, expected_sha256, mapping=GENE_MAPPING):
    """{set: [genes]} from a GMT whose sha256 matches, with B4's symbol mapping applied (exact symbols only)."""
    digest = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    if digest != expected_sha256:
        raise ValueError(f"{Path(path).name}: sha256 checksum {digest} != expected {expected_sha256}")
    sets = {}
    for line in Path(path).read_text().splitlines():
        fields = line.split("\t")
        sets[fields[0]] = [mapping.get(g, g) for g in fields[2:]]
    return sets


def coverage(collections, genes):
    """Per set: n_genes, n_found among `genes` (exact symbol), absent genes (';'-joined)."""
    present = set(genes)
    return pd.DataFrame([
        {"collection": c, "pathway": name, "n_genes": len(g), "n_found": sum(x in present for x in g),
         "absent": ";".join(x for x in g if x not in present)}
        for c, sets in collections.items() for name, g in sets.items()])


def run_gsea(ranks, collections, seed):
    """Write ranks + mapped GMTs to a temp dir, run R/m2c_fgsea.R, return its table (raises with R's output)."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        ranks.to_csv(tmp / "ranks.csv", index=False)
        gmts = []
        for name, sets in collections.items():
            gmts.append(tmp / f"{name}.gmt")
            gmts[-1].write_text("".join(f"{s}\tna\t" + "\t".join(g) + "\n" for s, g in sets.items()))
        proc = subprocess.run([RSCRIPT, str(R_SCRIPT), str(tmp / "ranks.csv"), *map(str, gmts),
                               str(tmp / "gsea.csv"), str(seed)], capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"m2c_fgsea.R failed:\n{proc.stdout}\n{proc.stderr}")
        return pd.read_csv(tmp / "gsea.csv")


def run_m2c(in_dir, out_dir, seed=0):
    """Whole module: <in_dir>/processed/m1b_de.csv (main run) + <in_dir>/msigdb GMTs -> m2c_* files in out_dir.

    Writes m2c_top_genes.csv, m2c_top_sharing.csv, m2c_jaccard.csv, m2c_ranks.csv, m2c_gene_set_coverage.csv,
    m2c_gsea.csv, m2c_programs.csv. All checks (checksums, stats) run before R or any write.
    """
    in_dir, out_dir = Path(in_dir), Path(out_dir)
    de = pd.read_csv(in_dir / "processed/m1b_de.csv")
    calls = de.loc[de["run"] == "main"]
    if calls.empty:
        raise ValueError("no 'main' run rows in m1b_de.csv")
    collections = {c: load_gmt(in_dir / "msigdb" / f, sha) for c, (f, sha) in GMTS.items()}
    ranks = rank_genes(calls)
    top = top_genes(calls)
    gsea = run_gsea(ranks, collections, seed)

    out_dir.mkdir(parents=True, exist_ok=True)
    top.to_csv(out_dir / "m2c_top_genes.csv", index=False)
    top_sharing(top).to_csv(out_dir / "m2c_top_sharing.csv", index=False)
    ranks.to_csv(out_dir / "m2c_ranks.csv", index=False)
    coverage(collections, calls["gene"].unique()).to_csv(out_dir / "m2c_gene_set_coverage.csv", index=False)
    gsea.to_csv(out_dir / "m2c_gsea.csv", index=False)
    label_programs(gsea, cell_types=sorted(calls["cell_type"].unique())).to_csv(out_dir / "m2c_programs.csv", index=False)
    pd.concat({d: jaccard_matrix(calls, d) for d in ("up", "down")}, names=["direction"]).to_csv(
        out_dir / "m2c_jaccard.csv")


def _check_unique(calls):
    dup = calls.loc[calls.duplicated(["cell_type", "gene"], keep=False), ["cell_type", "gene"]]
    if len(dup):
        raise ValueError(f"duplicated cell_type x gene: {sorted(set(zip(dup['cell_type'], dup['gene'])))}")
