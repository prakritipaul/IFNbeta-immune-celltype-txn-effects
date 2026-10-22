import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from src.m2c_programs import jaccard_matrix, label_programs, run_m2c, top_genes


REPO_ROOT = Path(__file__).resolve().parents[1]
RSCRIPT = "/usr/local/bin/Rscript"
FGSEA_SCRIPT = REPO_ROOT / "R/m2c_fgsea.R"
CELL_TYPES = [f"cell_{i}" for i in range(1, 9)]


def _call(gene, cell_type, log2fc, qvalue, de=None, **extra):
    if de is None:
        de = "up" if qvalue < 0.05 and log2fc >= 1 else (
            "down" if qvalue < 0.05 and log2fc <= -1 else "none"
        )
    return {
        "run": "main",
        "cell_type": cell_type,
        "gene": gene,
        "baseMean": extra.pop("baseMean", 10.0),
        "log2FC": log2fc,
        "lfcSE": extra.pop("lfcSE", 0.5),
        "qvalue": qvalue,
        "tested": extra.pop("tested", True),
        "de": de,
        **extra,
    }


def _run_fgsea(tmp_path, ranks, pathways, seed=0):
    tmp_path.mkdir(parents=True)
    ranks_path = tmp_path / "ranks.csv"
    gmt_path = tmp_path / "sets.gmt"
    out_path = tmp_path / "gsea.csv"
    ranks.to_csv(ranks_path, index=False)
    gmt_path.write_text("".join(
        f"{name}\tna\t" + "\t".join(genes) + "\n"
        for name, genes in pathways.items()
    ))
    subprocess.run(
        [RSCRIPT, str(FGSEA_SCRIPT), str(ranks_path), str(gmt_path),
         str(out_path), str(seed)],
        check=True,
        capture_output=True,
        text=True,
    )
    return pd.read_csv(out_path).sort_values(["cell_type", "pathway"]).reset_index(drop=True)


# Catches wrong sort keys, reversed directions, top-N off-by-one, and padded short lists.
def test_T1_top20_ordering():
    rows = []
    for i in range(25):
        q = 0.001 + (i // 2) * 0.001
        rows.append(_call(f"UP{i:02d}", "A", 1.0 + (i % 2), q))
    rows.extend(_call(f"DOWN{i}", "A", -2.0 - i / 10, 0.02 + i / 1000) for i in range(4))

    got = top_genes(pd.DataFrame(rows), n=20)
    up = got.loc[got["direction"] == "up"].sort_values("rank")
    down = got.loc[got["direction"] == "down"].sort_values("rank")
    expected_up = [f"UP{i:02d}" for i in range(20)]
    for i in range(0, 20, 2):
        expected_up[i], expected_up[i + 1] = expected_up[i + 1], expected_up[i]
    assert up["gene"].tolist() == expected_up
    assert up["rank"].tolist() == list(range(1, 21))
    assert down["gene"].tolist() == [f"DOWN{i}" for i in range(4)]
    assert down["rank"].tolist() == [1, 2, 3, 4]


# Catches ranking non-DE genes, inflated sharing joins, and loss of direction.
def test_T2_only_de_genes_and_sharing_joins():
    rows = [
        _call("A", "cell_1", 2, 0.01), _call("B", "cell_1", 2, 0.02),
        _call("A", "cell_2", 2, 0.01), _call("A", "cell_3", -2, 0.01),
        _call("NOT_Q", "cell_1", 100, 0.05, de="none"),
        _call("NOT_LFC", "cell_2", 0.999, 1e-12, de="none"),
    ]
    got = top_genes(pd.DataFrame(rows), n=20)
    memberships = set(map(tuple, got[["gene", "cell_type", "direction"]].to_numpy()))
    assert memberships == {
        ("A", "cell_1", "up"), ("A", "cell_2", "up"),
        ("A", "cell_3", "down"), ("B", "cell_1", "up"),
    }


# Catches an incorrect Jaccard formula, mixed directions, asymmetry, or omitted cell types.
def test_T3_jaccard():
    up = {
        "cell_1": {"a", "b", "c"}, "cell_2": {"b", "c", "d"},
        "cell_3": {"a", "b", "c"}, "cell_4": {"x"},
        "cell_5": {"e"}, "cell_6": {"f"}, "cell_7": {"g"}, "cell_8": {"h"},
    }
    rows = [_call(g, ct, 2, 0.01) for ct, genes in up.items() for g in genes]
    rows += [_call("down_only", ct, -2, 0.01) for ct in CELL_TYPES]
    calls = pd.DataFrame(rows)
    got = jaccard_matrix(calls, "up")
    assert got.shape == (8, 8)
    assert set(got.index) == set(CELL_TYPES) == set(got.columns)
    assert np.allclose(got.to_numpy(), got.to_numpy().T)
    assert np.allclose(np.diag(got), 1)
    assert got.loc["cell_1", "cell_2"] == 0.5
    assert got.loc["cell_1", "cell_3"] == 1
    assert got.loc["cell_1", "cell_4"] == 0
    assert jaccard_matrix(calls, "down").loc["cell_1", "cell_4"] == 1


# Catches baseMean boundary errors, ranking by LFC/q, and silent invalid-stat or duplicate removal.
def test_T4_ranking_input(tmp_path):
    base = pd.DataFrame([
        _call("LOW", "A", 9, 0.01, baseMean=0.99),
        _call("EDGE", "A", 8, 0.01, baseMean=1.00),
        _call("KEEP", "A", 2, 0.01, baseMean=1.01, lfcSE=0.5),
    ])
    in_dir = tmp_path / "data"
    (in_dir / "processed").mkdir(parents=True)
    (in_dir / "msigdb").mkdir()
    base.to_csv(in_dir / "processed/m1b_de.csv", index=False)
    for source in (REPO_ROOT / "data/msigdb").glob("*.gmt"):
        (in_dir / "msigdb" / source.name).write_bytes(source.read_bytes())
    run_m2c(in_dir, tmp_path / "out", seed=0)
    ranks = pd.read_csv(tmp_path / "out/m2c_ranks.csv")
    assert ranks[["gene", "stat"]].to_dict("records") == [{"gene": "KEEP", "stat": 4.0}]

    invalids = [
        (pd.concat([base, pd.DataFrame([_call("BAD", "A", np.nan, 0.01)])]), "A", "BAD"),
        (pd.concat([base, pd.DataFrame([_call("ZERO", "A", 2, 0.01, lfcSE=0)])]), "A", "ZERO"),
        (pd.concat([base, base.loc[base.gene == "KEEP"]]), "A", "KEEP"),
    ]
    for frame, cell_type, gene in invalids:
        frame.to_csv(in_dir / "processed/m1b_de.csv", index=False)
        with pytest.raises((ValueError, AssertionError), match=rf"(?i){cell_type}.*{gene}|{gene}.*{cell_type}"):
            run_m2c(in_dir, tmp_path / f"bad_{gene}", seed=0)


# Catches reversed NES signs, use of the wrong ranking column, or an ignored seed.
@pytest.mark.r
def test_T5_gsea_sign_oracle(tmp_path):
    genes = [f"G{i:04d}" for i in range(1000)]
    ranks = pd.DataFrame({"cell_type": "A", "gene": genes, "stat": np.arange(500, -500, -1)})
    sets = {"TOP": genes[:30], "BOTTOM": genes[-30:]}
    first = _run_fgsea(tmp_path / "first", ranks, sets, seed=17)
    second = _run_fgsea(tmp_path / "second", ranks, sets, seed=17)
    pdt.assert_frame_equal(first, second)
    by_set = first.set_index("pathway")
    assert by_set.loc["TOP", "NES"] > 0
    assert by_set.loc["BOTTOM", "NES"] < 0
    assert (by_set.loc[["TOP", "BOTTOM"], "padj"] < 0.05).all()


# Catches inclusive-size boundary errors, fuzzy symbol matching, or skipped GMT checksum validation.
@pytest.mark.r
def test_T6_gene_sets(tmp_path):
    genes = [f"G{i:04d}" for i in range(3000)]
    set_genes = np.random.default_rng(17).choice(genes, size=501, replace=False).tolist()
    ranks = pd.DataFrame({"cell_type": "A", "gene": genes, "stat": np.arange(3000, 0, -1)})
    sets = {
        "N14": set_genes[:14], "N15": set_genes[:15], "N500": set_genes[:500],
        "N501": set_genes[:501], "NEAR_CASE": set_genes[:14] + [set_genes[14].lower()],
        "NEAR_PUNCT": set_genes[:14] + [set_genes[14] + "."],
    }
    got = _run_fgsea(tmp_path / "r", ranks, sets)
    assert set(got["pathway"]) == {"N15", "N500"}

    in_dir = tmp_path / "data"
    (in_dir / "processed").mkdir(parents=True)
    (in_dir / "msigdb").mkdir()
    pd.DataFrame([_call("G0001", "A", 2, 0.01)]).to_csv(
        in_dir / "processed/m1b_de.csv", index=False
    )
    for source in (REPO_ROOT / "data/msigdb").glob("*.gmt"):
        (in_dir / "msigdb" / source.name).write_bytes(source.read_bytes())
    hallmark = in_dir / "msigdb/h.all.v2026.1.Hs.symbols.gmt"
    hallmark.write_bytes(hallmark.read_bytes() + b"X")
    with pytest.raises((ValueError, AssertionError), match=r"(?i)checksum|sha256"):
        run_m2c(in_dir, tmp_path / "out", seed=0)


# Catches global/per-collection BH scope and treating padj == 0.05 as affected.
def test_T7_bh_scope():
    gsea = pd.DataFrame([
        {"cell_type": "cell_1", "collection": "H", "pathway": "P1", "NES": 2, "pval": 0.01, "padj": 0.03},
        {"cell_type": "cell_1", "collection": "K", "pathway": "P2", "NES": 2, "pval": 0.02, "padj": 0.03},
        {"cell_type": "cell_1", "collection": "K", "pathway": "P3", "NES": 2, "pval": 0.0499, "padj": 0.0499},
        {"cell_type": "cell_2", "collection": "H", "pathway": "P1", "NES": -2, "pval": 0.05, "padj": 0.05},
    ])
    got = label_programs(
        gsea, cell_types=["cell_1", "cell_2"], padj_threshold=0.05
    ).set_index("pathway")
    assert set(got.index) == {"P1", "P2", "P3"}
    assert got.loc["P1", "cell_1"] == "↑"
    assert got.loc["P1", "cell_2"] == "–"
    assert got.loc["P3", "cell_1"] == "↑"


# Catches incorrect program labels or counting signs from non-significant results.
def test_T8_program_labels():
    rows = []
    patterns = {
        "SHARED": [(ct, 2, 0.01) for ct in CELL_TYPES],
        "PARTIAL": [(ct, 2, 0.01) for ct in CELL_TYPES[:2]],
        "SPECIFIC": [(CELL_TYPES[0], -2, 0.01)] +
                    [(ct, -2, 0.5) for ct in CELL_TYPES[1:]],
        "MIXED": [(CELL_TYPES[0], 2, 0.01), (CELL_TYPES[1], -2, 0.01)],
        "NONE": [(CELL_TYPES[0], 2, 0.05)],
        "IGNORE_SIGN": [(CELL_TYPES[0], 2, 0.01), (CELL_TYPES[1], -2, 0.5)],
        "UNTESTED": [(CELL_TYPES[0], 2, 0.01)] +
                    [(ct, 2, 0.5) for ct in CELL_TYPES[1:-1]],
    }
    for pathway, values in patterns.items():
        rows.extend({"cell_type": ct, "collection": "H", "pathway": pathway,
                     "NES": nes, "pval": padj, "padj": padj}
                    for ct, nes, padj in values)
    got = label_programs(
        pd.DataFrame(rows), cell_types=CELL_TYPES, padj_threshold=0.05
    ).set_index("pathway")
    assert got["label"].to_dict() == {
        "SHARED": "shared", "PARTIAL": "partial", "SPECIFIC": "specific",
        "MIXED": "mixed", "IGNORE_SIGN": "specific", "UNTESTED": "specific",
    }
    assert "NONE" not in got.index
    assert all(got.loc["SHARED", ct] == "↑" for ct in CELL_TYPES)
    assert got.loc["SPECIFIC", CELL_TYPES[0]] == "↓"
    assert all(got.loc["SPECIFIC", ct] == "–" for ct in CELL_TYPES[1:])
    assert got.loc["UNTESTED", CELL_TYPES[0]] == "↑"
    assert all(got.loc["UNTESTED", ct] == "–" for ct in CELL_TYPES[1:-1])
    assert got.loc["UNTESTED", CELL_TYPES[-1]] == "n/t"


# Catches an end-to-end pipeline that loses the established IFN response on repository data.
@pytest.mark.realdata
def test_T9_real_data_tripwire(tmp_path):
    run_m2c(REPO_ROOT / "data", tmp_path, seed=0)
    gsea = pd.read_csv(tmp_path / "m2c_gsea.csv")
    ifna = gsea.loc[gsea["pathway"] == "HALLMARK_INTERFERON_ALPHA_RESPONSE"]
    assert set(ifna["cell_type"]) == {
        "B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells",
        "Dendritic cells", "FCGR3A+ Monocytes", "Megakaryocytes", "NK cells",
    }
    assert (ifna["NES"] > 0).all()
    well_powered = {
        "B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells",
        "FCGR3A+ Monocytes", "NK cells",
    }
    assert (ifna.set_index("cell_type").loc[list(well_powered), "padj"] < 0.05).all()
