# Contract M3: IFN signature scores

**Status:** APPROVED 2026-10-01 (decisions in `audit/decisions.md`, Wave 1 section). **Wave 1**, parallel with M1a and Figs 1–2. Tier: produces results (all 8 steps + planted-bug check). Fig 3 itself is drawn in Wave 4.

## Goal
Score every cell for Hallmark IFN-α response (the signature) and Hallmark IFN-γ response (specificity control), and summarise per donor × cell type × condition. Checks P5 and P6.

## Inputs
| File | Check on load |
|---|---|
| `data/processed/m0_cells.h5ad` | As in M1a |
| `data/msigdb/h.all.v2026.1.Hs.symbols.gmt` (downloaded from `data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/`) | sha256 `eecaf6da…f596` (B2); both sets present |

## Outputs
| File | Contents |
|---|---|
| `src/m3_signature.py` | Load gene sets, normalise, score, summarise. Every parameter explicit |
| `data/processed/m3_cell_scores.csv` | Per cell: `barcode, donor, condition, cell_type, ifna_score, ifng_score, ifng_noshared_score` |
| `data/processed/m3_sample_scores.csv` | Per donor × cell type × condition: `n_cells`, summary of each score (B3) |
| `audit/m3_gene_sets.md` | Set sizes, genes found / missing in our data (coverage %), genes shared by both sets |
| `notebooks/m3_signature.py` | Executed: coverage, score violins per cell type × condition, paired stim − ctrl per cell type (P5), IFN-α vs IFN-γ increase (P6), score vs library size (I6) |

## Steps
1. Load the two sets; apply the 6 approved mappings (B4); report coverage. Any other missing gene = stop and report (coverage is 100% after mapping).
2. Normalise for scoring (B1); score with the chosen method (B1), fixed seed.
3. IFN-γ scored twice: full set, and with genes shared with IFN-α removed (P6).
4. Summarise per sample (B3); median donor-paired stim − ctrl per cell type. Numbers, no pass marks.

## Assumptions (asserted)
- One score per cell per set; no NaNs (a NaN = stop and report).
- Gene-set symbols match `var_names` exactly after the 6 approved mappings (B4); any other mismatch = stop and ask.

## Validation checks owned
P5, P6, I6 (score vs library size).

## Open decisions
None. B1–B4 decided (`audit/decisions.md`).

## Return package
As M1a.
