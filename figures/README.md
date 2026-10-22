# figures/

Every figure has a paper-style legend in [`CAPTIONS.md`](CAPTIONS.md), giving what is plotted, the input, the method and n. Each figure is embedded in [`../FINDINGS.md`](../FINDINGS.md) next to its results, and its main takeaway is in the "Figures at a glance" tables at the top.

| Figures | Made by | What they are |
|---|---|---|
| `fig1_composition`, `fig2_umap` | `notebooks/fig1_fig2` | Main Figs 1–2: cells by donor / condition / cell type; UMAP |
| `fig3_ifn_signature`, `fig4_de_heatmap`, `fig5_programs_heatmap` | `notebooks/w4_figures_table` | Main Figs 3–5: IFN signature scores; DE heatmap of shared and specific genes; GSEA programs |
| `fig4_readme_shared`, `fig4_readme_specific` | crops of `fig4_de_heatmap.png` | Fig 4 split into its two blocks, for the main README |
| `m0_paper_fig3a`, `m0_paper_fig3d` | exported from `notebooks/00_eda` outputs | Recreations of Kang 2018 Fig 3a (authors' t-SNE) and Fig 3d (cell-type proportions) |
| `m1b_concordance` | `notebooks/m1b_de` | Our DE genes compared with Kang 2018 Supp Table 1 |
| `m2b_*` | `notebooks/m2_validation` | Validation diagnostics: canonical ISGs (`p1`), expression-matched null (`p4`), label shuffles (`n1`), robustness runs (`r2_r3`), fold changes vs the paper (`v9b`), the recreated paper Fig 3b (`f3b`), effect size vs sample size (`i6`) |

The main figures are saved as `.png` and `.pdf`, and the diagnostics as `.png` only.
