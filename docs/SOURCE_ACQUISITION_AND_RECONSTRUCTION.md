# Source acquisition and reconstruction boundaries

This document gives official source routes without embedding unverified direct-download URLs. It is a retrieval guide, not evidence that source files are bundled in the public code tree or data deposit. Before retrieving or redistributing an upstream source, verify the exact file record and reuse terms.

## Activity CROP-seq baseline

- Official landing pages: GEO sample GSM8787485 (`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8787485`) and series GSE289235 (`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289235`). Retrieve the sample's listed supplementary file through the official GEO page. Do not substitute another GEO file or follow an inferred direct file URL.
- Expected local file: `GSM8787485_EXP224.h5ad`, 1,007,226,520 bytes, SHA-256 `d25a014c52cc27104fff4ad79b2efdb99130da671eda08c11f604c454cbbf23c` (reused from the frozen input manifest; no H5AD was downloaded or included in the data deposit).
- Frozen ancillary inputs: `ACTIVITY_LOEUF_COMPLETE_CASE.tsv` SHA-256 `3568f1e018c38e6d9400184cbffc337668e364349bddeb07e87f82039b7ebaf7`; `activity_annotations.tsv` SHA-256 `7d0fa5f89cd9b89511c95191ec3cbf0f74d4a434180f8d81b62109aebf021134`. Their source tree and exact paths are listed in the private input provenance and source dependency map.
- The H5AD file's public access does not by itself clear redistribution. It is not packaged. The Activity gene-level expression table combines measured expression with project membership, hit, and LOEUF fields; it is not a source-only expression export.
- Rebuild command after source rights and input provenance are cleared:

```sh
arch -x86_64 python3 scripts/build_activity_ntc_baseline.py \
  --h5ad /path/to/GSM8787485_EXP224.h5ad \
  --membership /path/to/ACTIVITY_LOEUF_COMPLETE_CASE.tsv \
  --annotations /path/to/activity_annotations.tsv \
  --outdir /path/to/new_activity_rebuild
```

The script checks the frozen SHA-256 values, sample shape and locked counts, then writes local outputs only. Do not distribute its combined derivative without separate file-level review. It does not run the Activity logistic model; that requires the withheld 1,267-row paired model input and the R runner.

## DepMap Public 24Q4

The linked Zenodo Dataset v1.0.0 includes the original `CRISPRInferredCommonEssentials.csv` from Figshare+ dataset item 27993248. The file is byte-identical to the preserved local source and its exact Figshare+ item record states CC BY 4.0. Retain the citation and original license notice in `THIRD_PARTY_NOTICES.md`.

## Sources that remain external or withheld

| Source | Official route | Dataset v1.0.0 / current scope |
|---|---|---|
| Gordon gene scores | Nature article/source-data page: `https://www.nature.com/articles/s41586-025-10047-5` | Download/review exact file identity and source-data scope; no score table included. |
| ASD-PPI membership/bait/prey | Science article DOI: `https://doi.org/10.1126/science.ady4523` | Exact supplementary item and reuse terms unconfirmed; do not bundle. |
| BioPlex degree/abundance/bait covariates | Official resource: `https://bioplex.hms.harvard.edu/data/` | Official source route only; local file redistribution terms not confirmed. |
| IGVF Perturb-seq | Preprint DOI `10.64898/2026.09.02.748925` and its official supplementary files | CC BY-NC 4.0 source context; exact mixed derived tables not included. |
| Neuron in-vivo scores | Article DOI: `https://doi.org/10.1016/j.neuron.2026.08.006` | Exact table redistribution scope unconfirmed; source scores not included. |
| Wamsley | Dryad record `https://datadryad.org/dataset/doi:10.5061/dryad.4b8gthtkr` | Dryad record terms are CC0; exact local member mapping is pending and no workbook is bundled. |
| Velmeshev | Official article/source repository must be reviewed at exact file level | Derived files and cohort/donor lineage are withheld; never publish donor crosswalks. |

## Figure source tree

The frozen renderer reads hash-bound files relative to the extracted `ASD_PPI_CRISPR_MP_FINAL_FREEZE_AND_ASTRA_HANDOFF_R1_1/` source root; the R1.2 `08_scientific_sources/` subset by itself is not a sufficient generic input tree. The exact script bindings plus the R6 main and R8.2 supplementary source-index paths, frozen hashes, row counts, schemas, and permission states are listed in `metadata/FIGURE_SOURCE_TREE_REQUIREMENTS.tsv` and the review-only module dependency manifest. Local validation found the R8.2 index entry `S12B_RETAINED_DONOR_SAMPLE_CELLTYPE_COUNTS.tsv` missing at its indexed source path. Figure redraw therefore remains blocked. The source files are not bundled because rights remain unresolved. Do not replace expected SHA values with hashes from alternate files.

## S11/S12

`scripts/reproduce_tables.py` is **COPY_ONLY**: it copies supplied clean S11/S12 TSVs byte-for-byte after explicit file-level rights and SHA checks; it does not regenerate tables or recalculate analyses. Zenodo Dataset v1.0.0 contains 12 selected aggregate S11 outputs, but not the complete 30-file S11/S12 source set, gene-level mixed-source tables, or submitted workbooks. The data DOI is https://doi.org/10.5281/zenodo.23230371. Obtain and rights-clear any additional source tables independently before using the copy helper.
