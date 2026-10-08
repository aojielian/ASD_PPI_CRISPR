# ASD_PPI_CRISPR

[![Zenodo dataset DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23230371.svg)](https://doi.org/10.5281/zenodo.23230371)

**Project code version:** 1.0.0 · **Software license:** [MIT](LICENSE) · **Repository:** [aojielian/ASD_PPI_CRISPR](https://github.com/aojielian/ASD_PPI_CRISPR)

This repository provides original analysis and plotting code, workflow entry points, and source/reproducibility documentation for investigating the functional associations of ASD risk-protein interaction partners across neuronal perturbation assays. The accompanying Zenodo deposit is a **separately versioned dataset**, not a software or manuscript DOI.

## Companion data deposit

**Data DOI:** [10.5281/zenodo.23230371](https://doi.org/10.5281/zenodo.23230371) · **Dataset version:** 1.0.0

The deposited archive contains 45 author-derived aggregate TSVs (12 S11 tables, 20 main-figure panel summaries, and 13 supplementary-figure panel summaries) plus one unchanged DepMap Public 24Q4 common-essential CSV. Original aggregate tables are described by the dataset's CC BY 4.0 terms insofar as the depositors hold the rights; the DepMap CSV keeps its original CC BY 4.0 attribution. This repository does **not** relicense source-study data. See [Data availability](DATA_AVAILABILITY.md), [source inventory](metadata/SOURCES.tsv), [external inputs](metadata/EXTERNAL_INPUTS.tsv), and [third-party notices](THIRD_PARTY_NOTICES.md).

**Scope of reproduction:** the Zenodo data deposit is a selected frozen-aggregate release, not a complete gene-level analysis-ready input bundle. It does not contain the raw single-cell datasets or the complete hash-bound input tree. Consequently, cloning this code and downloading that archive **alone cannot rerun every model or redraw every figure**.

## Quick start and limitations

1. Read [RUN_FIRST.md](RUN_FIRST.md), the [dependency record](environments/README.md), and [source acquisition requirements](docs/SOURCE_ACQUISITION_AND_RECONSTRUCTION.md).
2. Download the published Zenodo data archive, extract it, and verify its actual `SHA256SUMS.txt` manifest:

   ```bash
   python3 scripts/verify_inputs.py --data-dir /path/to/extracted/zenodo_archive
   ```

3. Full Activity/IGVF/in-vivo model runs require independent acquisition and validation of **additional** inputs, as specified in `metadata/EXTERNAL_INPUTS.tsv`. `scripts/run_analysis.py` fails explicitly when required source inputs are missing; this is intentional.
4. `scripts/reproduce_tables.py` is **COPY_ONLY**: it checks user-supplied file-level permission and hash lists and copies clean frozen S11/S12 TSVs. It does not rerun analyses or independently regenerate those tables. The public dataset has 12 selected S11 tables, not the full 30-file S11/S12 source set.
5. `scripts/reproduce_figures.py` depends on the complete excluded, hash-bound `ASD_PPI_CRISPR_MP_FINAL_FREEZE_AND_ASTRA_HANDOFF_R1_1/` input tree, upstream study files and plotting environment. The R8.2 S12B expected source path is unresolved; do not use a similarly named replacement. The deposited figure-panel summaries are not all inputs required for rendering every PDF.

The primary IGVF expression-plus-LOEUF strengthening analysis remains `NOT_RUN` (`P=NA`) because the required qualified baseline expression was unavailable. Activity and in-vivo additional sensitivities did not reach significance after their specified three-test Holm correction. Correlated assay readouts should not be described as independent cross-platform validations.

## Licensing and citation

Original project software is released under the **MIT License**; see [LICENSE](LICENSE) and [LICENSES.md](LICENSES.md). Dependencies and original study inputs have their own terms. The original DepMap data in Zenodo must keep its original attribution.

**For the deposited data, cite the version DOI** [10.5281/zenodo.23230371](https://doi.org/10.5281/zenodo.23230371). **For the software, reference this GitHub repository and its corresponding commit or version tag.** No separate software DOI is claimed. A formal `CITATION.cff` is not included because an exact project-specific author roster and citation metadata have not been verified; do not use the data DOI as a software DOI.
