# Run first

The original project code is released under the [MIT License](LICENSE). This does **not** grant any rights to download or redistribute source-study data beyond their original terms.

1. Read `metadata/REPRODUCIBILITY_LEVELS.tsv`, `metadata/EXTERNAL_INPUTS.tsv`, and `docs/SOURCE_ACQUISITION_AND_RECONSTRUCTION.md` before trying any analysis.
2. The separately released Zenodo Dataset v1.0.0 uses DOI [10.5281/zenodo.23230371](https://doi.org/10.5281/zenodo.23230371) and contains selected aggregate results, not every input needed to rerun all models.
3. To check an extracted copy of that ZIP, run `python3 scripts/verify_inputs.py --data-dir /path/to/zenodo_archive`. This validates its `SHA256SUMS.txt`. An older `metadata/MANIFEST_SHA256.tsv` checksum format is also recognized if provided.
4. `scripts/reproduce_tables.py` is **COPY_ONLY**. It requires supplied S11/S12 source tables, an explicit file-level rights allowlist and frozen hashes; it does not recompute results.
5. `scripts/reproduce_figures.py` requires the exact full hash-bound source tree and external inputs listed in metadata; those are not included. The R8.2 S12B expected path remains unresolved.
6. Model reruns require independently acquired, rights-cleared gene-level inputs. Frozen result tables alone are not model reruns. IGVF expression+LOEUF remains `NOT_RUN`/`NA`; three-platform coverage does not establish three-platform validation.
