# Runtime and dependency record

Observed validation runtime: Python 3.13.2, x86_64 under Rosetta, NumPy 2.5.1, pandas 3.0.5, matplotlib 3.11.1, Pillow 12.3.0, PyMuPDF 1.28.0, h5py 3.16.0; R 4.5.1, data.table 1.18.2.1, digest 0.6.39. These are observed versions, not a fresh-install guarantee.

On the current Apple Silicon machine, use an x86_64 compatibility interpreter because the installed compiled numerical packages are x86_64. For example:

```sh
arch -x86_64 python3 -m venv .venv
arch -x86_64 .venv/bin/python -m pip install -r environments/requirements.txt
```

Do not use native arm64 Python with these x86_64 extensions. The R scripts require `Rscript`; install the two R packages from `r-packages.txt` in a project/user library. This file does not install packages or change a global environment.

A fresh isolated package installation was not performed. Matplotlib may need a writable cache directory; set `MPLCONFIGDIR` for a single command if needed. Fonts are not bundled; renderer output may differ if DejaVu Sans is unavailable.
