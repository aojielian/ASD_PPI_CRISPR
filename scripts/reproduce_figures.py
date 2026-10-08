#!/usr/bin/env python3
"""Run the frozen figure renderer from an explicit data root."""
import argparse, subprocess, sys
from pathlib import Path
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--data-dir",required=True,type=Path);ap.add_argument("--outdir",required=True,type=Path);a=ap.parse_args()
 code=Path(__file__).resolve().parents[1]/"src/plotting/r5_render_figures.py"
 raise SystemExit(subprocess.call([sys.executable,str(code),"--input",str(a.data_dir.resolve()),"--output",str(a.outdir.resolve())]))
if __name__=="__main__":main()
