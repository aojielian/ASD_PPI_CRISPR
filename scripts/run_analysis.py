#!/usr/bin/env python3
"""Run one frozen analysis module using explicit paths only."""
import argparse, shutil, subprocess, sys
from pathlib import Path

def main():
 ap=argparse.ArgumentParser();ap.add_argument("module",choices=["activity","igvf","invivo"]);ap.add_argument("--data-dir",required=True,type=Path);ap.add_argument("--outdir",required=True,type=Path);a=ap.parse_args()
 data=a.data_dir.resolve(); out=a.outdir.resolve(); out.mkdir(parents=True,exist_ok=True)
 if a.module=="activity":
  inp=data/"input/activity/ACTIVITY_COMMON_COMPLETE_CASE_MODEL_INPUT.tsv"; exe=shutil.which("Rscript")
  script=Path(__file__).resolve().parents[1]/"src/activity/fit_activity_expression_loeuf.R"
  if not inp.is_file():raise SystemExit("Blocked: license-cleared Activity complete-case model input is not in the current data package.")
  if not exe:raise SystemExit("Rscript is required (see environments/README.md).")
  raise SystemExit(subprocess.call([exe,str(script),str(inp),str(out)]))
 if a.module=="igvf":
  one=data/"input/igvf/IGVF_623.tsv";two=data/"input/igvf/IGVF_610_LOEUF.tsv";exe=shutil.which("Rscript")
  if not one.is_file() or not two.is_file():raise SystemExit("Blocked: IGVF model inputs include BioPlex-derived covariates whose redistribution permission is unresolved.")
  if not exe:raise SystemExit("Rscript is required (see environments/README.md).")
  script=Path(__file__).resolve().parents[1]/"src/igvf/run_igvf_models.R"
  raise SystemExit(subprocess.call([exe,str(script),str(one),str(two),str(out)]))
 scores=data/"input/invivo/INVIVO_COMPOSITE_SCORES.tsv"; depmap=data/"input/common_essential_exclusion/CRISPRInferredCommonEssentials.csv"
 if not scores.is_file():raise SystemExit("Blocked: the frozen in-vivo score input is withheld pending source-file redistribution review.")
 script=Path(__file__).resolve().parents[1]/"src/invivo/common_essential_descriptive.py"
 raise SystemExit(subprocess.call([sys.executable,str(script),"--scores",str(scores),"--depmap-list",str(depmap),"--outdir",str(out)]))
if __name__=="__main__":main()
