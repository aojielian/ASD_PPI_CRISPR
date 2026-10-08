#!/usr/bin/env python3
"""COPY_ONLY: copy supplied clean frozen S11/S12 TSVs from a hash-verified rights allowlist.

This is a file-copy utility. It does not rebuild tables or recalculate models.
"""
import argparse, csv, hashlib, shutil
from pathlib import Path, PurePosixPath

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser(description="COPY_ONLY: copies existing clean S11/S12 TSVs; does not rebuild or analyze them.")
    ap.add_argument("--source-root",required=True,type=Path)
    ap.add_argument("--outdir",required=True,type=Path)
    ap.add_argument("--rights-manifest",required=True,type=Path,help="TSV with relative_path,sha256,decision; all files must be PUBLISH")
    a=ap.parse_args(); src=a.source_root.resolve(strict=True); out=a.outdir.resolve()
    with a.rights_manifest.open(encoding="utf-8",newline="") as f:
        allow={r["relative_path"]:r for r in csv.DictReader(f,delimiter="\t")}
    copied=0
    for section in ("S11","S12"):
        folder=src/f"{section}_source_tables"
        files=sorted(folder.glob("*.tsv"))
        if not files: raise SystemExit(f"Missing clean {section} source TSVs")
        for p in files:
            rel=p.relative_to(src).as_posix(); row=allow.get(rel)
            if row is None or row.get("decision")!="PUBLISH":
                raise SystemExit(f"Blocked pending explicit file-level rights clearance: {rel}")
            if sha(p)!=row["sha256"]: raise SystemExit(f"Frozen source hash mismatch: {rel}")
            dest=out/section/p.name; dest.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(p,dest); copied+=1
    print(f"Copied {copied} hash-verified rights-cleared S11/S12 TSVs; no tables or analyses were rebuilt.")
if __name__=="__main__":main()
