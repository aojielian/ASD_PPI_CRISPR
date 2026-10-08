#!/usr/bin/env python3
"""Recompute frozen in-vivo group means and condition directions after exclusion.

The eight supplied condition columns are already frozen normal scores. This
script never reranks or renormalizes genes. It does not perform resampling.
"""
import argparse, csv, re
from pathlib import Path

CONDITIONS = ["Z_Syn1-GW-2M", "Z_Emx1-GW-2M", "Z_Vglut2-GW-2M", "Z_Vgat-GW-2M",
              "Z_Syn1-GW-4M", "Z_Emx1-GW-4M", "Z_Vglut2-GW-4M", "Z_Vgat-GW-4M"]

def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--scores", required=True, type=Path)
    ap.add_argument("--depmap-list", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    a=ap.parse_args()
    with a.depmap_list.open(encoding="utf-8-sig", newline="") as f:
        positive=set()
        for r in csv.DictReader(f):
            label=(r.get("Essentials") or "").strip()
            m=re.fullmatch(r"(.+?)\s*\((\d+)\)",label)
            if not m: raise ValueError("DepMap list has an unexpected row")
            positive.add(m.group(1).strip())
    rows=read_tsv(a.scores)
    required={"human_gene_symbol_ensembl","primary_group",*CONDITIONS,"composite_vulnerability_8screen"}
    if not rows or required-set(rows[0]): raise ValueError("Input schema mismatch")
    original=rows
    kept=[r for r in rows if r["human_gene_symbol_ensembl"] not in positive]
    if len(original)!=9949 or len(kept)!=8772: raise ValueError("Frozen input/exclusion membership counts differ")
    out=[]
    for label, data in [("original",original),("depmap_excluded",kept)]:
        for col in ["composite_vulnerability_8screen",*CONDITIONS]:
            prey=[float(r[col]) for r in data if r["primary_group"]=="prey_only"]
            comp=[float(r[col]) for r in data if r["primary_group"]=="nonppi_nonbait"]
            delta=sum(prey)/len(prey)-sum(comp)/len(comp)
            out.append({"population":label,"score":col,"N":len(data),"prey_N":len(prey),"comparator_N":len(comp),
                        "prey_mean":sum(prey)/len(prey),"comparator_mean":sum(comp)/len(comp),"delta":delta,
                        "direction":"positive" if delta>0 else "negative" if delta<0 else "zero",
                        "scale_note":"frozen original condition ranks/normal scores; no reranking"})
    a.outdir.mkdir(parents=True,exist_ok=True)
    with (a.outdir/"INVIVO_DETERMINISTIC_GROUP_SUMMARIES.tsv").open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(out)
if __name__=="__main__": main()
