#!/usr/bin/env python3
"""Build the locked R1.3 Activity NTC pseudobulk baseline, streaming sparse X.

This script writes a project-annotated expression table and model input. Those
outputs are not cleared for redistribution. See docs/SOURCE_ACQUISITION_AND_RECONSTRUCTION.md.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import subprocess
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

EXPECTED_H5AD_SHA256 = "d25a014c52cc27104fff4ad79b2efdb99130da671eda08c11f604c454cbbf23c"
EXPECTED_MEMBERSHIP_SHA256 = "3568f1e018c38e6d9400184cbffc337668e364349bddeb07e87f82039b7ebaf7"
EXPECTED_ANNOTATIONS_SHA256 = "7d0fa5f89cd9b89511c95191ec3cbf0f74d4a434180f8d81b62109aebf021134"
NTC_GUIDES = ["non-targeting_00777", "non-targeting_00873", "non-targeting_00876", "non-targeting_00908", "non-targeting_00982"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def decode(x) -> str:
    return x.decode("utf-8", "replace") if isinstance(x, (bytes, np.bytes_)) else str(x)


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, na_rep="NA", lineterminator="\n")


def qsummary(values: np.ndarray) -> tuple[float, float, float]:
    return (float(np.median(values)), float(np.quantile(values, .25)), float(np.quantile(values, .75)))


def main() -> None:
    parser = argparse.ArgumentParser(description="Rebuild the frozen R1.3 Activity NTC baseline from exact source files.")
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--membership", required=True, type=Path)
    parser.add_argument("--annotations", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    H5AD = args.h5ad.resolve(strict=True)
    MEMBERSHIP = args.membership.resolve(strict=True)
    ANNOTATIONS = args.annotations.resolve(strict=True)
    OUT = args.outdir.resolve()
    for path, expected, label in [
        (H5AD, EXPECTED_H5AD_SHA256, "H5AD"),
        (MEMBERSHIP, EXPECTED_MEMBERSHIP_SHA256, "membership"),
        (ANNOTATIONS, EXPECTED_ANNOTATIONS_SHA256, "annotations"),
    ]:
        actual = sha256(path)
        if actual != expected:
            raise RuntimeError(f"Frozen {label} SHA-256 mismatch: expected {expected}; actual {actual}")
    qdir = OUT / "01_activity_ntc_qc"
    bdir = OUT / "02_activity_baseline"
    qdir.mkdir(parents=True, exist_ok=True)
    bdir.mkdir(parents=True, exist_ok=True)
    membership = pd.read_csv(MEMBERSHIP, sep="\t", dtype={"gene": str})
    annotation = pd.read_csv(ANNOTATIONS, sep="\t", dtype={"gene": str})
    annotation = annotation.drop_duplicates("gene", keep=False).set_index("gene", drop=False)
    if len(membership) != 1271 or membership.gene.duplicated().any():
        raise RuntimeError("Frozen Activity LOEUF membership is not unique N=1,271")
    if not set(membership.gene).issubset(annotation.index):
        raise RuntimeError("Frozen membership does not fully join the frozen annotation table")

    with h5py.File(H5AD, "r") as f:
        if "X" not in f or not isinstance(f["X"], h5py.Group) or not {"data", "indices", "indptr"}.issubset(f["X"].keys()):
            raise RuntimeError("H5AD X is not the expected CSR sparse matrix")
        x = f["X"]
        data, indices, indptr = x["data"], x["indices"], x["indptr"]
        data_dtype = str(data.dtype)
        n_stored_values = len(data)
        n_obs, n_vars = int(x.attrs["shape"][0]), int(x.attrs["shape"][1])
        if n_obs != 1_514_144 or n_vars != 36_543:
            raise RuntimeError(f"H5AD shape mismatch: {n_obs} x {n_vars}")
        obs_index = decode(f["obs"].attrs.get("_index", "_index"))
        guide_group = f["obs"]["guide"]
        guide_codes = np.asarray(guide_group["codes"][:], dtype=np.int32)
        guide_categories = [decode(z) for z in guide_group["categories"][:]]
        cat_to_code = {z: i for i, z in enumerate(guide_categories)}
        if not set(NTC_GUIDES).issubset(cat_to_code):
            raise RuntimeError("One or more five locked NTC guide categories are missing")
        ntc_codes = {cat_to_code[g]: g for g in NTC_GUIDES}
        matched_var_index = decode(f["var"].attrs.get("_index", "_index"))
        var_ids = [decode(z) for z in f["var"][matched_var_index][:]]
        if len(var_ids) != n_vars or len(set(var_ids)) != n_vars:
            raise RuntimeError("H5AD var identifiers are not unique or do not match X columns")
        var_lookup = {name: i for i, name in enumerate(var_ids)}
        exp_names = list(membership.gene)
        exact_map = {g: var_lookup[g] for g in exp_names if g in var_lookup}
        mt_mask = np.fromiter((v.startswith("MT-") for v in var_ids), dtype=bool, count=n_vars)
        ribo_mask = np.fromiter((v.startswith("RPL") or v.startswith("RPS") for v in var_ids), dtype=bool, count=n_vars)

        # Verify every stored sparse value once, in bounded chunks. The full matrix is never densified.
        noninteger = 0
        min_value, max_value = np.inf, -np.inf
        for start in range(0, len(data), 5_000_000):
            values = np.asarray(data[start:min(start + 5_000_000, len(data))], dtype=np.float64)
            if values.size:
                min_value = min(min_value, float(values.min()))
                max_value = max(max_value, float(values.max()))
                noninteger += int(np.count_nonzero(np.abs(values - np.rint(values)) > 1e-6))
        if noninteger:
            raise RuntimeError(f"X integer-valued check contradicted R1.2: {noninteger} noninteger values")

        guide_counts = {g: int(np.count_nonzero(guide_codes == cat_to_code[g])) for g in NTC_GUIDES}
        assigned_rows: list[int] = []
        assigned_groups: list[str] = []
        assigned_totals: list[float] = []
        assigned_detected: list[int] = []
        assigned_mt_counts: list[float] = []
        assigned_ribo_counts: list[float] = []
        ntc_rows: list[int] = []
        ntc_counts: list[float] = []
        ntc_detected: list[int] = []
        ntc_mt_counts: list[float] = []
        ntc_ribo_counts: list[float] = []
        ntc_group_counts: dict[str, int] = {g: 0 for g in NTC_GUIDES}
        pb_all = np.zeros(n_vars, dtype=np.float64)
        pb_qc = np.zeros(n_vars, dtype=np.float64)
        block_rows = 50_000
        for row0 in range(0, n_obs, block_rows):
            row1 = min(row0 + block_rows, n_obs)
            cats = guide_codes[row0:row1]
            selected = np.flatnonzero(cats >= 0)
            if selected.size == 0:
                continue
            ptr0 = np.asarray(indptr[row0:row1 + 1], dtype=np.int64)
            val0, val1 = int(ptr0[0]), int(ptr0[-1])
            vals = np.asarray(data[val0:val1], dtype=np.float64)
            cols = np.asarray(indices[val0:val1], dtype=np.int32)
            offsets = ptr0 - val0
            pref = np.empty(len(vals) + 1, dtype=np.float64)
            pref[0] = 0
            np.cumsum(vals, out=pref[1:])
            mt_pref = np.empty(len(vals) + 1, dtype=np.float64)
            mt_pref[0] = 0
            np.cumsum(vals * mt_mask[cols], out=mt_pref[1:])
            ribo_pref = np.empty(len(vals) + 1, dtype=np.float64)
            ribo_pref[0] = 0
            np.cumsum(vals * ribo_mask[cols], out=ribo_pref[1:])
            local_rows = selected
            local_codes = cats[local_rows]
            starts, ends = offsets[local_rows], offsets[local_rows + 1]
            totals = pref[ends] - pref[starts]
            mt_totals = mt_pref[ends] - mt_pref[starts]
            ribo_totals = ribo_pref[ends] - ribo_pref[starts]
            detected = np.diff(ptr0)[local_rows]
            for j, local_row in enumerate(local_rows):
                code = int(local_codes[j])
                group = ntc_codes.get(code)
                global_row = row0 + int(local_row)
                assigned_rows.append(global_row)
                assigned_groups.append(group if group else decode(guide_categories[code]))
                assigned_totals.append(float(totals[j]))
                assigned_detected.append(int(detected[j]))
                assigned_mt_counts.append(float(mt_totals[j]))
                assigned_ribo_counts.append(float(ribo_totals[j]))
                if group:
                    ntc_rows.append(global_row)
                    ntc_counts.append(float(totals[j]))
                    ntc_detected.append(int(detected[j]))
                    ntc_mt_counts.append(float(mt_totals[j]))
                    ntc_ribo_counts.append(float(ribo_totals[j]))
                    ntc_group_counts[group] += 1
                    # Also retain all NTC and QC-passing NTC profile sums while sparse blocks are in memory.
                    lo, hi = int(ptr0[local_row] - val0), int(ptr0[local_row + 1] - val0)
                    np.add.at(pb_all, cols[lo:hi], vals[lo:hi])
                    if totals[j] >= 2000 and detected[j] >= 100:
                        np.add.at(pb_qc, cols[lo:hi], vals[lo:hi])

        assigned = pd.DataFrame({
            "guide": assigned_groups,
            "total_counts": assigned_totals,
            "detected_genes": assigned_detected,
            "identified_MT_labeled_feature_counts": assigned_mt_counts,
            "identified_RPL_RPS_labeled_feature_counts": assigned_ribo_counts,
        })
        ntc = assigned[assigned.guide.isin(NTC_GUIDES)].copy()
        if len(assigned) != 17211 or len(ntc) != 847:
            raise RuntimeError(f"Guide-assigned or NTC counts differ from R1.2: {len(assigned)}, {len(ntc)}")
        ntc["qc_pass_reconstructable"] = (ntc.total_counts >= 2000) & (ntc.detected_genes >= 100)
        n_pass = int(ntc.qc_pass_reconstructable.sum())
        if n_pass == 0:
            raise RuntimeError("No NTC cells pass the reconstructable source QC")
        # Select guide-category cell sets already retained in pb_qc; count consistency is checked below.
        guide_rows = []
        for g in NTC_GUIDES:
            z = ntc.loc[ntc.guide == g]
            guide_rows.append({
                "guide_category": g,
                "n_before_qc": len(z),
                "n_pass_total_counts_ge_2000_and_detected_genes_ge_100": int(z.qc_pass_reconstructable.sum()),
                "n_excluded_reconstructable_qc": int((~z.qc_pass_reconstructable).sum()),
                "median_total_counts": float(z.total_counts.median()),
                "median_detected_genes": float(z.detected_genes.median()),
            })
        if sum(x["n_before_qc"] for x in guide_rows) != 847 or sum(x["n_pass_total_counts_ge_2000_and_detected_genes_ge_100"] for x in guide_rows) != n_pass:
            raise RuntimeError("Per-guide QC totals do not reconcile")

    # QC summaries are aggregate only; individual barcode-level records are not exported.
    summaries = []
    for label, frame in [("all_guide_assigned", assigned), ("five_explicit_NTC_guides", ntc)]:
        summary = {
            "cell_set": label,
            "n_cells_before_qc": len(frame),
            "n_total_counts_lt_2000": int((frame.total_counts < 2000).sum()),
            "n_detected_genes_lt_100": int((frame.detected_genes < 100).sum()),
            "n_passing_both_reconstructable_thresholds": int(((frame.total_counts >= 2000) & (frame.detected_genes >= 100)).sum()),
            "total_counts_median": float(frame.total_counts.median()),
            "total_counts_q25": float(frame.total_counts.quantile(.25)),
            "total_counts_q75": float(frame.total_counts.quantile(.75)),
            "detected_genes_median": float(frame.detected_genes.median()),
            "detected_genes_q25": float(frame.detected_genes.quantile(.25)),
            "detected_genes_q75": float(frame.detected_genes.quantile(.75)),
            "mitochondrial_fraction": "NA: mixed ENSG/symbol var axis; exact MT-labeled numerator only; full reconstruction ambiguous",
            "ribosomal_fraction": "NA: mixed ENSG/symbol var axis; exact RPL/RPS-labeled numerator only; full reconstruction ambiguous",
            "mean_identified_MT_labeled_fraction_of_total": float((frame.identified_MT_labeled_feature_counts / frame.total_counts.replace(0, np.nan)).mean()),
            "mean_identified_RPL_RPS_labeled_fraction_of_total": float((frame.identified_RPL_RPS_labeled_feature_counts / frame.total_counts.replace(0, np.nan)).mean()),
        }
        summaries.append(summary)
    write_tsv(pd.DataFrame(guide_rows), qdir / "NTC_GUIDE_COUNTS.tsv")
    write_tsv(pd.DataFrame(summaries), qdir / "NTC_CELL_QC_SUMMARY.tsv")

    pd.DataFrame([{
        "input": str(H5AD), "matrix": "/X", "representation": "CSR sparse", "dtype": data_dtype,
        "shape": f"{n_obs} x {n_vars}", "stored_values_n": n_stored_values, "min_stored_value": min_value,
        "max_stored_value": max_value, "noninteger_stored_values_n": noninteger,
        "decision": "PASS_COUNT_LIKE_OPERATIONAL_PSEUDOBULK",
        "qualification_note": "H5AD has no embedded normalization provenance; stored X values are integer/count-like. Following the taskbook, X is operationally treated as counts for sum pseudobulk.",
    }]).to_csv(qdir / "EXPRESSION_SCALE_DECISION.tsv", sep="\t", index=False, na_rep="NA", lineterminator="\n")
    with (qdir / "EXPRESSION_SCALE_DECISION.md").open("w") as out:
        out.write("# Expression scale decision\n\n")
        out.write("`/X` is CSR sparse and was rechecked across all stored values. It contains only integer-valued positive entries (minimum %.0f, maximum %.0f; noninteger entries %d/%d). The H5AD itself has no embedded normalization provenance. Following the R1.3 taskbook, the matrix is operationally treated as counts for sum pseudobulk; this is a representation decision supported by observed values, not a claim about undocumented upstream processing.\n" % (min_value, max_value, noninteger, n_stored_values))

    pd.DataFrame([{
        "sample": "GSM8787485", "series": "GSE289235", "sample_title": "iPSC-Neuron, Day 21, CROP-seq, rep1",
        "cell_line": "WTC11", "sample_level_cell_type": "iPSC-neurons", "neuronal_identity_basis": "GEO sample-level experimental annotation",
        "identity_gate": "PASS", "treatment": "sgRNA library", "baseline_wording": "closely matched same-study Day-21 WTC11 iPSC-neuron CROP-seq NTC reference; exact CaMPARI2 screening batch not established",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8787485",
    }]).to_csv(OUT / "00_provenance/GEO_SAMPLE_PROVENANCE.tsv", sep="\t", index=False, lineterminator="\n")

    # Expression is derived only from the QC-passing NTC subset; all five guide groups are non-targeting.
    pb_counts = pb_qc
    libsize = float(pb_counts.sum())
    if libsize <= 0:
        raise RuntimeError("NTC pseudobulk library size is zero")
    cpm = pb_counts / libsize * 1_000_000.0
    logcpm = np.log1p(cpm)
    all_pb_cpm = pb_all / float(pb_all.sum()) * 1_000_000.0
    all_pb_logcpm = np.log1p(all_pb_cpm)
    expr_rows = []
    for _, row in membership.iterrows():
        gene = row.gene
        h5_index = var_lookup.get(gene)
        if h5_index is None:
            pc, cv, lv = np.nan, np.nan, np.nan
            matched = "NO_EXACT_MATCH"
        else:
            pc, cv, lv = float(pb_counts[h5_index]), float(cpm[h5_index]), float(logcpm[h5_index])
            matched = "EXACT_STRING_MATCH"
        ann = annotation.loc[gene]
        group = "prey" if int(row.exposure) == 1 else ("bait" if int(ann.bait_status) == 1 else "non_PPI_nonbait")
        expr_rows.append({
            "frozen_activity_gene_identifier": gene, "h5ad_identifier": gene if h5_index is not None else "NA",
            "exact_match_status": matched, "pseudobulk_count": pc, "CPM": cv, "log1p_CPM": lv,
            "exposure": int(row.exposure), "group_label": group, "hit_outcome": int(ann.hit),
            "loeuf_available": bool(np.isfinite(float(ann.loeuf_z))), "loeuf_z": float(ann.loeuf_z),
        })
    expr = pd.DataFrame(expr_rows)
    cc = expr.exact_match_status.eq("EXACT_STRING_MATCH") & expr.log1p_CPM.notna() & expr.loeuf_available
    n_cc = int(cc.sum())
    if n_cc == 0:
        raise RuntimeError("No Activity complete cases after expression join")
    expr["baseline_expression_z"] = np.nan
    expr.loc[cc, "baseline_expression_z"] = (expr.loc[cc, "log1p_CPM"] - expr.loc[cc, "log1p_CPM"].mean()) / expr.loc[cc, "log1p_CPM"].std(ddof=1)
    write_tsv(expr, bdir / "ACTIVITY_BASELINE_EXPRESSION.tsv")

    group_counts = expr.loc[cc].groupby("group_label").size().to_dict()
    coverage = pd.DataFrame([
        {"metric": "frozen_LOEUF_complete_case_N", "value": len(membership), "definition": "R1.1 frozen ACTIVITY_LOEUF_COMPLETE_CASE population"},
        {"metric": "direct_exact_match_N", "value": int(expr.exact_match_status.eq("EXACT_STRING_MATCH").sum()), "definition": "H5AD var/_index exact string equality only"},
        {"metric": "unmatched_genes_N", "value": int(expr.exact_match_status.ne("EXACT_STRING_MATCH").sum()), "definition": "no exact match; left unmatched, no alias conversion"},
        {"metric": "direct_exact_match_coverage_fraction", "value": float(expr.exact_match_status.eq("EXACT_STRING_MATCH").mean()), "definition": "direct exact matches / 1,271 frozen members"},
        {"metric": "final_expression_complete_case_N", "value": n_cc, "definition": "exact expression match plus finite frozen LOEUF"},
        {"metric": "prey_N_final_complete_cases", "value": int(group_counts.get("prey", 0)), "definition": "exposure=1"},
        {"metric": "non_PPI_nonbait_N_final_complete_cases", "value": int(group_counts.get("non_PPI_nonbait", 0)), "definition": "exposure=0 and bait_status=0"},
        {"metric": "bait_N_final_complete_cases", "value": int(group_counts.get("bait", 0)), "definition": "exposure=0 and bait_status=1"},
        {"metric": "NTC_source_QC_passing_cells", "value": n_pass, "definition": "only five explicit NTC categories; counts>=2000 and detected genes>=100"},
        {"metric": "pseudobulk_library_size", "value": int(round(libsize)), "definition": "sum of all features over qualified NTC cells"},
        {"metric": "zero_pseudobulk_count_exact_matched_genes_N", "value": int((expr.loc[expr.exact_match_status.eq("EXACT_STRING_MATCH"), "pseudobulk_count"] == 0).sum()), "definition": "zero is measured and retained, not missing"},
        {"metric": "zero_pseudobulk_count_exact_matched_genes_fraction", "value": float((expr.loc[expr.exact_match_status.eq("EXACT_STRING_MATCH"), "pseudobulk_count"] == 0).mean()), "definition": "fraction of exact-matched frozen Activity genes with measured zero pseudobulk count"},
    ])
    write_tsv(coverage, bdir / "ACTIVITY_BASELINE_COVERAGE_QA.tsv")

    descriptives = []
    matched = expr.loc[cc].copy()
    for group in ["prey", "non_PPI_nonbait", "bait"]:
        vals = matched.loc[matched.group_label == group, "log1p_CPM"].to_numpy(float)
        median, q25, q75 = qsummary(vals)
        descriptives.append({"group_label": group, "n": len(vals), "median_log1p_CPM": median, "q25": q25, "q75": q75, "zero_pseudobulk_count_fraction": float((matched.loc[matched.group_label == group, "pseudobulk_count"] == 0).mean())})
    rho = float(pd.Series(matched.baseline_expression_z).corr(pd.Series(matched.loeuf_z), method="spearman"))
    profile_idx = np.array([var_lookup[g] for g in exp_names if g in var_lookup], dtype=int)
    profile_corr = float(pd.Series(np.log1p(pb_all[profile_idx] / pb_all.sum() * 1e6)).corr(pd.Series(np.log1p(pb_qc[profile_idx] / pb_qc.sum() * 1e6)), method="spearman"))
    if n_pass == 847:
        profile_status = "COMPUTED; all-NTC and QC-passing NTC are the same 847 cells"
    else:
        profile_status = "COMPUTED across exact-matched Activity genes"
    write_tsv(pd.DataFrame(descriptives), bdir / "ACTIVITY_BASELINE_GROUP_DESCRIPTIVES.tsv")
    pd.DataFrame([{
        "n_final_complete_cases": n_cc, "spearman_baseline_expression_vs_loeuf": rho,
        "spearman_all_NTC_vs_source_QC_NTC_expression_profile": profile_corr, "profile_comparison_status": profile_status,
    }]).to_csv(bdir / "ACTIVITY_BASELINE_CORRELATION_QA.tsv", sep="\t", index=False, lineterminator="\n")

    # Frozen model input in the exact R1.1 membership order.
    model = annotation.loc[membership.gene].reset_index(drop=True).copy()
    # Remove the R1.1 placeholder NA field before joining the measured R1.3 expression.
    model = model.drop(columns=[c for c in ["baseline_expression_z", "baseline_expression_status"] if c in model.columns])
    model = model.merge(expr[["frozen_activity_gene_identifier", "h5ad_identifier", "exact_match_status", "pseudobulk_count", "CPM", "log1p_CPM", "baseline_expression_z"]], left_on="gene", right_on="frozen_activity_gene_identifier", how="left", sort=False, validate="one_to_one")
    model = model.drop(columns=["frozen_activity_gene_identifier"])
    model["model_row_order"] = np.arange(1, len(model) + 1)
    model["complete_case"] = model.baseline_expression_z.notna() & model.loeuf_z.notna()
    model = model.loc[model.complete_case].copy().reset_index(drop=True)
    model["model_row_order"] = np.arange(1, len(model) + 1)
    model["group_label"] = np.where(model.exposure.astype(int).eq(1), "prey", np.where(model.bait_status.astype(int).eq(1), "bait", "non_PPI_nonbait"))
    write_tsv(model, OUT / "03_activity_primary/ACTIVITY_COMMON_COMPLETE_CASE_MODEL_INPUT.tsv")
    write_tsv(model[["model_row_order", "gene"]], OUT / "03_activity_primary/ACTIVITY_MODEL_POPULATION_ORDER.tsv")
    if len(model) != n_cc or int(model.exposure.sum()) != int(group_counts.get("prey", 0)):
        raise RuntimeError("Activity model input population does not reconcile with expression coverage")

    # Portable file provenance; all input bytes were hashed before this build began.
    input_manifest = pd.DataFrame([
        {"logical_role": "GEO_GSM8787485_Day21_WTC11_iPSC_neuron_NTC_candidate", "path": H5AD.name, "bytes": H5AD.stat().st_size, "sha256": EXPECTED_H5AD_SHA256, "sha_basis": "fresh SHA256 verified against frozen R1.2 value", "sample_id": "EXP224_10X"},
        {"logical_role": "frozen_R1.1_Activity_annotations", "path": ANNOTATIONS.name, "bytes": ANNOTATIONS.stat().st_size, "sha256": EXPECTED_ANNOTATIONS_SHA256, "sha_basis": "fresh SHA256 verified against frozen R1.1 value"},
        {"logical_role": "frozen_R1.1_Activity_LOEUF_complete_case_membership", "path": MEMBERSHIP.name, "bytes": MEMBERSHIP.stat().st_size, "sha256": EXPECTED_MEMBERSHIP_SHA256, "sha_basis": "fresh SHA256 verified against frozen R1.1 value"},
    ])
    write_tsv(input_manifest, OUT / "00_provenance/INPUT_MANIFEST_R1_3.tsv")

    with (qdir / "NTC_CELL_QC_DECISION.md").open("w") as out:
        out.write("# NTC cell QC decision\n\n")
        out.write("Neuronal identity passes on GEO sample-level experimental annotation: GSM8787485 is a Day-21 WTC11 iPSC-neuron CROP-seq sample. The absence of a redundant per-cell neuron label does not fail the sample identity gate.\n\n")
        out.write(f"The H5AD has {n_obs:,} observations and {n_vars:,} features; {len(assigned):,} cells have an assigned guide. The five locked explicit non-targeting categories contain {len(ntc):,} observations. The source count/gene filters (total counts >=2,000 and detected genes >=100) are reconstructable from sparse X; {n_pass:,} NTC observations pass both. Primary NTC pseudobulk uses only those {n_pass:,} cells.\n\n")
        out.write("GEO describes additional ribosomal (3–8%) and mitochondrial (<6.5%) QC. The H5AD feature axis mixes Ensembl IDs and symbol-like identifiers and supplies no separate symbol annotation, so full mitochondrial/ribosomal fractions cannot be unambiguously reconstructed. Exact `MT-`, `RPL`, and `RPS` labeled-feature counts were calculated as lower-bound diagnostics but were not used to exclude cells.\n\n")
        out.write("All five explicit NTC categories were included; no perturbed-guide cells entered either NTC profile. See `NTC_GUIDE_COUNTS.tsv` and aggregate `NTC_CELL_QC_SUMMARY.tsv`; cell-level barcodes and values are not exported.\n")

    print(json.dumps({"status": "BASELINE_BUILT", "n_obs": n_obs, "n_vars": n_vars, "guide_assigned": len(assigned), "ntc_before": len(ntc), "ntc_after": n_pass, "noninteger_x": noninteger, "activity_membership_n": len(membership), "exact_matches": len(exact_map), "model_n": len(model), "prey_n": int(model.exposure.astype(int).sum()), "pseudobulk_library_size": int(round(libsize)), "output": str(OUT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
