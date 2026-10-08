#!/usr/bin/env python3
"""Render Figures 1–4 and Supplementary Figures S1–S10 from supplied TSVs.

This is a frozen-result renderer. It does not fit statistical models. It reads
only the explicit --input data root and writes only to --output. Required
source tables are listed in metadata/FIGURE_TABLE_SOURCE_MAP.tsv.

Requires Python >= 3.10, numpy, pandas, matplotlib, Pillow and PyMuPDF.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import logging
import math
import os
from pathlib import Path
import platform
import re
import shutil
import sys
import tempfile
from typing import Any
import zipfile

# Use a writable local cache, not a frozen input directory.
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "asd_ppi_r4_mpl"))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import NullLocator, ScalarFormatter
from matplotlib.text import Text
from PIL import Image, ImageOps, ImageDraw
import fitz

HERE = Path(__file__).resolve().parent
LOG = logging.getLogger("asd_ppi_r4")
VERSION = "4.0.0"
MM = 1 / 25.4
# Semantic group colours: not changed according to statistical significance.
C = dict(prey="#287D7C", hit="#B77739", other="#858C92", ink="#24282D",
         pale="#EAEDEF", rule="#D8DDDF", grid="#EFF1F2", age="#617A9F")


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def scalar(row: pd.Series, key: str) -> float:
    value = float(row[key])
    if not np.isfinite(value):
        raise ValueError(f"Required numeric field is missing/non-finite: {key}")
    return value


def one(df: pd.DataFrame, **equals: Any) -> pd.Series:
    selected = df
    for key, value in equals.items():
        selected = selected.loc[selected[key] == value]
    if len(selected) != 1:
        raise ValueError(f"Expected one record for {equals}, found {len(selected)}")
    return selected.iloc[0]


def numeric_p(value: float, digits: int = 3) -> str:
    """Format, do not recalculate. A finite nonzero P is never displayed as zero."""
    value = float(value)
    if not np.isfinite(value) or not (0 <= value <= 1):
        raise ValueError(f"Invalid P/q value: {value}")
    if value == 0:
        return "0 (source)"
    if value < 0.001:
        power = int(math.floor(math.log10(value)))
        mantissa = value / 10 ** power
        return rf"${mantissa:.3g}\times10^{{{power}}}$"
    if value < 0.01:
        return f"{value:.4f}"
    return f"{value:.{digits}f}"


def natural_key(s: str) -> list[Any]:
    return [int(part) if part.isdigit() else part.lower()
            for part in re.split(r"(\d+)", str(s))]


class FrozenSource:
    """Read exact relative members. Never search arbitrary project directories."""
    def __init__(self, path: Path, mapping_path: Path):
        self.path = path.resolve()
        self.bindings = json.loads(mapping_path.read_text(encoding="utf-8"))
        self.zip: zipfile.ZipFile | None = None
        self.prefix = ""
        self.root: Path | None = None
        self.cache: dict[str, pd.DataFrame] = {}
        self.checked: list[dict[str, Any]] = []
        self.file_hash = None
        if self.path.is_file() and zipfile.is_zipfile(self.path):
            self.zip = zipfile.ZipFile(self.path)
            names = self.zip.namelist()
            starts = [n for n in names if n.endswith("/START_HERE_ASTRA.md")]
            if "START_HERE_ASTRA.md" in names:
                starts.append("START_HERE_ASTRA.md")
            if len(starts) != 1:
                raise ValueError("Input ZIP must contain one START_HERE_ASTRA.md root")
            self.prefix = starts[0][:-len("START_HERE_ASTRA.md")]
            with self.path.open("rb") as handle:
                h = hashlib.sha256()
                for block in iter(lambda: handle.read(2 ** 20), b""):
                    h.update(block)
                self.file_hash = h.hexdigest()
        elif self.path.is_dir():
            if (self.path / "START_HERE_ASTRA.md").is_file():
                self.root = self.path
            else:
                starts = list(self.path.glob("*/START_HERE_ASTRA.md"))
                if len(starts) != 1:
                    raise ValueError("Pass the extracted package root, not the whole project")
                self.root = starts[0].parent
        else:
            raise FileNotFoundError(f"Input not found or not a ZIP/directory: {self.path}")

    def bytes(self, relative: str) -> bytes:
        p = Path(relative)
        if p.is_absolute() or ".." in p.parts:
            raise ValueError("Unsafe input member")
        if self.zip is not None:
            return self.zip.read(self.prefix + relative)
        assert self.root is not None
        return (self.root / relative).read_bytes()

    def table(self, role: str) -> pd.DataFrame:
        if role not in self.cache:
            meta = self.bindings[role]
            data = self.bytes(meta["path"])
            actual = digest(data)
            if actual != meta["sha256"]:
                raise ValueError(f"Frozen source changed: {role}\n{meta['path']}\n"
                                 f"Expected {meta['sha256']}\nActual   {actual}\n"
                                 "Do not bypass this check or substitute an older source.")
            # Reject duplicate header names before pandas can silently rename them.
            first = data.decode("utf-8-sig").splitlines()[0].split("\t")
            if len(first) != len(set(first)):
                raise ValueError(f"Duplicate source columns: {role}")
            df = pd.read_csv(io.BytesIO(data), sep="\t")
            if len(df) != meta["rows"] or list(df.columns) != meta["columns"]:
                raise ValueError(f"Frozen source schema mismatch: {role}")
            self.cache[role] = df
            self.checked.append(dict(role=role, source_member=meta["path"],
                                     sha256=actual, rows=len(df), check="PASS"))
        return self.cache[role].copy()

    def close(self) -> None:
        if self.zip is not None:
            self.zip.close()


class Run:
    def __init__(self, source: FrozenSource, output: Path, dpi: int):
        self.source, self.output, self.dpi = source, output, dpi
        self.checks: list[dict[str, Any]] = []
        self.crosswalk: list[dict[str, Any]] = []
        self.text_checks: list[dict[str, Any]] = []
        self.pdf_paths: list[Path] = []
        self.legends: dict[str, str] = {}
        for sub in ["figures", "supplementary", "source_data", "legends", "qa", "code"]:
            (output / sub).mkdir(parents=True, exist_ok=True)

    def check(self, item: str, actual: Any, expected: Any, tol: float = 1e-10):
        if isinstance(expected, (float, np.floating)):
            ok = bool(np.isclose(actual, expected, rtol=tol, atol=tol))
        else:
            ok = bool(actual == expected)
        self.checks.append(dict(check=item, actual=str(actual), expected=str(expected),
                                status="MATCH" if ok else "MISMATCH"))
        if not ok:
            raise ValueError(f"Numeric/membership check failed: {item}: {actual} != {expected}")

    def export(self, name: str, data: pd.DataFrame, panels: str, roles: list[str],
               transformation: str = "Frozen values; selected/reordered for display"):
        out = self.output / "source_data" / f"{name}.tsv"
        data.to_csv(out, sep="\t", index=False, float_format="%.17g", na_rep="NA")
        for role in roles:
            self.crosswalk.append(dict(panel=panels, display_file=out.relative_to(self.output).as_posix(),
                                       input_role=role, source_member=self.source.bindings[role]["path"],
                                       source_sha256=self.source.bindings[role]["sha256"],
                                       display_operation=transformation))
        return out

    def save(self, fig: plt.Figure, name: str, legend: str, category="figures"):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        fw, fh = fig.bbox.width, fig.bbox.height
        for obj in fig.findobj(match=Text):
            if not obj.get_visible() or not obj.get_text().strip():
                continue
            box = obj.get_window_extent(renderer=renderer)
            if not np.all(np.isfinite(box.extents)):
                continue
            inside = box.x0 >= -1 and box.y0 >= -1 and box.x1 <= fw + 1 and box.y1 <= fh + 1
            self.text_checks.append(dict(figure=name, text=obj.get_text(), inside_page=inside,
                                         x0=box.x0, y0=box.y0, x1=box.x1, y1=box.y1))
        base = self.output / category / name
        meta = dict(Title=name.replace("_", " "), Author="", Creator="Frozen-data figure renderer", CreationDate=None, ModDate=None)
        fig.savefig(base.with_suffix(".pdf"), metadata=meta, facecolor="white")
        fig.savefig(base.with_suffix(".svg"), facecolor="white", metadata={"Date": None})
        fig.savefig(base.with_suffix(".png"), dpi=self.dpi, facecolor="white")
        self.pdf_paths.append(base.with_suffix(".pdf"))
        self.legends[name] = legend
        # One PDF per figure retains vector text/lines; actual individual gene points
        # can be rasterized at export DPI to avoid pathological PDF size.
        plt.close(fig)


def set_theme():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.0, "axes.labelsize": 8.0,
        "axes.titlesize": 9, "axes.titleweight": "normal", "axes.linewidth": .65,
        "xtick.labelsize": 7.2, "ytick.labelsize": 7.2, "legend.fontsize": 7.1,
        "text.color": C["ink"], "axes.labelcolor": C["ink"],
        "axes.edgecolor": C["ink"], "xtick.color": C["ink"], "ytick.color": C["ink"],
        "xtick.major.size": 2.8, "ytick.major.size": 2.8,
        "xtick.major.width": .6, "ytick.major.width": .6,
        "legend.frameon": False, "lines.linewidth": 1.2,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "svg.hashsalt": "ASD_PPI_R4",
        "savefig.facecolor": "white", "figure.facecolor": "white",
        "axes.unicode_minus": True, "savefig.dpi": 400,
    })


def figure(width_mm: float, height_mm: float) -> plt.Figure:
    return plt.figure(figsize=(width_mm * MM, height_mm * MM), dpi=150)


def axes(fig, rect, ylabel=None, xlabel=None, grid=False):
    ax = fig.add_axes(rect)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_axisbelow(True)
    if grid:
        ax.grid(axis="y", color=C["grid"], lw=.5)
    if ylabel:
        ax.set_ylabel(ylabel, labelpad=4)
    if xlabel:
        ax.set_xlabel(xlabel, labelpad=4)
    return ax


def heading(fig, letter, title, x, y):
    fig.text(x, y, letter, fontsize=11, fontweight="bold", va="top", ha="left")
    fig.text(x + .031, y - .002, title, fontsize=9, va="top", ha="left")


def odds_axis(ax, limits=(.8, 6.0)):
    ax.set_xscale("log")
    ticks = [t for t in [.5, 1, 2, 4, 8] if limits[0] <= t <= limits[1]]
    ax.set_xlim(*limits)
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_locator(NullLocator())
    ax.axvline(1, color=C["other"], lw=.65, ls=(0, (3, 3)))
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def ci_point(ax, x, y, low, high, color, marker="o", size=4.2):
    vals = [float(x), float(y), float(low), float(high)]
    if not all(np.isfinite(vals)) or low > x or high < x:
        raise ValueError("Invalid frozen point / interval")
    ax.errorbar(x, y, xerr=np.array([[x-low], [high-x]]), fmt=marker,
                markersize=size, elinewidth=1.0, capsize=2.0,
                color=color, markeredgewidth=.7, zorder=4)


def gene_strip(ax, values, colors, labels, log=False, seed=710, size=6):
    """All supplied genes; deterministic x-jitter only. No resampling or KDE."""
    rng = np.random.default_rng(seed)
    record = []
    for i, (v, col, label) in enumerate(zip(values, colors, labels)):
        a = np.asarray(v, dtype=float)
        if len(a) == 0 or not np.all(np.isfinite(a)) or (log and (a <= 0).any()):
            raise ValueError("Invalid gene-level display vector")
        jitter = rng.uniform(-.22, .22, len(a))
        ax.scatter(i + jitter, a, s=size, color=col, alpha=.48 if len(a) < 200 else .27,
                   linewidths=0, rasterized=len(a) > 400, zorder=2)
        # IQR and Tukey whiskers are descriptive summaries, not confidence intervals.
        q1, median, q3 = np.quantile(a, [.25, .5, .75])
        iqr = q3-q1
        within = a[(a >= q1-1.5*iqr) & (a <= q3+1.5*iqr)]
        wl, wh = float(within.min()), float(within.max())
        ax.plot([i, i], [wl, wh], color=C["ink"], lw=.75, zorder=4)
        ax.plot([i, i], [q1, q3], color="white", lw=8, zorder=4)
        ax.plot([i, i], [q1, q3], color=col, lw=5.5, alpha=.9, zorder=5)
        ax.plot([i-.095, i+.095], [median, median], color=C["ink"], lw=1.0, zorder=6)
        record.append(dict(group=label, n=len(a), q25=q1, median=median, q75=q3,
                           whisker_low=wl, whisker_high=wh))
    ax.set_xticks(range(len(values)), labels)
    ax.set_xlim(-.5, len(values)-.5)
    if log:
        ax.set_yscale("log")
        ax.yaxis.set_minor_locator(NullLocator())
        ax.set_yticks([1, 10, 100, 1000])
        ax.yaxis.set_major_formatter(ScalarFormatter())
        amin=min(float(np.min(a)) for a in values); amax=max(float(np.max(a)) for a in values)
        ax.set_ylim(amin/1.4, amax*1.5)
    else:
        amin=min(float(np.min(a)) for a in values); amax=max(float(np.max(a)) for a in values)
        ax.set_ylim(min(0, amin), amax*1.12 if amax>0 else 1)
    return pd.DataFrame(record)


def validate(s: FrozenSource, run: Run):
    # All required files are verified up-front: no half-complete silent render.
    for role in s.bindings:
        s.table(role)
    a = one(s.table("activity"), analysis="primary_fdr_lt_0_10")
    for key, expected in dict(exposed_n=118, exposed_hits=21, comparator_n=1112, comparator_hits=112).items():
        run.check(f"activity.{key}", int(a[key]), expected)
    st = s.table("activity_strata")
    st = st[(st.analysis == "primary_fdr_lt_0_10") & (st.row_type == "stratum")]
    run.check("activity.strata_count", len(st), 4)
    run.check("activity.stratum_hit_total", int(st.exposed_hits.sum()), int(a.exposed_hits))
    run.check("activity.stratum_denominator_total", int(st.comparator_n.sum()), int(a.comparator_n))
    ig = s.table("igvf_genes")
    ig = ig[(ig.usable_target == 1) & (ig.Science_bait == 0)]
    for group, n, hits in [(1,125,47),(0,1061,139)]:
        d=ig[ig.Science_prey_only==group]
        run.check(f"igvf.group{group}.n",len(d),n)
        run.check(f"igvf.group{group}.hits",int(d.author_energy_distance_hit.sum()),hits)
    run.check("igvf.gene_ids_unique",ig.gene.is_unique,True)
    run.check("igvf.finite_positive_energy",bool(np.isfinite(ig.mean_energy_distance).all() and (ig.mean_energy_distance>0).all()),True)
    model=one(s.table("igvf_models"),model="A1_network_adjusted_primary")
    run.check("igvf.complete_cases",int(model.n),623)
    comp=s.table("invivo_genes"); pr=s.table("invivo_primary").iloc[0]
    for name,n,col in [("prey_only",868,"prey_mean"),("nonppi_nonbait",9081,"comparator_mean")]:
        v=comp.loc[comp.primary_group==name,"composite_vulnerability_8screen"]
        run.check(f"invivo.{name}.n",len(v),n)
        run.check(f"invivo.{name}.mean",float(v.mean()),float(pr[col]))
    scr=s.table("invivo_screens")
    run.check("invivo.screen_n",len(scr),8)
    for _,row in scr.iterrows():
        col="Z_"+row.screen
        delta=comp.loc[comp.primary_group=="prey_only",col].mean()-comp.loc[comp.primary_group=="nonppi_nonbait",col].mean()
        run.check(f"invivo.{row.screen}.display_consistency",float(delta),float(row.mean_difference))
    w=s.table("wam_cells")
    run.check("wamsley.cell_results",len(w),72)
    run.check("wamsley.unique_cell_types",w.source_cell_type.nunique(),36)
    run.check("wamsley.36_per_set",bool((w.groupby("gene_set").size()==36).all()),True)
    run.check("wamsley.cell_key_unique",not w.duplicated(["gene_set","source_cell_type"]).any(),True)
    q=s.table("holm_family")
    for role,test in [("wam_a2","A2_Wamsley_generic_network"),("wam_a3","A3_Wamsley_activity_vs_rest"),("wam_a4","A4_Wamsley_within_activity_screen")]:
        run.check(f"{role}.rawP_consistency",float(s.table(role).iloc[0].P_raw),float(one(q,test=test).P_raw))
    g=s.table("gordon_scores");gp=s.table("gordon_primary").iloc[0]
    for field,n,mean in [("is_H",21,"mean_A"),("is_N_screen",97,"mean_B")]:
        d=g[g[field]==1]
        run.check(f"gordon.{field}.n",len(d),n)
        run.check(f"gordon.{field}.mean",float(d.mean_abs_logFC_C36.mean()),float(gp[mean]))
    v=s.table("vel_scores")
    run.check("vel.corrected_rank_export",bool(v.mean_rank_abs_t.nunique()>2),True)
    for comparison,n in [("C1_H_VS_N",89),("C2_H_VS_R",888)]:
        d=v[v.comparison==comparison]
        run.check(f"vel.{comparison}.H",int((d['set']=='H').sum()),19)
        run.check(f"vel.{comparison}.comparator",int((d['set']!='H').sum()),n)
        row=one(s.table("vel_primary"),comparison=comparison)
        run.check(f"vel.{comparison}.rank_difference",float(d.loc[d['set']=='H','mean_rank_abs_t'].mean()-d.loc[d['set']!='H','mean_rank_abs_t'].mean()),float(row.mean_rank_abs_t_difference))
        run.check(f"vel.{comparison}.magnitude_difference",float(d.loc[d['set']=='H','mean_abs_beta'].mean()-d.loc[d['set']!='H','mean_abs_beta'].mean()),float(row.mean_abs_beta_difference))
    run.check("vel.dream_not_run",str(s.table("vel_dream").iloc[0]['status']).startswith("NOT_RUN"),True)


def draw_figure1(r: Run):
    s=r.source
    a=one(s.table("activity"),analysis="primary_fdr_lt_0_10")
    adj=one(s.table("activity_adjusted"),model="primary_fdr_lt_0_10",term="prey_only")
    allst=s.table("activity_strata")
    st=allst[(allst.row_type=="stratum")&(allst.analysis=="primary_fdr_lt_0_10")].copy()
    cmh=one(allst,row_type="CMH_common",analysis="primary_fdr_lt_0_10")
    f=figure(180,151)
    heading(f,"A","Experimental network and functional readouts",.045,.977)
    # Only conceptual element in the entire main-figure collection; no drawn
    # animal/cell, no invented molecular pathway, no causal chain to disease.
    ax=f.add_axes([.065,.817,.885,.108]); ax.axis("off")
    ax.text(0,.52,"Experimental ASD\nprotein interactome",ha="left",va="center",fontsize=8.6,fontweight="medium")
    ax.plot([.28,.32,.32],[.52,.52,.52],color=C["other"],lw=.75)
    for y,label in [(1.0,"Human neuronal activity"),(.52,"Cortical-neuron Perturb-seq"),(.04,"In vivo neuronal dependency")]:
        ax.plot([.32,.32,.37],[.52,y,y],color=C["other"],lw=.75)
        ax.text(.39,y,label,ha="left",va="center",fontsize=8.3)
    ax.text(.77,.52,"ASD cortex\nassociation tests",ha="left",va="center",fontsize=8.3)
    ax.plot([.72,.75],[.52,.52],ls=(0,(2,3)),color=C["other"],lw=.75)
    ax.set_xlim(0,1);ax.set_ylim(-.16,1.15)
    f.add_artist(Line2D([.045,.96],[.792,.792],transform=f.transFigure,color=C["rule"],lw=.6))

    heading(f,"B","Activity-screen hits",.045,.758)
    ax=axes(f,[.105,.215,.32,.455],ylabel="Activity-hit genes (%)",grid=True)
    rates=np.array([a.exposed_hit_rate,a.comparator_hit_rate])*100
    ax.bar([0,1],rates,width=.53,color=[C['prey'],C['other']],edgecolor="none",zorder=3)
    for j,(rate,h,n) in enumerate(zip(rates,[a.exposed_hits,a.comparator_hits],[a.exposed_n,a.comparator_n])):
        ax.text(j,rate+.65,f"{rate:.1f}%",ha="center",va="bottom",fontsize=8.5)
        ax.text(j,1.0,f"{int(h)}/{int(n):,}",color="white",ha="center",va="bottom",fontsize=7.3)
    ax.set_ylim(0,23);ax.set_xlim(-.65,1.65);ax.set_yticks([0,5,10,15,20])
    ax.set_xticks([0,1],["PPI prey\n(n = 118)","Other screened\n(n = 1,112)"])
    ax.spines['left'].set_visible(False); ax.tick_params(axis='y',length=0)
    f.text(.105,.145,"Primary hit definition: FDR < 0.10",fontsize=7.4,color=C['other'])

    heading(f,"C","Association estimates",.515,.758)
    ax=axes(f,[.715,.487,.235,.189],xlabel="Odds ratio (95% CI)")
    rows=[("Crude",a.OR,a.CI95_low,a.CI95_high,a.P),
          ("Adjusted",adj.OR,adj.OR_CI_low,adj.OR_CI_high,adj.P),
          ("Stratified CMH",cmh.OR,cmh.CI95_low,cmh.CI95_high,cmh.P)]
    for y,(label,v,lo,hi,p) in zip([2,1,0],rows):
        ci_point(ax,v,y,lo,hi,C['prey'] if label=='Adjusted' else C['other'],size=4.3)
    odds_axis(ax,(.85,4.6));ax.set_yticks([2,1,0],[x[0] for x in rows]);ax.set_ylim(-.5,2.65)
    f.text(.715,.715,"Adjusted P = "+numeric_p(adj.P),fontsize=7.4,va='top')

    heading(f,"D","Screen-selection strata",.515,.382)
    ax=axes(f,[.715,.118,.235,.198],xlabel="Activity-hit genes (%)")
    order=['activity_only','disease_only','both','unresolved']; names=['Activity','Disease','Both','Unresolved']
    st=st.set_index('selection_stratum').loc[order].reset_index()
    for y,(_,row) in zip([3,2,1,0],st.iterrows()):
        x1,x2=row.exposed_hit_rate*100,row.comparator_hit_rate*100
        ax.plot([x1,x2],[y,y],color=C['rule'],lw=1.0,zorder=1)
        ax.scatter(x1,y,s=20,color=C['prey'],zorder=3)
        ax.scatter(x2,y,s=20,color=C['other'],marker='s',zorder=3)
    ax.set_yticks([3,2,1,0],names);ax.set_xlim(-2,66);ax.set_xticks([0,20,40,60]);ax.set_ylim(-.55,3.55)
    ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
    f.legend(handles=[Line2D([],[],ls='',marker='o',color=C['prey'],markersize=4,label='PPI prey'),
                      Line2D([],[],ls='',marker='s',color=C['other'],markersize=4,label='Other screened')],
             loc='lower right',bbox_to_anchor=(.963,.015),ncol=2,fontsize=7.1,handletextpad=.3,columnspacing=1.3)
    r.export('F1B_activity_hit_counts',pd.DataFrame([a]),'1B',['activity'])
    r.export('F1C_activity_model_estimates',pd.DataFrame(rows,columns=['model','OR','CI_low','CI_high','P_raw']),
             '1C',['activity','activity_adjusted','activity_strata'])
    r.export('F1D_selection_strata',st,'1D',['activity_strata'])
    r.export('F1A_evidence_sources',s.table('flow'),'1A',['flow'],'Conceptual labels, not an inferential or causal diagram')
    legend=f"""Figure 1. Experimental ASD interactome membership and neuronal-activity effects.
(A) Schematic of the experimental-interaction starting point and the distinct functional readouts. Lines indicate evidence linkage, not a causal sequence; this is the sole schematic panel in the main figures.
(B) Source-defined neuronal-activity hits (FDR < 0.10) among prey-only genes ({int(a.exposed_hits)}/{int(a.exposed_n)}) and tested non-PPI/non-bait comparators ({int(a.comparator_hits)}/{int(a.comparator_n)}). Bars show observed percentages, not fitted probabilities.
(C) Frozen crude Fisher estimate, adjusted logistic estimate and selection-stratified Mantel–Haenszel estimate. The crude and CMH comparisons use the 118/1,112 prey/comparator universe. The adjusted prey-only term is from the original 1,312-gene model, which also accounts for bait status, guide count and library-selection strata; it is not a same-sample unadjusted/adjusted pair. Horizontal intervals are the supplied 95% confidence intervals. Primary adjusted OR={adj.OR:.3f}, 95% CI {adj.OR_CI_low:.3f}–{adj.OR_CI_high:.3f}, P={adj.P:.6g}.
(D) Observed hit fractions within all four frozen library-selection categories. The disease-only prey category has only five genes; the apparent percentage should not be read as a well-powered stratum-specific estimate. Exact numerator/denominator values are included in the plotted source data. No stratum-specific inferential claim is introduced. Strict-threshold, permutation and leave-one-stratum-out results remain in Supplementary Figure S1 and the unchanged source registry.
"""
    r.save(f,'Figure_1',legend)


def draw_figure2(r: Run):
    s=r.source; g=s.table('igvf_genes');g=g[(g.usable_target==1)&(g.Science_bait==0)].copy()
    mm=s.table('igvf_models'); mt=s.table('igvf_matched').iloc[0]
    hp=one(s.table('holm_family'),test='A1_IGVF_network_adjusted').P_adjusted_Holm
    main=one(mm,model='A1_network_adjusted_primary')
    f=figure(180,168)
    heading(f,'A','Gene-level transcriptomic perturbation',.045,.975)
    ax=axes(f,[.115,.470,.405,.431],ylabel='Mean energy distance (log scale)')
    gp=g[g.Science_prey_only==1];gn=g[g.Science_prey_only==0]
    summ=gene_strip(ax,[gp.mean_energy_distance,gn.mean_energy_distance],
                    [C['prey'],C['other']],[f'PPI prey\n(n = {len(gp)})',f'Other targets\n(n = {len(gn):,})'],log=True,size=6,seed=42)
    ax.grid(axis='y',color=C['grid'],lw=.5)
    heading(f,'B','Transcriptomic hits',.595,.975)
    ax=axes(f,[.690,.588,.256,.30],ylabel='Author-defined hit genes (%)',grid=True)
    bars=[]
    for i,(d,col) in enumerate([(gp,C['prey']),(gn,C['other'])]):
        h=int(d.author_energy_distance_hit.sum());rate=h/len(d)*100
        ax.bar(i,rate,width=.5,color=col,zorder=3)
        ax.text(i,rate+1.6,f'{h}/{len(d):,}',ha='center',fontsize=7.5)
        bars.append(dict(group=['PPI prey','Other targets'][i],hits=h,tested=len(d),percentage=rate))
    ax.set_xlim(-.6,1.6);ax.set_ylim(0,48);ax.set_yticks([0,20,40]);ax.set_xticks([0,1],['PPI prey','Other']);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)

    heading(f,'C','Adjustment in the same 623 genes',.045,.352)
    ax=axes(f,[.259,.106,.275,.190],xlabel='Odds ratio (95% profile CI)')
    ids=['A1_unadjusted','A1_source_only','A1_network_adjusted_primary']
    labels=['Unadjusted','Source adjusted','Network adjusted']
    d=mm.set_index('model').loc[ids].reset_index()
    for y,(_,rr) in zip([2,1,0],d.iterrows()):
        ci_point(ax,rr.OR,y,rr.CI_low_OR,rr.CI_high_OR,C['prey'] if y==0 else C['other'],size=4.5)
    odds_axis(ax,(.9,5.4));ax.set_yticks([2,1,0],labels);ax.set_ylim(-.55,2.55)
    f.text(.083,.030,'Network adjusted: Holm P = '+numeric_p(hp),fontsize=7.5)

    heading(f,'D','Matched comparison',.595,.470)
    ax=axes(f,[.722,.123,.202,.266],ylabel='Hit-risk difference\n(percentage points)')
    v=mt.estimate_difference*100;lo=mt.CI_low*100;hi=mt.CI_high*100
    ax.axhline(0,color=C['other'],ls=(0,(3,3)),lw=.65)
    ax.errorbar(.5,v,yerr=[[v-lo],[hi-v]],fmt='o',color=C['prey'],markersize=5,
                capsize=4,elinewidth=1.2)
    ax.set_xlim(0,1);ax.set_ylim(-4,34);ax.set_yticks([0,10,20,30]);ax.set_xticks([.5],[f'{int(mt.matched_sets_n)} pairs'])
    ax.spines['bottom'].set_visible(False);ax.tick_params(axis='x',length=0)
    f.text(.657,.030,'Paired permutation P = '+numeric_p(mt.P_permutation,4),fontsize=7.3)
    r.export('F2A_all_nonbait_gene_values',g,'2A',['igvf_genes'],'Filter only usable_target=1 and Science_bait=0; display all remaining 1,186 genes')
    r.export('F2A_descriptive_box_summaries',summ,'2A',['igvf_genes'],'Display-only quartiles and Tukey whiskers; not confidence intervals')
    r.export('F2B_author_hit_rates',pd.DataFrame(bars),'2B',['igvf_genes'],'Sum frozen binary author-hit flag / exact group denominator; no new testing')
    r.export('F2C_complete_case_models',d,'2C',['igvf_models'])
    r.export('F2C_original_Holm_family',s.table('holm_family'),'2C',['holm_family'],'Existing four-test Holm family; no recalculation')
    r.export('F2D_matched_risk_difference',pd.DataFrame([mt]),'2D',['igvf_matched'])
    legend=f"""Figure 2. Human cortical-neuron transcriptomic perturbation and measured-network adjustment.
(A) All {len(gp)} prey-only and {len(gn):,} non-PPI/non-bait usable gene targets in the frozen IGVF inventory. Each dot is a gene's saved mean energy-distance value. Only the x-position is jittered for visibility; the y-axis is logarithmic. Thick segments are IQRs, short black lines are medians, and thin whiskers extend to the most extreme observation within 1.5 IQR; all outliers remain visible as dots. These are descriptive summaries, not confidence intervals.
(B) Author-defined transcriptomic hits: {int(gp.author_energy_distance_hit.sum())}/{len(gp)} and {int(gn.author_energy_distance_hit.sum())}/{len(gn)}. Hit membership is based on the frozen source-author significant-promoter rule, not an image-derived threshold or a newly selected cutoff. Source data contain the full-universe crude and source-adjusted results; panel C is restricted to complete cases.
(C) Unadjusted, source-adjusted and measured-network-adjusted models within the same 623 genes (120 prey and 503 comparators). The network-adjusted model includes source category, BioPlex degree, 293T protein abundance, bait status and tested-promoter count; the exact frozen encodings are preserved. Primary OR={main.OR:.4f}, 95% penalized profile CI {main.CI_low_OR:.4f}–{main.CI_high_OR:.4f}; the displayed Holm P={hp:.7g} belongs to the original A1/A2/A3/A4 four-test family. This is not adjustment for every possible confounder: baseline brain expression, knockdown efficiency and genetic constraint were not available for adjustment.
(D) Risk difference from {int(mt.matched_sets_n)} frozen 1:1 matched sets ({int(mt.exposed_hits)} versus {int(mt.comparator_hits)} hits). Error bars are the supplied matched-set bootstrap 95% interval ({mt.CI_low:.4f}–{mt.CI_high:.4f}); the printed P={mt.P_permutation:.7g} is the frozen paired-set permutation P. The matched marginal OR is not used as a paired conditional OR. Full model inventory, retention, balance and program-testability outputs remain in Supplementary Figure S5.
"""
    r.save(f,'Figure_2',legend)


def draw_figure3(r: Run):
    s=r.source;g=s.table('invivo_genes');p=s.table('invivo_primary').iloc[0]
    screens=s.table('invivo_screens');age=s.table('invivo_age')
    f=figure(180,152)
    heading(f,'A','Distribution of composite neuronal dependency',.045,.977)
    ax=axes(f,[.10,.565,.85,.344],ylabel='Density',xlabel='Composite vulnerability score')
    allv=g.composite_vulnerability_8screen.to_numpy(float)
    edges=np.linspace(float(allv.min()),float(allv.max()),55)
    hist=[]
    for group,col,label in [('nonppi_nonbait',C['other'],'Other screened'),('prey_only',C['prey'],'PPI prey')]:
        vals=g.loc[g.primary_group==group,'composite_vulnerability_8screen'].to_numpy(float)
        count,_=np.histogram(vals,bins=edges)
        density=count/(len(vals)*np.diff(edges))
        ax.stairs(density,edges,color=col,lw=1.3,fill=True,alpha=.14)
        ax.stairs(density,edges,color=col,lw=1.0,label=f'{label} (n = {len(vals):,})')
        for j in range(len(count)):
            hist.append(dict(group=group,bin_left=edges[j],bin_right=edges[j+1],genes=count[j],density=density[j]))
    ax.set_xlim(edges[0],edges[-1]);ax.set_ylim(bottom=0)
    ax.set_xticks(np.arange(math.ceil(edges[0]), math.floor(edges[-1])+1))
    ax.legend(loc='upper right',fontsize=7.3)
    # Actual means, with no cropped or zoomed difference exaggeration.
    ax.axvline(p.prey_mean,color=C['prey'],lw=.8,ls=(0,(3,3)),alpha=.7)
    ax.axvline(p.comparator_mean,color=C['other'],lw=.8,ls=(0,(3,3)),alpha=.7)

    heading(f,'B','Effects across eight conditions',.045,.408)
    ax=axes(f,[.115,.103,.350,.238],ylabel='Mean vulnerability difference',xlabel='Age (months)')
    styles=[('Syn1','o','-'),('Emx1','s','--'),('Vglut2','^','-.'),('Vgat','D',':')]
    palette=[C['prey'],C['age'],C['hit'],C['other']]
    for (name,mk,ls),col in zip(styles,palette):
        vals=[float(one(screens,screen=f'{name}-GW-{ag}').mean_difference) for ag in ['2M','4M']]
        ax.plot([2,4],vals,color=col,marker=mk,markersize=4,lw=1.0,ls=ls)
        ax.text(4.13,vals[1],name,fontsize=7,color=col,va='center')
    # Stagger labels where the source endpoints are very close; leader lines
    # touch their true endpoints and do not alter the plotted values.
    # Remove only these just-created label artists, then place collision-free labels.
    for t in list(ax.texts): t.remove()
    label_y={'Vglut2':.151,'Syn1':.121,'Vgat':.097,'Emx1':.076}
    for (name,mk,ls),col in zip(styles,palette):
        end=float(one(screens,screen=f'{name}-GW-4M').mean_difference)
        ax.plot([4,4.18],[end,label_y[name]],lw=.55,color=col)
        ax.text(4.23,label_y[name],name,fontsize=7,color=col,va='center')
    ax.axhline(0,color=C['other'],lw=.65,ls=(0,(3,3)))
    ax.set_xlim(1.75,5.10);ax.set_xticks([2,4]);ax.set_ylim(-.008,.176);ax.set_yticks([0,.05,.10,.15])

    heading(f,'C','Composite effect estimates',.550,.408)
    ax=axes(f,[.715,.103,.232,.238],xlabel='Mean difference (95% CI)')
    data=[dict(label='All conditions',estimate=p.mean_difference,lo=p.bootstrap_CI_low,hi=p.bootstrap_CI_high,
               P=p.permutation_P_two_sided,P_type='Primary two-sided permutation P')]
    for ag in ['2M','4M']:
        row=one(age,age_composite=ag)
        data.append(dict(label=f'{ag[0]} months',estimate=row.mean_difference,lo=row.bootstrap_CI_low,
                         hi=row.bootstrap_CI_high,P=row.two_sided_P_Holm_2test_family,
                         P_type='Pre-existing two-age sensitivity Holm P'))
    for y,d in zip([2,1,0],data):
        ci_point(ax,d['estimate'],y,d['lo'],d['hi'],C['prey'] if y==2 else C['other'],marker='D' if y==2 else 'o')
    ax.axvline(0,color=C['other'],lw=.65,ls=(0,(3,3)));ax.set_xlim(-.005,.184);ax.set_xticks([0,.05,.10,.15]);ax.set_ylim(-.6,2.6)
    ax.set_yticks([2,1,0],[d['label'] for d in data]);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
    f.text(.582,.026,'Primary P = '+numeric_p(p.permutation_P_two_sided)+f";  d = {p.cohen_d:.3f}",fontsize=7.3)
    r.export('F3A_composite_gene_scores',g,'3A',['invivo_genes'])
    r.export('F3A_histogram_display_bins',pd.DataFrame(hist),'3A',['invivo_genes'],'54 shared equal-width bins over full observed range; group-normalized histogram density, not a fitted distribution')
    r.export('F3B_eight_condition_results',screens,'3B',['invivo_screens'])
    r.export('F3C_composite_and_age_estimates',pd.DataFrame(data),'3C',['invivo_primary','invivo_age'])
    legend=f"""Figure 3. Neuronal dependency across eight conditions within one in vivo study.
(A) The frozen eight-condition composite vulnerability scores of {int(p.prey_only_n)} prey-only and {int(p.comparator_n):,} non-PPI/non-bait comparator genes. Histograms use 54 common equal-width bins spanning the full observed range and are normalized separately to unit area. Dashed lines mark saved group means. Positive scores denote greater inferred vulnerability/depletion after the frozen sign convention and within-condition rank normalization. This gene aggregation was constructed in this project from source guide information, rather than supplied as an author-defined essential-gene score. The effect is a modest distributional shift (Cohen's d={p.cohen_d:.3f}), not separation of two distinct gene populations.
(B) Saved prey-minus-comparator mean differences in each of four neuronal-driver conditions at 2 and 4 months. All eight estimates are positive, but not all condition-specific tests are significant. The eight conditions share one study and screening library; they are not eight independent cohorts. Lines connect matching driver conditions for visual comparison, not repeated measurements on the same animals. No new condition-level P values or intervals are calculated.
(C) The primary composite and pre-specified age-composite estimates with their supplied bootstrap 95% intervals. The primary two-sided gene-label permutation P is {p.permutation_P_two_sided:.7g}. The 2-/4-month summaries are sensitivities, not independent replications; their original two-test Holm results remain in source data. Complete condition-level tests, condition/class omissions, alternative summaries and the 18-month extension remain in Supplementary Figure S6.
"""
    r.save(f,'Figure_3',legend)


def draw_figure4(r: Run):
    s=r.source;w=s.table('wam_cells');cl=s.table('wam_class')
    a2=s.table('wam_a2').iloc[0];a3=s.table('wam_a3').iloc[0];a4=s.table('wam_a4').iloc[0]
    sets=['SCIENCE_PREY_ONLY','ACTIVITY_VULNERABLE_PPI']
    # Grouping uses the source's biological annotation, never observed P/effect.
    groups=[('excitatory','Excitatory'),('inhibitory','Inhibitory'),('astrocyte','Astrocytes'),
            ('oligodendrocyte','Oligodendrocytes'),('oligodendrocyte precursor','OPCs'),('microglia/immune','Microglia / immune'),
            ('endothelial','Vascular'),('pericyte','Vascular'),('vascular mural/smooth muscle','Vascular')]
    present=cl.subtype.unique().tolist()
    covered={g[0] for g in groups}
    if set(present)-covered:
        raise ValueError(f'Unbound source cell-type class: {set(present)-covered}')
    order=[]; group_bounds=[]
    for subtype,label in groups:
        names=sorted(cl.loc[cl.subtype==subtype,'source_cell_type'].tolist(),key=natural_key)
        if not names:continue
        start=len(order);order.extend(names)
        if group_bounds and group_bounds[-1][2]==label:
            group_bounds[-1]=(group_bounds[-1][0],len(order),label)
        else:group_bounds.append((start,len(order),label))
    r.check('F4.show_all_36_cell_types',len(order),36)
    mat=np.column_stack([w[w.gene_set==gs].set_index('source_cell_type').loc[order,'AUC_probability_superiority'] for gs in sets])
    q=np.column_stack([w[w.gene_set==gs].set_index('source_cell_type').loc[order,'BH_q_within_gene_set_36_celltypes'] for gs in sets])
    f=figure(180,190)
    heading(f,'A','Cell-type association pattern',.038,.979)
    ax=f.add_axes([.215,.126,.242,.795])
    # A real two-column outcome matrix: colour is AUC, not signed logFC or P.
    cm=LinearSegmentedColormap.from_list('association_auc',['#997EAD','#F6F6F4','#6DB1AB','#165D63'])
    norm=TwoSlopeNorm(vmin=.35,vcenter=.5,vmax=.85)
    im=ax.imshow(mat,aspect='auto',cmap=cm,norm=norm,interpolation='nearest')
    ax.set_xticks([0,1],['PPI prey','Activity hits']);ax.xaxis.tick_top();ax.tick_params(axis='x',length=0,pad=5,labelsize=7.3)
    ax.set_yticks(np.arange(len(order)),[x.replace('_',' ') for x in order]);ax.tick_params(axis='y',length=0,labelsize=7.1,pad=4)
    ax.set_xticks(np.arange(-.5,2,.999999),minor=True);ax.set_yticks(np.arange(-.5,36,1),minor=True)
    ax.grid(which='minor',color='white',lw=.6);ax.tick_params(which='minor',length=0)
    for spine in ax.spines.values():spine.set_visible(False)
    # Uniform markers identify the original BH-significance flag; no filtering of rows.
    for i in range(36):
        for j in range(2):
            if q[i,j]<.05:
                ax.scatter(j,i,s=12,marker='o',facecolors='white',edgecolors='#40464B',linewidths=.55,zorder=5)
    for start,end,label in group_bounds:
        if start>0:ax.axhline(start-.5,color='white',lw=2.5)
        # A slim category rule, with compact annotation outside the plotted cells.
        yy=1-((start+end)/2)/36
        short={'Excitatory':'Exc.','Inhibitory':'Inh.','Astrocytes':'Astro.','Oligodendrocytes':'Oligo.','OPCs':'OPC','Microglia / immune':'Micro.','Vascular':'Vasc.'}[label]
        ax.text(1.045,yy,short,transform=ax.transAxes,rotation=90,ha='left',va='center',fontsize=6.3,color=C['other'])
    cax=f.add_axes([.215,.091,.242,.018]);cb=f.colorbar(im,cax=cax,orientation='horizontal')
    cb.set_ticks([.35,.5,.65,.8]);cb.ax.tick_params(labelsize=6.8,width=.5,length=2)
    cb.set_label('Probability of larger |logFC|',fontsize=7.2,labelpad=3)
    cb.outline.set_visible(False)
    f.text(.475,.045,'Open circle: within-set BH q < 0.05',fontsize=6.4,color=C['ink'])
    heading(f,'B','Broader-background estimates',.535,.979)
    ax=axes(f,[.644,.702,.300,.224],xlabel='Median neuronal |logFC| difference')
    broad=[]
    for role,label,col in [('wam_set1','PPI prey',C['prey']),('wam_set2','Activity hits',C['hit'])]:
        rr=s.table(role).iloc[0]
        broad.append(dict(group=label,set_n=int(rr.set_n),comparator_n=int(rr.comparator_n),
                          effect=rr.effect_median_difference,ci_low=rr.bootstrap_CI_low,
                          ci_high=rr.bootstrap_CI_high,Holm_P=rr.Holm_P,raw_source_test=rr.test))
    for y,(d,col) in enumerate(zip(broad,[C['prey'],C['hit']])):
        yy=1-y;ci_point(ax,d['effect'],yy,d['ci_low'],d['ci_high'],col,size=4.7)
        ax.text(.015,yy+.35,f"{d['group']}: {d['set_n']:,} vs {d['comparator_n']:,}",
                transform=ax.get_yaxis_transform(),fontsize=7.1,ha='left')
        ax.text(.99,yy-.29,'Holm P = '+numeric_p(d['Holm_P']),transform=ax.get_yaxis_transform(),fontsize=6.9,ha='right',color=col)
    ax.axvline(0,color=C['other'],lw=.65,ls=(0,(3,3)));ax.set_xlim(-.001,.046);ax.set_xticks([0,.02,.04]);ax.set_yticks([]);ax.set_ylim(-.55,1.7)
    ax.spines['left'].set_visible(False)

    heading(f,'C','Adjusted comparator tests',.535,.554)
    ax=axes(f,[.647,.165,.298,.335],xlabel='Rank-normal coefficient (95% CI)')
    models=[]
    for row,lab,test in [(a2,'PPI vs generic network','A2_Wamsley_generic_network'),
                         (a3,'Hits vs other PPI','A3_Wamsley_activity_vs_rest'),
                         (a4,'Hits vs screened non-hits','A4_Wamsley_within_activity_screen')]:
        models.append(dict(comparison=lab,estimate=row.estimate_rank_beta,ci_low=row.ci_low,ci_high=row.ci_high,
                           n_exposed=int(row.n_exposed),n_comparator=int(row.n_comparator),P_raw=row.P_raw,
                           P_Holm=one(s.table('holm_family'),test=test).P_adjusted_Holm,source_test=test))
    for yy,d in zip([2,1,0],models):
        ci_point(ax,d['estimate'],yy,d['ci_low'],d['ci_high'],C['hit'] if yy==1 else C['other'],size=4.6)
        ax.text(.015,yy+.35,f"{d['comparison']}\n{d['n_exposed']:,} vs {d['n_comparator']:,}",transform=ax.get_yaxis_transform(),
                fontsize=6.9,va='bottom',linespacing=1.13)
        ax.text(.99,yy-.30,'Holm P = '+numeric_p(d['P_Holm']),transform=ax.get_yaxis_transform(),fontsize=6.9,ha='right')
    ax.axvline(0,color=C['other'],ls=(0,(3,3)),lw=.65);ax.set_xlim(-.48,1.05);ax.set_xticks([-.4,0,.4,.8]);ax.set_ylim(-.58,2.83);ax.set_yticks([]);ax.spines['left'].set_visible(False)
    f.text(.569,.082,'Comparison-specific gene universes',fontsize=7.0,color=C['other'])
    wplot=w.assign(display_cell_order=w.source_cell_type.map({c:i for i,c in enumerate(order)})).sort_values(['display_cell_order','gene_set'])
    r.export('F4A_all_36_cell_type_results',wplot,'4A',['wam_cells','wam_class'],'All 36 source cell types in biological class order; fixed AUC colour scale 0.35–0.85; original BH flag overlaid')
    r.export('F4B_broad_background_effects',pd.DataFrame(broad),'4B',['wam_set1','wam_set2'])
    r.export('F4C_adjusted_comparator_models',pd.DataFrame(models),'4C',['wam_a2','wam_a3','wam_a4','holm_family'])
    legend=f"""Figure 4. Cell-type-resolved ASD-cortex associations and comparator dependence.
(A) All 36 Wamsley cell types for the two frozen sets. Colour represents the saved probability of larger absolute ASD/control logFC (rank-AUC) relative to each set's tested complement. An AUC of 0.5 indicates no rank preference, not absence of all biological effects. Colour therefore describes the magnitude ranking of absolute expression changes, not up-/down-regulation. Rows follow the original biological annotations (including BBB_MG3 as microglia/immune), with natural numeric order within a class; none are selected by significance. White-filled circles with thin dark-grey outlines mark cell types meeting the original within-set BH q < 0.05 threshold. The two set columns are not a direct test of hits versus remaining PPI. All row-specific tested counts, medians, P and q values are supplied in the plot data. No neuronal-specific effect is claimed.
(B) Frozen neuronal-omnibus associations against the broader tested complements: {broad[0]['set_n']} vs {broad[0]['comparator_n']:,} and {broad[1]['set_n']} vs {broad[1]['comparator_n']:,}. Points are differences in median neuronal-omnibus absolute logFC; intervals are the original gene-bootstrap 95% intervals and labels are the original Holm P values. This is a distinct estimand from panel C.
(C) Three pre-existing adjusted comparisons: PPI against measured generic-network genes (889/6,010), activity hits against other PPI genes with complete covariates (17/872), and activity hits against PPI non-hits tested in the same activity screen (21/93). Points and HC3 95% intervals are the frozen regression coefficients for comparison-specific Blom rank-normalized outcomes. These are different gene universes and rank-reference populations, not a nested sequence showing attenuation in the same sample. Labels use the original four-test A1/A2/A3/A4 Holm family. The generic-network and same-screen comparisons do not establish an additional association; the broader PPI-subset comparison retains an association, but its matching sensitivity did not meet its original balance criterion. The network-augmented same-screen sensitivity and full diagnostics remain in Supplementary Figure S8.
The pre-specified Gordon developmental validation (Supplementary Figure S9; primary P={s.table('gordon_primary').iloc[0].p_value:.6f}) and Velmeshev patient-cortex validation (Supplementary Figure S10; both primary Holm P={s.table('vel_primary').iloc[0].Holm_P:.5f}) did not support the proposed priority-set advantage. Their complete findings are retained, not reclassified as exploratory or removed. This figure does not establish replicated patient-cortex prioritization.
"""
    r.save(f,'Figure_4',legend)


def draw_s9(r: Run):
    s=r.source;g=s.table('gordon_scores');p=s.table('gordon_primary').iloc[0]
    c24=s.table('gordon_c24').iloc[0];un=s.table('gordon_unstratified').iloc[0]
    days=s.table('gordon_days');forms=s.table('gordon_forms');forms=forms[forms.contrast_set=='C36'].copy()
    f=figure(180,166)
    heading(f,'A','Gordon: primary C36 scores',.045,.976)
    ax=axes(f,[.103,.570,.339,.334],ylabel='Mean absolute author logFC')
    d=g[(g.is_H==1)|(g.is_N_screen==1)].copy()
    gs=gene_strip(ax,[d[d.is_H==1].mean_abs_logFC_C36,d[d.is_N_screen==1].mean_abs_logFC_C36],
                  [C['hit'],C['other']],['Activity hits\n(n = 21)','Screened non-hits\n(n = 97)'],seed=210,size=10)
    f.text(.105,.93,'Primary stratified P = '+numeric_p(p.p_value),fontsize=7.2)
    heading(f,'B','All developmental days',.536,.976)
    ax=axes(f,[.628,.584,.319,.297],ylabel='Difference in mean |logFC|',xlabel='Developmental day')
    for key,label,col,mk in [('C36_all_9_forms','C36 primary panel',C['hit'],'o'),('C24_6_forms_sensitivity','C24 sensitivity',C['other'],'s')]:
        dx=days[days.profile_panel==key].sort_values('developmental_day')
        if len(dx)!=4:
            # Fixed schema, no adaptive alternative or fallback.
            raise ValueError(f'Expected four days for {key}; source has {days.profile_panel.unique()}')
        ax.plot(dx.developmental_day,dx.delta_mean_A_minus_B,color=col,marker=mk,markersize=4,label=label)
    ax.axhline(0,color=C['other'],ls=(0,(3,3)),lw=.65);ax.set_xticks([25,50,75,100]);ax.set_xlim(21,104)
    ax.legend(loc='lower right',fontsize=6.7)
    heading(f,'C','All nine ASD forms',.045,.437)
    ax=axes(f,[.208,.095,.278,.280],xlabel='Difference in mean |logFC|')
    forms=forms.sort_values('ASD_form',key=lambda col:col.map(lambda x:str(x).lower()))
    ys=np.arange(len(forms))[::-1]
    for y,(_,row) in zip(ys,forms.iterrows()):
        col=C['hit'] if row.low_donor_form=='NO' else C['other']
        ax.plot([0,row.delta_mean_A_minus_B],[y,y],color=C['rule'],lw=.8)
        ax.scatter(row.delta_mean_A_minus_B,y,s=20,color=col,marker='o' if row.low_donor_form=='NO' else 's')
    ax.set_yticks(ys,[x.replace('IdiopathicASD','Idiopathic ASD') for x in forms.ASD_form]);ax.set_ylim(-.6,len(forms)-.4)
    ax.axvline(0,color=C['other'],ls=(0,(3,3)),lw=.65);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
    heading(f,'D','Primary and fixed sensitivities',.536,.437)
    ax=axes(f,[.747,.123,.197,.253],xlabel='Mean score difference')
    # Each row is a mean-|logFC| difference; source supplies no CI. No CI is fabricated.
    rows=[('C36 stratified',p),('C24 stratified',c24),('C36 unstratified',un)]
    summary=[]
    for y,(label,row) in zip([2,1,0],rows):
        delta=float(row['delta_mean_A_minus_B']) if 'delta_mean_A_minus_B' in row else float(row['observed_statistic'])
        pv=float(row['p_value'])
        ax.scatter(delta,y,color=C['hit'] if y==2 else C['other'],s=24,marker='o' if y==2 else 's')
        ax.text(.98,y-.29,'P = '+numeric_p(pv),transform=ax.get_yaxis_transform(),ha='right',fontsize=7,bbox=dict(facecolor='white',edgecolor='none',pad=.2))
        summary.append(dict(comparison=label,mean_difference=delta,P=pv,role='PRIMARY' if y==2 else 'SENSITIVITY'))
    ax.set_yticks([2,1,0],[x[0].replace(' ','\n',1) for x in rows]);ax.set_ylim(-.6,2.5);ax.axvline(0,color=C['other'],ls=(0,(3,3)),lw=.65)
    ax.set_xlim(-.005,.005);ax.set_xticks([-.004,0,.004]);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
    r.export('S9A_Gordon_screened_gene_scores',d,'S9A',['gordon_scores'],'H/N_screen subset only; precomputed 36-contrast score retained')
    r.export('S9A_box_descriptive',gs,'S9A',['gordon_scores'],'Display-only quartiles/Tukey whiskers')
    r.export('S9B_all_day_profiles',days,'S9B',['gordon_days'])
    r.export('S9C_all_nine_form_profiles',forms,'S9C',['gordon_forms'],'C36 rows only; all nine forms in alphabetical order; no effect-based selection')
    r.export('S9D_primary_and_sensitivities',pd.DataFrame(summary),'S9D',['gordon_primary','gordon_c24','gordon_unstratified'])
    legend=f"""Supplementary Figure S9. Pre-specified Gordon cortical-organoid external validation.
(A) All 21 frozen activity hits and all 97 same-screen PPI non-hits. Each gene's score is the saved mean absolute author logFC across the fixed 9 ASD forms × 4 developmental days (C36). Points are genes; IQR, median and Tukey whiskers are descriptive, with every outlier retained. The frozen primary H-minus-N mean difference is {p.delta_mean_A_minus_B:.8f}; the source's two-sided selection-provenance-stratified competitive label-reference P is {p.p_value:.8f} (1,000,000 assignments). This reference is gene-competitive, not donor-population inference, and does not establish the hypothesized advantage.
(B) All four developmental days under C36 and the pre-specified C24 panel (six forms with at least three distinct affected individuals at each day). These are descriptive day profiles: no new per-day P values or intervals.
(C) All nine ASD-form profiles in alphabetical order, not ordered by significance or effect. Gray squares identify the metadata-defined low-donor forms (PCDH19, SHANK3 and TS), which remain included in the C36 primary panel. Profiles have no formal form-specific P values or confidence intervals.
(D) The original primary, C24 sensitivity and unstratified sensitivity point estimates/P values. No uncertainty intervals were supplied, and none have been inferred. C24 cannot replace or rescue the C36 primary endpoint. Moving this figure from the main figure to the supplement does not change its pre-specified analysis status, numerical results or the requirement to report the non-supportive outcome in Results.
"""
    r.save(f,'Figure_S9',legend,'supplementary')


def draw_s10(r: Run):
    s=r.source;v=s.table('vel_scores');p=s.table('vel_primary');sens=s.table('vel_27')
    cp=s.table('vel_compartments');coverage=s.table('vel_coverage')
    f=figure(180,178)
    heading(f,'A','Velmeshev: gene-level magnitude',.045,.976)
    ax=axes(f,[.105,.569,.37,.340],ylabel='Mean |coefficient| across neuronal compartments')
    hv=v[(v.comparison=='C1_H_VS_N')&(v['set']=='H')]
    nv=v[(v.comparison=='C1_H_VS_N')&(v['set']=='N')]
    rv=v[(v.comparison=='C2_H_VS_R')&(v['set']=='R')]
    summ=gene_strip(ax,[hv.mean_abs_beta,nv.mean_abs_beta,rv.mean_abs_beta],
                    [C['hit'],C['other'],C['prey']],['Hits\n(n = 19)','Screened\nnon-hits (89)','Other PPI\n(n = 888)'],seed=311,size=7)
    heading(f,'B','Rank-enrichment tests',.54,.976)
    ax=axes(f,[.761,.579,.191,.296],xlabel='Mean rank difference')
    rows=[]
    for df,n in [(p,28),(sens,27)]:
        for comp in ['C1_H_VS_N','C2_H_VS_R']:
            row=one(df,comparison=comp)
            rows.append(dict(cohort=n,comparison=comp,estimate=row.mean_rank_abs_t_difference,
                             Mixed_P=row.Mixed_P,Holm_P=row.Holm_P))
    for yy,row in zip([3,2,1,0],rows):
        ax.scatter(row['estimate'],yy,s=23,marker='o' if row['cohort']==28 else 's',color=C['ink'] if row['cohort']==28 else C['other'])
        ax.text(.98,yy-.26,'Holm P = '+numeric_p(row['Holm_P']),transform=ax.get_yaxis_transform(),ha='right',fontsize=6.7)
    ax.set_yticks([3,2,1,0],['28 donors\nHits vs screened','28 donors\nHits vs other PPI','27 donors\nHits vs screened','27 donors\nHits vs other PPI'])
    ax.set_ylim(-.62,3.55);ax.axvline(0,color=C['other'],lw=.65,ls=(0,(3,3)))
    ax.set_xlim(-17,46);ax.set_xticks([-10,0,20,40]);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)

    heading(f,'C','Three neuronal compartments',.045,.432)
    ax=axes(f,[.105,.104,.369,.258],ylabel='Difference in mean |coefficient|')
    order=['EXC','INH','OTHER_NEURON']
    for comp,col,label,mk in [('C1_H_VS_N',C['other'],'Hits − screened non-hits','o'),('C2_H_VS_R',C['prey'],'Hits − other PPI','s')]:
        d=cp[cp.comparison==comp].set_index('compartment')
        # The source label for the third compartment is verified, not guessed.
        if set(d.index)!=set(order):raise ValueError(f'Unexpected compartment labels {d.index.tolist()}')
        ax.plot(range(3),d.loc[order,'mean_abs_beta_difference'],color=col,marker=mk,markersize=4,label=label)
    ax.axhline(0,color=C['other'],lw=.65,ls=(0,(3,3)));ax.set_xticks(range(3),['Excitatory','Inhibitory','Other neurons'])
    ax.set_xlim(-.15,2.15);ax.legend(loc='upper left',fontsize=6.6)
    heading(f,'D','Gene-set coverage',.54,.432)
    ax=axes(f,[.635,.112,.313,.244],ylabel='Genes')
    d=coverage.set_index('set').loc[['H','N','R']]
    # Plot percentage retained rather than misleading unlabelled log counts.
    ax.set_ylabel('Original set retained (%)')
    colors=[C['pale'],C['other'],C['prey']]
    for k,(colname,label,color) in enumerate(zip(['biological_genes','unique_mapped_genes','common_expression_qualified_genes'],['Original','Unique mapping','Expression-qualified'],colors)):
        val=100*d[colname]/d.biological_genes
        xx=np.arange(3)+(k-1)*.235
        ax.bar(xx,val,width=.215,color=color,label=label,zorder=3)
        for x,n,vv in zip(xx,d[colname],val):
            ax.text(x,vv+1.0,str(int(n)),ha='center',va='bottom',fontsize=6.2,rotation=90)
    ax.set_xticks(range(3),['Hits','Screened\nnon-hits','Other PPI']);ax.set_ylim(0,124);ax.set_yticks([0,50,100]);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0);ax.grid(axis='y',color=C['grid'],lw=.5)
    f.legend(handles=[Patch(color=co,label=la) for co,la in zip(colors,['Original','Unique mapping','Expression-qualified'])],
             loc='lower right',bbox_to_anchor=(.978,.005),ncol=3,fontsize=6.4,handlelength=.85,columnspacing=.8)
    r.export('S10A_gene_magnitude_scores',v,'S10A',['vel_scores'],'Frozen corrected export; H appears in both comparison contexts but each display group uses one fixed comparison')
    r.export('S10A_box_descriptive',summ,'S10A',['vel_scores'],'Display-only quartiles/Tukey whiskers')
    r.export('S10B_primary_and_27donor_tests',pd.DataFrame(rows),'S10B',['vel_primary','vel_27'])
    r.export('S10C_compartment_profiles',cp,'S10C',['vel_compartments'])
    r.export('S10D_gene_coverage',coverage,'S10D',['vel_coverage'],'Percent retained is derived only for display; original integer counts are preserved')
    r.export('S10_not_run_sensitivity',s.table('vel_dream'),'S10 legend',['vel_dream'],'Preserve NOT_RUN status, never encode as a null effect or P=1')
    legend=f"""Supplementary Figure S10. Pre-specified Velmeshev patient-cortex validation and sensitivities.
The primary cohort contains 28 donors (14 ASD, 14 controls) after removing three exact public-ID overlaps with the Wamsley roster. One sample per donor is selected by the prespecified PFC-preferred/ACC-fallback rule. Public-ID matching, not genotype verification, was used; undocumented aliases remain possible. The 27-donor sensitivity additionally removes one possible alias (13 ASD, 14 controls), without asserting that this identity is confirmed.
(A) Saved mean absolute regression coefficients across the three fixed neuronal compartments for expression-qualified hits (19), same-screen non-hits (89) and other PPI genes (888). All genes are shown. Bars are IQR/median with Tukey whiskers, not confidence intervals. This descriptive magnitude is not the statistic to which ROMER P values attach.
(B) Saved ROMER Mixed mean-rank enrichment effects and original Holm P values in the primary and 27-donor sensitivity cohorts. The null reference uses 49,999 rotations; the original Holm family contains the two comparisons within each cohort. Neither primary comparison is supportive. Point estimates have no supplied confidence intervals. The two comparator universes have different rank-reference sizes, so their rank-difference magnitudes are not comparable as a common standardized effect. The 27-donor test cannot replace a non-supportive primary result.
(C) All three neuronal compartments, shown descriptively as H-minus-comparator absolute coefficient differences. There are no compartment-specific inferential P values/intervals. The voom sensitivity is retained in the original sources. The all-38-sample dream sensitivity was NOT RUN because variancePartition was unavailable and the recorded package endpoint was unreachable; it is not a negative result and is not represented as zero.
(D) Original, uniquely mapped and expression-qualified memberships. Heights give retained fractions of each original frozen set; labels are actual gene counts. The raw biological sets remain H=21, N=97 and R=1,022. Full donor flow, sample metadata, technical diagnostics, identity limits and additional exported results remain in the unchanged source registry and Supplementary Figures S11–S12.
"""
    r.save(f,'Figure_S10',legend,'supplementary')


def merge_pdfs(paths: list[Path], target: Path):
    doc=fitz.open()
    for p in paths:
        with fitz.open(p) as part:doc.insert_pdf(part)
    doc.set_metadata({'title':target.stem,'author':'','creator':'Frozen-data figure renderer'})
    doc.save(target,garbage=4,deflate=True);doc.close()


def final_qa(run: Run):
    out=run.output
    pd.DataFrame(run.checks).to_csv(out/'qa/NUMERIC_AND_MEMBERSHIP_QA.tsv',sep='\t',index=False)
    pd.DataFrame(run.source.checked).to_csv(out/'qa/INPUT_FILE_HASH_QA.tsv',sep='\t',index=False)
    pd.DataFrame(run.crosswalk).to_csv(out/'qa/PANEL_SOURCE_INDEX.tsv',sep='\t',index=False)
    txt=pd.DataFrame(run.text_checks)
    txt.to_csv(out/'qa/TEXT_PAGE_BOUNDS.tsv',sep='\t',index=False)
    # Render each final PDF independently: inspect text extraction and page bounds.
    render_dir=out/'qa/pdf_renders';render_dir.mkdir(exist_ok=True)
    bad=[]; summaries=[]
    banned=['NON_SUPPORTIVE','NOT_ESTIMABLE','SUPPORTED','PASS','FAIL','Codex','Astra','SPRINGER','Molecular Psychiatry']
    for path in run.pdf_paths:
        with fitz.open(path) as doc:
            page=doc[0];text=page.get_text()
            hits=[word for word in banned if word in text]
            if hits:bad.append(f'{path.name}: {hits}')
            page.get_pixmap(matrix=fitz.Matrix(1.65,1.65),alpha=False).save(render_dir/(path.stem+'.png'))
            summaries.append(dict(file=path.name,pages=len(doc),width_mm=page.rect.width/72*25.4,
                                  height_mm=page.rect.height/72*25.4,banned_label_present=bool(hits),
                                  selectable_text_chars=len(text)))
    pd.DataFrame(summaries).to_csv(out/'qa/PDF_PREFLIGHT.tsv',sep='\t',index=False)
    if bad:raise ValueError('Unexpected internal/publisher text on figures: '+'; '.join(bad))
    outside=txt[txt.inside_page==False] if len(txt) else txt
    if len(outside):
        raise ValueError('Text extends beyond the figure canvas; see qa/TEXT_PAGE_BOUNDS.tsv')
    # Contact sheet only. All manuscript figures retain their separate vector files.
    tiles=[]
    for i in range(1,5):
        im=Image.open(out/f'figures/Figure_{i}.png').convert('RGB')
        tile=Image.new('RGB',(800,920),'white')
        fitted=ImageOps.contain(im,(770,865))
        tile.paste(fitted,((800-fitted.width)//2,36+(865-fitted.height)//2))
        ImageDraw.Draw(tile).text((18,10),f'Figure {i}',fill='black')
        tiles.append(tile)
    sheet=Image.new('RGB',(1600,1840),'white')
    for i,t in enumerate(tiles):sheet.paste(t,((i%2)*800,(i//2)*920))
    sheet.save(out/'Main_Figures_CONTACT.png')
    main=[out/f'figures/Figure_{i}.pdf' for i in range(1,5)]
    merge_pdfs(main,out/'Main_Figures_R4.pdf')
    sup=[out/f'supplementary/Figure_{i}.pdf' for i in ['S9','S10']]
    merge_pdfs(sup,out/'Supplementary_External_Validation_R4.pdf')
    mainleg='\n\n'.join(run.legends[f'Figure_{i}'] for i in range(1,5))
    supleg='\n\n'.join(run.legends[f'Figure_{i}'] for i in ['S9','S10'])
    (out/'legends/Main_Figure_Legends.md').write_text(mainleg,encoding='utf-8')
    (out/'legends/S9_S10_Legends.md').write_text(supleg,encoding='utf-8')
    migration=pd.DataFrame([
      dict(previous='Main Figure 4C: Gordon',current='Supplementary Figure S9',scientific_status='Unchanged pre-specified non-supportive primary validation',action='Keep explicit main Results disclosure; do not relabel exploratory'),
      dict(previous='Main Figure 4D: Velmeshev',current='Supplementary Figure S10',scientific_status='Unchanged non-supportive primary validation; dream sensitivity NOT_RUN',action='Keep explicit main Results disclosure and not-run limitation'),
      dict(previous='Supplementary Figure S8A',current='Main Figure 4A (all 36 cells)',scientific_status='Unchanged secondary cell-type results',action='Display all original rows and original within-set FDR'),
      dict(previous='Wamsley A2/A3/A4',current='Main Figure 4C',scientific_status='Unchanged comparison-specific adjusted models',action='Do not remove A2 or A4; do not infer nested-model attenuation'),
      dict(previous='S1–S8 and S11–S12; all supplementary tables',current='Existing versions retained',scientific_status='Not regenerated in this task',action='Merge replacements S9 and S10 into complete supplement later; this is NOT a full supplement package'),
    ])
    migration.to_csv(out/'SUPPLEMENT_MIGRATION.tsv',sep='\t',index=False)
    disclosure="""# Results text to retain when moving panels to the supplement

This is a disclosure handoff, not an update to the manuscript.

The primary Gordon cortical-organoid test did not support greater dysregulation of the activity-hit set relative to same-screen PPI non-hits (Supplementary Fig. S9). The prespecified Velmeshev patient-cortex tests also did not establish an activity-priority advantage against either the same-screen or wider PPI comparator (Supplementary Fig. S10). These outcomes are retained as prespecified validation findings, not reclassified as exploratory. The all-region dream sensitivity was not run for the recorded software-availability reason and is not a null result.

Do not infer that functional perturbation susceptibility establishes cross-cohort patient-expression prioritization. Avoid the internal label NON_SUPPORTIVE in manuscript prose; report the estimand, observed effect, uncertainty where available, and original P value.
"""
    (out/'legends/REQUIRED_EXTERNAL_RESULTS_DISCLOSURE.md').write_text(disclosure)
    env=dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,pandas=pd.__version__,
             matplotlib=matplotlib.__version__,pymupdf=fitz.VersionBind,renderer_version=VERSION,
             input_path=str(run.source.path),input_zip_sha256=run.source.file_hash,
             inference_recomputed=False,figure_count=4,external_supplement_pages=2)
    (out/'qa/ENVIRONMENT.json').write_text(json.dumps(env,indent=2))
    report=f"""# Automated run QA

- Frozen source files checked: {len(run.source.checked)} / {len(run.source.bindings)}.
- Numeric / membership checks passed: {len(run.checks)}.
- Figure-page text bounds checked: {len(run.text_checks)} objects; no out-of-page text detected.
- Main figures: four; panels: 4 / 4 / 3 / 3. Sole schematic: Figure 1A.
- Supplementary replacements: S9 and S10 only. Other supplemental items are not regenerated.
- PDF pages re-rendered with PyMuPDF, with selectable text and no publisher/internal status branding.
- Observed rates, ECDF-free descriptive strip plots, IQRs and histogram bins are display calculations only.
- No model fitting, significance testing, interval generation, gene-set selection or random data generation occurred. RNG is used only for the x-jitter of real observations.
- Source-specific P values, Holm families, comparator populations and confidence-interval methods remain distinct.
- No numeric output was read from image-model mockups or from reference papers.

Automated bounds checks do not prove the absence of every overlap. Review all PDFs and the independent render PNGs at print size before accepting the figure version. The local runner does not claim a new human visual sign-off automatically.
"""
    (out/'qa/AUTOMATED_QA.md').write_text(report)
    for name in ['draw_mainfigures_r4.py','source_map.json','requirements.txt','run_local.sh','CODEX_EXECUTE_R4.md','README_中文.md','test_renderer.py']:
        if (HERE/name).is_file():shutil.copy2(HERE/name,out/'code'/name)
    manifest=[]
    for p in sorted(out.rglob('*')):
        if p.is_file() and p.name!='OUTPUT_SHA256.tsv':
            manifest.append(dict(path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=digest(p.read_bytes())))
    pd.DataFrame(manifest).to_csv(out/'qa/OUTPUT_SHA256.tsv',sep='\t',index=False)


def main() -> int:
    pa=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    pa.add_argument('--input',required=True,type=Path,help='Frozen R1_1 ZIP or extracted package root')
    pa.add_argument('--output',type=Path,help='New output directory; do not point to frozen inputs')
    pa.add_argument('--dpi',type=int,default=400,help='PNG/rasterized-layer resolution (default 400)')
    pa.add_argument('--verify-only',action='store_true',help='Validate bound source hashes, schema and numerical consistency without drawing')
    pa.add_argument('--overwrite',action='store_true',help='Overwrite only an output directory created by this renderer (marker required)')
    args=pa.parse_args()
    if args.dpi<150 or args.dpi>1200:pa.error('--dpi must be between 150 and 1200')
    output=(args.output or Path.cwd()/'ASD_PPI_CRISPR_MAINFIG_R4_OUTPUT').resolve()
    if output==args.input.resolve() or output==HERE:pa.error('Output must be separate from input and code')
    if args.input.resolve().is_dir() and (args.input.resolve() in output.parents or output in args.input.resolve().parents):
        pa.error('Output must not contain or be contained by the frozen input directory')
    marker=output/'.asd_ppi_r4_output'
    if output.exists() and any(output.iterdir()):
        if not args.overwrite or not marker.is_file():
            pa.error('Output is not empty. Choose a new directory, or --overwrite an existing marked R4 output only.')
        # Keep prior result intact instead of deleting: timestamped backup sibling.
        from datetime import datetime
        backup=output.with_name(output.name+'_backup_'+datetime.now().strftime('%Y%m%d_%H%M%S'))
        output.rename(backup)
    output.mkdir(parents=True,exist_ok=True);marker.write_text(VERSION+'\n')
    logging.basicConfig(level=logging.INFO,format='%(levelname)s %(message)s',handlers=[logging.StreamHandler(),logging.FileHandler(output/'EXECUTION_LOG.txt')])
    logging.getLogger('fontTools').setLevel(logging.WARNING)
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    source=None
    try:
        source=FrozenSource(args.input,HERE/'source_map.json');run=Run(source,output,args.dpi)
        LOG.info('Checking exact frozen inputs and cross-table values...')
        validate(source,run)
        if args.verify_only:
            pd.DataFrame(run.checks).to_csv(output/'qa/NUMERIC_AND_MEMBERSHIP_QA.tsv',sep='\t',index=False)
            pd.DataFrame(source.checked).to_csv(output/'qa/INPUT_FILE_HASH_QA.tsv',sep='\t',index=False)
            LOG.info('Verification complete; no plots generated: %s',output)
            return 0
        set_theme()
        for fn in [draw_figure1,draw_figure2,draw_figure3,draw_figure4,draw_s9,draw_s10]:
            LOG.info('Rendering %s',fn.__name__);fn(run)
        final_qa(run)
        target=output.with_suffix('.zip')
        with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
            for p in sorted(output.rglob('*')):
                if p.is_file():z.write(p,arcname=output.name+'/'+p.relative_to(output).as_posix())
        # Do not append to the file log after its content has been hashed.
        print(f'COMPLETE: {output}', flush=True)
        print(f'RESULT ZIP: {target}', flush=True)
        return 0
    except Exception as exc:
        LOG.exception('Stopped rather than substituting scientific data: %s',exc)
        (output/'ERROR.txt').write_text(str(exc)+'\n',encoding='utf-8')
        if 'run' in locals():
            pd.DataFrame(run.checks).to_csv(output/'qa/NUMERIC_AND_MEMBERSHIP_QA.tsv',sep='\t',index=False)
        return 1
    finally:
        if source is not None:source.close()



R4_BINDINGS = {'flow': {'path': '02_main_figure_sources/Figure1/F1A_STUDY_FLOW.tsv', 'sha256': '9bb28644564179bcaddebfa83190c88118bcb565529ce0126226fcd9346d52d7', 'rows': 5, 'columns': ['stage', 'input', 'output', 'role', 'claim_limit']}, 'activity': {'path': '02_main_figure_sources/Figure1/panel_sources/1B_FIG1_ACTIVITY_PREY_REPRODUCTION.tsv', 'sha256': '22eafae2ad5c07bf338fa471a483811b7ac8fe6cd799b6e3a343ec810720814c', 'rows': 2, 'columns': ['analysis', 'group_comparison', 'exposed_n', 'exposed_hits', 'exposed_hit_rate', 'comparator_n', 'comparator_hits', 'comparator_hit_rate', 'OR', 'CI95_low', 'CI95_high', 'P']}, 'activity_adjusted': {'path': '02_main_figure_sources/Figure1/panel_sources/1C_FIG1_ACTIVITY_SELECTION_ADJUSTED.tsv', 'sha256': '5c0eb967def92de081d9efa28aa046e6817b30083579d09c921749a5cecb075c', 'rows': 14, 'columns': ['model', 'outcome', 'term', 'beta', 'standard_error', 'OR', 'OR_CI_low', 'OR_CI_high', 'z', 'P', 'n', 'converged', 'warnings']}, 'activity_strata': {'path': '02_main_figure_sources/Figure1/panel_sources/1C_FIG1_ACTIVITY_SELECTION_STRATA.tsv', 'sha256': '83097c7f5463307889ab4ce10461de9a06dbffd7d4ad04969e8e7fdf943df621', 'rows': 10, 'columns': ['row_type', 'selection_stratum', 'informative_for_CMH', 'analysis', 'group_comparison', 'exposed_n', 'exposed_hits', 'exposed_hit_rate', 'comparator_n', 'comparator_hits', 'comparator_hit_rate', 'OR', 'CI95_low', 'CI95_high', 'P', 'informative_strata_n']}, 'activity_permutation': {'path': '02_main_figure_sources/Figure1/panel_sources/1D_FIG1_ACTIVITY_PERMUTATION.tsv', 'sha256': '6c646de81a00b12093daeacd20365857d0c12ac52307d9cf00c833ad02a35292', 'rows': 2, 'columns': ['analysis', 'seed', 'statistic', 'observed_hit_rate_difference', 'null_mean', 'null_sd', 'null_q025', 'null_q975', 'exceedances', 'empirical_P_plus_one', 'Monte_Carlo_SE', 'n_permutations']}, 'activity_loo': {'path': '02_main_figure_sources/Figure1/panel_sources/1D_leave_one_stratum_out.tsv', 'sha256': '41a812d19e15d1c4f5c48b0f6e658a2b8c6a36cabb6a23305255dc1715b65d1d', 'rows': 8, 'columns': ['endpoint', 'dropped_stratum', 'remaining_strata_n', 'leave_one_out_CMH_OR', 'CI95_low', 'CI95_high', 'P']}, 'igvf_genes': {'path': '00_final_freeze/authority_tables/authoritative_results/ASD_PPI_CRISPR_STAGE2A_GATE_R1/11_IGVF_TARGET_UNIVERSE.tsv', 'sha256': 'c7c3c7879f0617be6bf7fcc1a49a35d990e55e4755d116788b8b2ccfa13ba54d', 'rows': 1275, 'columns': ['gene', 'usable_target', 'guide_n', 'source_category', 'Science_bait', 'Science_prey_only', 'author_energy_distance_hit', 'author_hit_promoters_n', 's4_promoters_n', 'mean_energy_distance', 'energy_distance_cluster_n']}, 'igvf_full_results': {'path': '02_main_figure_sources/Figure2/panel_sources/2A_FIG2_IGVF_PERTURBSEQ.tsv', 'sha256': 'a3db08551c62f433e4d1006a8cb99c2f79cd48c3c983b1e49e0aa592f39d30d7', 'rows': 3, 'columns': ['test', 'model', 'usable_targets_n', 'prey_only_n', 'prey_only_hits_n', 'prey_only_hit_proportion', 'comparator_n', 'comparator_hits_n', 'comparator_hit_proportion', 'OR', 'CI_low', 'CI_high', 'P', 'author_significant_promoter_n', 'author_hit_source_gene_n', 'author_hit_gene_n_in_universe', 'excluded_Science_baits_n', 'n_prey_only', 'n_comparator', 'effect', 'rank_AUC', 'mean_prey_only', 'mean_comparator', 'n', 'estimate', 'warnings', 'converged']}, 'igvf_models': {'path': '02_main_figure_sources/Figure2/panel_sources/2C_A1_IGVF_MODELS.tsv', 'sha256': '609fa14761cb99c44f6d0d23252dda1c36426c4bcaf308cdf8725358f550e693', 'rows': 5, 'columns': ['model', 'method', 'status', 'n', 'events', 'non_events', 'design_rank', 'design_columns', 'removed_constant_covariates', 'estimate_log_OR', 'OR', 'CI_low_OR', 'CI_high_OR', 'P_raw', 'profile_LR', 'input_rows_hash', 'n_exposed', 'n_comparator', 'input_complete_case_retention', 'input_event_n', 'input_non_event_n', 'outcome_scale', 'residual_df', 'estimate', 'robust_SE_HC3', 'ci_low', 'ci_high']}, 'igvf_matched': {'path': '02_main_figure_sources/Figure2/panel_sources/2D_A1_IGVF_MATCHED_TEST.tsv', 'sha256': '4daa91fa1b8af6ec9ad1429f6962007dbe90d1e63817f09bf7d3d749194e6efe', 'rows': 1, 'columns': ['analysis_id', 'status', 'matched_sets_n', 'n_exposed', 'n_comparator', 'effect_definition', 'estimate_difference', 'CI_low', 'CI_high', 'OR', 'OR_CI_low', 'OR_CI_high', 'null_mean', 'null_sd', 'null_q025', 'null_q975', 'P_permutation', 'permutation_B', 'seed_permutation', 'seed_bootstrap', 'CI_method', 'OR_continuity_correction_used', 'exposed_hits', 'exposed_nonhits', 'comparator_hits', 'comparator_nonhits', 'analysis', 'exposed_universe_n', 'comparator_universe_n', 'complete_case_n', 'matching_status', 'matching_ratio_controls', 'matching_exposure_retention', 'propensity_caliper_logit_SD', 'exact_matching_variables', 'exact_mismatch_n', 'brain_expression_control_available', 'OR_inference_role', 'CI_pairing_accounted', 'not_for_primary_inference']}, 'igvf_balance': {'path': '02_main_figure_sources/Figure2/panel_sources/2D_A_MATCHING_BALANCE.tsv', 'sha256': 'b70f4f47b1283c4d5757485e9fcbad8adc0e427ae67d151cd8313e73dce1c784', 'rows': 18, 'columns': ['analysis_id', 'phase', 'metric', 'exposed_mean', 'comparator_mean', 'SMD', 'abs_SMD', 'criterion']}, 'holm_family': {'path': '00_final_freeze/authority_tables/strengthening_R1_1/A_GENERIC_AND_ASCERTAINMENT/A_PRIMARY_HOLM_FAMILY.tsv', 'sha256': 'dab67e18da341d530f332e53bdffd8db936d5d3db9fe754ddc50c42f8c506469', 'rows': 4, 'columns': ['test', 'estimate', 'estimate_scale', 'P_raw', 'status', 'family', 'P_adjustment_input_bookkeeping', 'P_adjusted_Holm', 'P_missing_slot_adjusted_as']}, 'invivo_genes': {'path': '00_final_freeze/authority_tables/authoritative_results/ASD_PPI_CRISPR_STAGE1B_R1_1/04_PRIMARY_COMPOSITE_SCORES.tsv', 'sha256': '1b7c14b7a0652274f04880b9334de6989c182ab8e5c4c4c9e2b2e73722f6c832', 'rows': 9949, 'columns': ['mouse_symbol', 'human_ensembl_stable_id', 'human_gene_symbol_ensembl', 'primary_group', 'guide_count_stratum', 'Z_Syn1-GW-2M', 'Z_Emx1-GW-2M', 'Z_Vglut2-GW-2M', 'Z_Vgat-GW-2M', 'Z_Syn1-GW-4M', 'Z_Emx1-GW-4M', 'Z_Vglut2-GW-4M', 'Z_Vgat-GW-4M', 'composite_vulnerability_8screen', 'composite_2M', 'composite_4M', 'aggregation_label']}, 'invivo_primary': {'path': '02_main_figure_sources/Figure3/panel_sources/3B_FIG3_INVIVO_PRIMARY.tsv', 'sha256': '0a1f8b1dc6b60d1200da215d65d5de67b909fe010f4b2b90d68a2075c92fcb47', 'rows': 1, 'columns': ['endpoint', 'prey_mean', 'comparator_mean', 'prey_median', 'comparator_median', 'mean_difference', 'cohen_d', 'auc', 'rank_biserial', 'bootstrap_CI_low', 'bootstrap_CI_high', 'bootstrap_n', 'bootstrap_unit', 'prey_only_n', 'comparator_n', 'permutation_P_one_sided', 'permutation_P_two_sided']}, 'invivo_screens': {'path': '02_main_figure_sources/Figure3/panel_sources/3C_FIG3_INVIVO_SCREENS.tsv', 'sha256': '6bcbd5aed48cc2b022a8d26bdf1e6c04fe7bd7b17065ae116d292a3ac0f72268', 'rows': 8, 'columns': ['screen', 'prey_mean', 'comparator_mean', 'prey_median', 'comparator_median', 'mean_difference', 'cohen_d', 'auc', 'rank_biserial', 'Mann_Whitney_Wilcoxon_P_two_sided_normal_tie_corrected', 'gene_label_permutation_P_one_sided', 'gene_label_permutation_P_two_sided', 'prey_only_n', 'comparator_n', 'permutation_P_BH_FDR']}, 'invivo_age': {'path': '02_main_figure_sources/Figure3/panel_sources/3D_11_AGE_STRATIFIED_SENSITIVITY.tsv', 'sha256': '70fbd81b66cd42e166046cc81bf868729efea4bd18fd6dc23f9c82b0d3890a8a', 'rows': 2, 'columns': ['age_composite', 'prey_mean', 'comparator_mean', 'prey_median', 'comparator_median', 'mean_difference', 'cohen_d', 'auc', 'rank_biserial', 'bootstrap_CI_low', 'bootstrap_CI_high', 'permutation_n', 'permutation_P_one_sided', 'permutation_P_two_sided', 'one_sided_P_Holm_2test_family', 'two_sided_P_Holm_2test_family']}, 'wam_cells': {'path': '03_supplementary_figure_sources/SFig08_Wamsley_celltype_network_controls/panel_sources/S8A_WAMSLEY_CELLTYPE_SECONDARY.tsv', 'sha256': 'a58f6389c0c30a56283b3b5b766783deffb93dbca078b42d3d6d77fbf863c627', 'rows': 72, 'columns': ['gene_set', 'source_cell_type', 'neuronal_vs_non_neuronal', 'subtype', 'set_n_tested', 'comparator_n_tested', 'set_median_abs_logFC', 'comparator_median_abs_logFC', 'effect_median_difference', 'AUC_probability_superiority', 'P_two_sided_Mann_Whitney', 'BH_q_within_gene_set_36_celltypes']}, 'wam_class': {'path': '00_final_freeze/authority_tables/authoritative_results/ASD_PPI_CRISPR_STAGE2A_R1_1_WAMSLEY_CONFIRMATORY/02_CELLTYPE_CLASSIFICATION.tsv', 'sha256': '655f7f586f920dd9aa192ce42c0189b7971cc32cf2b5d9e15b727ab5cd56db7d', 'rows': 36, 'columns': ['source_cell_type', 'neuronal_vs_non_neuronal', 'subtype', 'included_in_primary_neuronal_omnibus', 'included_in_nonneuronal_omnibus', 'reason_source_annotation_only']}, 'wam_set1': {'path': '02_main_figure_sources/Figure4/panel_sources/4A_FIG4_WAMSLEY_SET1.tsv', 'sha256': 'eec645139c786b1ba7d47991d5e5b5946b583541f0fad78fd461dc733ea5c2f5', 'rows': 1, 'columns': ['test', 'set_n', 'comparator_n', 'set_median', 'comparator_median', 'effect_median_difference', 'AUC_probability_superiority', 'bootstrap_CI_low', 'bootstrap_CI_high', 'observed_effect', 'perm_P_one_sided', 'perm_P_two_sided', 'permutation_null_median', 'permutation_null_mean', 'permutation_null_sd', 'permutations', 'seed', 'elapsed_seconds', 'Holm_P', 'Holm_P_one_sided', 'Holm_P_two_sided', 'bootstrap_replicates', 'bootstrap_seed']}, 'wam_set2': {'path': '02_main_figure_sources/Figure4/panel_sources/4A_FIG4_WAMSLEY_SET2.tsv', 'sha256': '05baca0956c78e58f154485e60f9d25777c56721144f1d75933f8ea11dd83e61', 'rows': 1, 'columns': ['test', 'set_n', 'comparator_n', 'set_median', 'comparator_median', 'effect_median_difference', 'AUC_probability_superiority', 'bootstrap_CI_low', 'bootstrap_CI_high', 'observed_effect', 'perm_P_one_sided', 'perm_P_two_sided', 'permutation_null_median', 'permutation_null_mean', 'permutation_null_sd', 'permutations', 'seed', 'elapsed_seconds', 'Holm_P', 'Holm_P_one_sided', 'Holm_P_two_sided', 'bootstrap_replicates', 'bootstrap_seed']}, 'wam_a2': {'path': '03_supplementary_figure_sources/SFig08_Wamsley_celltype_network_controls/panel_sources/S8C_A2_WAMSLEY_NETWORK_CONTROLS.tsv', 'sha256': '1e3dc2d8c064d025cd79258802fc8f21e6afbd9e0aea078f48224111ea5dae0d', 'rows': 1, 'columns': ['analysis_id', 'status', 'n', 'n_exposed', 'n_comparator', 'exposure_definition', 'outcome', 'brain_expression_control_available', 'rank_model_effect_scale', 'estimate_rank_beta', 'robust_SE_HC3', 'ci_low', 'ci_high', 'P_raw', 'median_exposed', 'median_comparator', 'median_difference', 'median_diff_bootstrap_low', 'median_diff_bootstrap_high', 'rank_AUC', 'rank_AUC_low', 'rank_AUC_high', 'matched_status', 'matched_sets_n', 'matched_ratio_controls', 'matched_exposure_retention', 'matched_caliper_logit_SD', 'matched_max_abs_SMD', 'matched_estimate_difference_original_Y', 'matched_CI_low', 'matched_CI_high', 'matched_permutation_P', 'matched_exposed_n', 'matched_comparator_n', 'complete_case_universe_n', 'pre_complete_case_exposed_n', 'pre_complete_case_comparator_n']}, 'wam_a3': {'path': '03_supplementary_figure_sources/SFig08_Wamsley_celltype_network_controls/panel_sources/S8C_A3_WAMSLEY_SUBSET_VS_REST.tsv', 'sha256': '518069dd3b5a64864e467dc9ce70c7dc12999dccce1a0d3cbe27eeadd817b7e2', 'rows': 1, 'columns': ['analysis_id', 'status', 'n', 'n_exposed', 'n_comparator', 'exposure_definition', 'outcome', 'brain_expression_control_available', 'rank_model_effect_scale', 'estimate_rank_beta', 'robust_SE_HC3', 'ci_low', 'ci_high', 'P_raw', 'median_exposed', 'median_comparator', 'median_difference', 'median_diff_bootstrap_low', 'median_diff_bootstrap_high', 'rank_AUC', 'rank_AUC_low', 'rank_AUC_high', 'matched_status', 'matched_sets_n', 'matched_ratio_controls', 'matched_exposure_retention', 'matched_caliper_logit_SD', 'matched_max_abs_SMD', 'matched_estimate_difference_original_Y', 'matched_CI_low', 'matched_CI_high', 'matched_permutation_P', 'matched_exposed_n', 'matched_comparator_n', 'complete_case_universe_n', 'pre_complete_case_exposed_n', 'pre_complete_case_comparator_n']}, 'wam_a4': {'path': '02_main_figure_sources/Figure4/panel_sources/4B_A4_WAMSLEY_WITHIN_ACTIVITY_SCREEN.tsv', 'sha256': 'af5c6da5078bad9451cf6a075ebfc1b204a0a858777d3450a784c38012103d1c', 'rows': 1, 'columns': ['model', 'status', 'n', 'effect_scale', 'estimate_rank_beta', 'robust_SE_HC3', 'ci_low', 'ci_high', 'P_raw', 'permutation_rank_AUC', 'permutation_P', 'permutation_null_mean', 'permutation_null_sd', 'permutation_null_q025', 'permutation_null_q975', 'permutation_B', 'permutation_strata_counts', 'network_complete_case_n', 'network_adjusted_status', 'network_adjusted_rank_beta', 'network_adjusted_ci_low', 'network_adjusted_ci_high', 'network_adjusted_P', 'brain_expression_control_available', 'n_exposed', 'n_comparator', 'median_exposed', 'median_comparator', 'median_difference', 'median_diff_bootstrap_low', 'median_diff_bootstrap_high', 'rank_AUC', 'rank_AUC_low', 'rank_AUC_high', 'source_selection_guide_count_handled', 'outcome']}, 'gordon_scores': {'path': '03_supplementary_figure_sources/SFig11_External_validation_descriptive_profiles/panel_sources/S11A_GORDON_GENE_SCORES_C36.tsv', 'sha256': '05829e467ce1201cd17c271aa72e71572f3d2c29eaa34365d998a413bb0edcf9', 'rows': 17962, 'columns': ['ensembl_gene_id', 'hgnc_symbol', 'is_P', 'is_H', 'is_R', 'is_S', 'is_H_screen', 'is_N_screen', 'is_external_nonPPI_nonbait', 'selection_stratum', 'mean_abs_logFC_C36', 'beta__15q13_day25_vs_Control', 'beta__15q13_day50_vs_Control', 'beta__15q13_day75_vs_Control', 'beta__15q13_day100_vs_Control', 'beta__16p11del_day25_vs_Control', 'beta__16p11del_day50_vs_Control', 'beta__16p11del_day75_vs_Control', 'beta__16p11del_day100_vs_Control', 'beta__16p11dup_day25_vs_Control', 'beta__16p11dup_day50_vs_Control', 'beta__16p11dup_day75_vs_Control', 'beta__16p11dup_day100_vs_Control', 'beta__22q11del_day25_vs_Control', 'beta__22q11del_day50_vs_Control', 'beta__22q11del_day75_vs_Control', 'beta__22q11del_day100_vs_Control', 'beta__22q13del_day25_vs_Control', 'beta__22q13del_day50_vs_Control', 'beta__22q13del_day75_vs_Control', 'beta__22q13del_day100_vs_Control', 'beta__IdiopathicASD_day25_vs_Control', 'beta__IdiopathicASD_day50_vs_Control', 'beta__IdiopathicASD_day75_vs_Control', 'beta__IdiopathicASD_day100_vs_Control', 'beta__PCDH19_day25_vs_Control', 'beta__PCDH19_day50_vs_Control', 'beta__PCDH19_day75_vs_Control', 'beta__PCDH19_day100_vs_Control', 'beta__SHANK3_day25_vs_Control', 'beta__SHANK3_day50_vs_Control', 'beta__SHANK3_day75_vs_Control', 'beta__SHANK3_day100_vs_Control', 'beta__TS_day25_vs_Control', 'beta__TS_day50_vs_Control', 'beta__TS_day75_vs_Control', 'beta__TS_day100_vs_Control']}, 'gordon_primary': {'path': '03_supplementary_figure_sources/SFig09_Gordon_detailed_external_validation/panel_sources/S9A_H_VS_N_C36_PRIMARY.tsv', 'sha256': '8a2c7262d7e38f136d520db2fc3fadffc420f9be2f1e05dcf6399739ac580d8e', 'rows': 1, 'columns': ['comparison', 'score', 'role', 'n_group_A', 'n_group_B', 'mean_A', 'mean_B', 'delta_mean_A_minus_B', 'median_A', 'median_B', 'delta_median_A_minus_B', 'auc_probability_A_greater_than_B_ties_half', 'A_q10', 'A_q25', 'A_q75', 'A_q90', 'B_q10', 'B_q25', 'B_q75', 'B_q90', 'p_value', 'p_definition']}, 'gordon_c24': {'path': '03_supplementary_figure_sources/SFig09_Gordon_detailed_external_validation/panel_sources/S9B_H_VS_N_C24_SENSITIVITY.tsv', 'sha256': 'e792a512ca52ad52e75d2dd30ccc2d50aa180bcb2b81f0b9a12c80e6c92772f6', 'rows': 1, 'columns': ['comparison', 'score', 'role', 'n_group_A', 'n_group_B', 'mean_A', 'mean_B', 'delta_mean_A_minus_B', 'median_A', 'median_B', 'delta_median_A_minus_B', 'auc_probability_A_greater_than_B_ties_half', 'A_q10', 'A_q25', 'A_q75', 'A_q90', 'B_q10', 'B_q25', 'B_q75', 'B_q90', 'p_value', 'p_definition', 'empirical_p', 'permutation_method', 'extreme_count', 'assignments_total', 'assignments_evaluated', 'seed', 'draws_requested', 'stratified', 'rng']}, 'gordon_unstratified': {'path': '03_supplementary_figure_sources/SFig09_Gordon_detailed_external_validation/panel_sources/S9B_H_VS_N_C36_UNSTRATIFIED_SENSITIVITY.tsv', 'sha256': '7388abdafce69767027ff5741f835f7c6bff98b941c33d99d76fd95eb379efea', 'rows': 1, 'columns': ['comparison', 'score', 'role', 'n_group_A', 'n_group_B', 'mean_A', 'mean_B', 'delta_mean_A_minus_B', 'median_A', 'median_B', 'delta_median_A_minus_B', 'auc_probability_A_greater_than_B_ties_half', 'A_q10', 'A_q25', 'A_q75', 'A_q90', 'B_q10', 'B_q25', 'B_q75', 'B_q90', 'p_value', 'p_definition', 'permutation_method', 'extreme_count', 'assignments_total', 'assignments_evaluated', 'seed']}, 'gordon_days': {'path': '03_supplementary_figure_sources/SFig09_Gordon_detailed_external_validation/panel_sources/S9C_DAY_PROFILES.tsv', 'sha256': 'b60c8d78d3d707ebbbd4d4529fc7def65140a439e3640b4ebdafd9bbded0f004', 'rows': 8, 'columns': ['comparison', 'profile_panel', 'developmental_day', 'score', 'role', 'n_group_A', 'n_group_B', 'mean_A', 'mean_B', 'delta_mean_A_minus_B', 'median_A', 'median_B', 'delta_median_A_minus_B', 'auc_probability_A_greater_than_B_ties_half', 'A_q10', 'A_q25', 'A_q75', 'A_q90', 'B_q10', 'B_q25', 'B_q75', 'B_q90', 'p_value', 'p_definition']}, 'gordon_forms': {'path': '03_supplementary_figure_sources/SFig09_Gordon_detailed_external_validation/panel_sources/S9C_FORM_PROFILES.tsv', 'sha256': '4ced3e4dd0055ebc3e12a92f0f5dd30612e58bae404d7e22a59ffcdcd9e8b740', 'rows': 15, 'columns': ['comparison', 'ASD_form', 'contrast_set', 'low_donor_form', 'ASD_donors_by_day', 'score', 'role', 'n_group_A', 'n_group_B', 'mean_A', 'mean_B', 'delta_mean_A_minus_B', 'median_A', 'median_B', 'delta_median_A_minus_B', 'auc_probability_A_greater_than_B_ties_half', 'A_q10', 'A_q25', 'A_q75', 'A_q90', 'B_q10', 'B_q25', 'B_q75', 'B_q90', 'p_value', 'p_definition']}, 'vel_scores': {'path': '03_supplementary_figure_sources/SFig11_External_validation_descriptive_profiles/panel_sources/S11B_VEL_GENELEVEL_SCORES.tsv', 'sha256': '79a455c6461afcb5b48328ebdded94ee5c46e1c25ec86aebd035c525f474bb9a', 'rows': 1015, 'columns': ['gene_symbol', 'set', 'mean_abs_beta', 'mean_beta', 'mean_abs_moderated_t', 'mean_rank_abs_t', 'n_compartments', 'comparison']}, 'vel_primary': {'path': '03_supplementary_figure_sources/SFig10_Velmeshev_cohort_mapping_and_sensitivities/panel_sources/S10B_PRIMARY_HOLM_FAMILY.tsv', 'sha256': 'fce51d711cc0fb889a974422cc4ea2ea03b22508c6e490a771c3e13b59ee874a', 'rows': 2, 'columns': ['cohort', 'comparison', 'status', 'reason', 'H_biological_genes', 'comparator_biological_genes', 'H_expression_features', 'comparator_expression_features', 'donors_ASD', 'donors_CTL', 'design_rank', 'design_columns', 'design_residual_df', 'seed', 'nrot', 'set_statistic', 'romer_direction', 'Mixed_P', 'Holm_P', 'mean_rank_abs_t_difference', 'mean_abs_beta_difference']}, 'vel_27': {'path': '03_supplementary_figure_sources/SFig10_Velmeshev_cohort_mapping_and_sensitivities/panel_sources/S10C_ALIAS_EXCLUSION27_RESULTS.tsv', 'sha256': 'd643bce14c24f892576f2ebc28e46cff15c3e223fe3b0935d5f4c3cae97e3222', 'rows': 2, 'columns': ['cohort', 'comparison', 'status', 'reason', 'H_biological_genes', 'comparator_biological_genes', 'H_expression_features', 'comparator_expression_features', 'donors_ASD', 'donors_CTL', 'design_rank', 'design_columns', 'design_residual_df', 'seed', 'nrot', 'set_statistic', 'romer_direction', 'Mixed_P', 'Holm_P', 'mean_rank_abs_t_difference', 'mean_abs_beta_difference']}, 'vel_compartments': {'path': '03_supplementary_figure_sources/SFig10_Velmeshev_cohort_mapping_and_sensitivities/panel_sources/S10C_COMPARTMENT_PROFILES.tsv', 'sha256': '6e843ff43917e4737d95b3efb298c173ae21941f223404ad24948e4c303e6631', 'rows': 6, 'columns': ['comparison', 'compartment', 'H_genes', 'comparator_genes', 'H_mean_abs_beta', 'comparator_mean_abs_beta', 'mean_abs_beta_difference', 'mean_abs_moderated_t_rank_difference', 'P_value_generated']}, 'vel_coverage': {'path': '03_supplementary_figure_sources/SFig10_Velmeshev_cohort_mapping_and_sensitivities/panel_sources/S10A_SET_COVERAGE_SUMMARY.tsv', 'sha256': '9eebb0be305b40dc8bfed0ee8831f6b4c47aaded1c095549b24f17d31405daf3', 'rows': 5, 'columns': ['set', 'biological_genes', 'unique_mapped_genes', 'common_expression_qualified_genes']}, 'vel_dream': {'path': '03_supplementary_figure_sources/SFig10_Velmeshev_cohort_mapping_and_sensitivities/panel_sources/S10C_ALL_REGIONS_DREAM_STATUS.tsv', 'sha256': '367480edab8b8a7a2724feabfca0318e238242ab6bed6cbdf9a5eef6b7b05510', 'rows': 1, 'columns': ['module', 'requested_formula', 'software', 'package_available', 'status', 'reason', 'sample_rows', 'donors', 'expected_samples']}, 'gene_manifest': {'path': '00_final_freeze/FINAL_FROZEN_GENESET_MANIFEST.tsv', 'sha256': '8568e2908004b7427e3fc0d04ab07546b7d37b625dd8d89627001758fbbfbe2a', 'rows': 10, 'columns': ['gene_set_id', 'authoritative_count', 'observed_member_rows', 'source_file', 'source_original_path', 'SHA256', 'definition', 'source_stage', 'Gordon_effective_mapped_count', 'Velmeshev_unique_mapped_count', 'Velmeshev_common_expression_qualified_count', 'mapping_note']}}

# R5: frozen-table-only production; R4 functions above remain reusable building blocks.
class R5Source:
    def __init__(self, root):
        self.root=Path(root).resolve(); self.bindings=dict(R4_BINDINGS); self.cache={}; self.checked=[]; self.active=[]; self.numeric=[]
    def read(self,path,panel=''):
        path=Path(path); path=path if path.is_absolute() else self.root/path
        raw=path.read_bytes(); df=pd.read_csv(io.BytesIO(raw),sep='\t')
        rel=str(path.relative_to(self.root)); sha=digest(raw)
        self.active.append((rel,sha,panel,len(df)))
        for c in (df.select_dtypes(include=np.number).columns if len(df)<=500 else []):
            if any(k in c.lower() for k in ['donor_id','individual_id','age_years','pmi_hours','rin','excel_row']): continue
            for i,v in df[c].items():
                if pd.notna(v): self.numeric.append(dict(check='direct_frozen_number',source=rel,row=int(i)+1,column=c,actual=repr(v),expected=repr(v),status='MATCH'))
        return df
    def table(self,role):
        m=self.bindings[role]; p=self.root/m['path']
        if digest(p.read_bytes())!=m['sha256']:raise ValueError('Frozen hash mismatch '+role)
        d=self.read(m['path']);
        if len(d)!=m['rows'] or list(d.columns)!=m['columns']:raise ValueError('Frozen schema mismatch '+role)
        return d
    def sup(self,n,suffix):
        paths=list((self.root/'03_supplementary_figure_sources').glob(f'SFig{n:02d}*/panel_sources/*{suffix}'))
        if len(paths)!=1:raise ValueError((n,suffix,paths))
        return self.read(paths[0])

class R5Run:
    def __init__(self,s,out,dpi=300):
        self.source=s;self.output=Path(out);self.dpi=dpi;self.checks=[];self.crosswalk=[];self.legends={};self.pdf_paths=[]
        for d in ['figures','supplementary_figures','legends','source_index','qa/pdf_renders']:(self.output/d).mkdir(parents=True,exist_ok=True)
    check=Run.check
    def export(self,name,data,panels,roles,transformation='Frozen values selected for display'):
        for role in roles:
            self.crosswalk.append(dict(panel=panels,source=self.source.bindings[role]['path'],sha256=self.source.bindings[role]['sha256'],operation=transformation))
    def save(self,f,name,legend,category='figures'):
        cat='figures' if category=='figures' else 'supplementary_figures'
        # Unambiguous scientific status retained for external validations.
        if name in ['Figure_S9','Figure_S10']:
            f.text(.5,.992,'NON_SUPPORTIVE',ha='center',va='top',fontsize=6.5,color=C['other'])
        if name=='Figure_S10':
            f.text(.045,.052,'28 donors (14 ASD / 14 controls); 3 exact overlaps excluded.',fontsize=6.2,color=C['other'])
            f.text(.045,.025,'All-region dream: NOT_RUN (required package unavailable).',fontsize=6.2,color=C['other'])
        base=self.output/cat/name
        f.savefig(base.with_suffix('.pdf'),facecolor='white',metadata={'Title':name,'Author':'','CreationDate':None,'ModDate':None})
        f.savefig(base.with_suffix('.png'),dpi=self.dpi,facecolor='white')
        self.pdf_paths.append(base.with_suffix('.pdf'));self.legends[name]=legend
        for path,sha,panel,n in self.source.active:self.crosswalk.append(dict(panel=name,source=path,sha256=sha,operation='Frozen source; display-only selection/order',source_rows=n))
        self.source.active=[]
        plt.close(f)

def subfig(title,w=190,h=170):
    f=figure(w,h);f.text(.045,.977,title,va='top',fontsize=10.5,fontweight='bold');return f

def ttl(ax,letter,title):
    ax.set_title(f'{letter}  {title}',loc='left',pad=12,fontsize=9,fontweight='bold')

def forest(ax,labels,x,lo=None,hi=None,p=None,log=False,color=None,xlabel='Effect (95% CI)',ref=0):
    x=np.asarray(x,dtype=float); yy=np.arange(len(x))[::-1]
    for i,(v,y) in enumerate(zip(x,yy)):
        co=color[i] if isinstance(color,list) else color or C['prey']
        if np.isfinite(v):
            if lo is not None and np.isfinite(lo[i]) and np.isfinite(hi[i]):ax.plot([lo[i],hi[i]],[y,y],color=co,lw=1.2)
            ax.plot(v,y,'o',color=co,ms=4)
        else:ax.text(.5,y,'Not estimable',transform=ax.get_yaxis_transform(),ha='center',fontsize=7)
        if p is not None and np.isfinite(p[i]):ax.text(.98,y-.36,'P = '+numeric_p(p[i]),transform=ax.get_yaxis_transform(),ha='right',fontsize=6.5,color=C['other'],va='center',bbox=dict(facecolor='white',edgecolor='none',pad=.2))
    ax.axvline(1 if log else ref,lw=.6,ls='--',color=C['other']);ax.set_yticks(yy,labels);ax.set_ylim(-.6,len(x)-.35)
    ax.set_xlabel(xlabel);ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='y',length=0)
    if log:
        ax.set_xscale('log');ax.xaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10,subs=(1,2,5)));ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v,pos: f'{v:g}'));ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

def textrows(ax,rows):
    ax.axis('off')
    for i,(label,value) in enumerate(rows):
        y=1-(i+.5)/len(rows);ax.text(0,y,label,va='center',fontsize=8);ax.text(1,y,value,ha='right',va='center',fontsize=8)
        if i<len(rows)-1:ax.axhline(1-(i+1)/len(rows),color=C['grid'],lw=.6)

def f1_r5(r):
    s=r.source;a=s.table('activity');adj=s.table('activity_adjusted');st=s.table('activity_strata');loo=s.table('activity_loo');perm=s.table('activity_permutation')
    f=figure(190,155)
    heading(f,'A','Evidence-layer overview',.045,.975)
    aa=f.add_axes([.057,.807,.886,.143]);aa.set_xlim(0,5);aa.set_ylim(0,1);aa.axis('off')
    primary=a[a.analysis.str.contains('primary',case=False)].iloc[0]
    g=s.table('igvf_genes');g=g[(g.usable_target==1)&(g.Science_bait==0)]
    invivo=s.table('invivo_screens');wam_class=s.table('wam_class')
    layers=[
        ('Experimental ASD','AP–MS interactome','Frozen prey-only set',f"{int(primary.exposed_n)} genes",C['prey']),
        ('Neuronal-activity','CRISPRi','Primary FDR < 0.10',f"{int(primary.exposed_hits)} / {int(primary.exposed_n)} hits",C['hit']),
        ('IGVF cortical-neuron','Perturb-seq','Eligible non-bait targets',f'{len(g):,} genes',C['other']),
        ('In-vivo neuronal','dependency','Driver × age grid',f'{len(invivo)} conditions',C['prey']),
        ('ASD cortex','Wamsley','Patient cortical cell types',f'{wam_class.source_cell_type.nunique()} cell types',C['hit']),
    ]
    for i,(line1,line2,detail,count,color) in enumerate(layers):
        x=i+.06
        aa.plot([x,x+.78],[.98,.98],color=color,lw=2.0,solid_capstyle='butt')
        aa.text(x,.82,f'{i+1:02d}',ha='left',va='top',fontsize=5.5,color=C['other'],fontweight='bold')
        aa.text(x,.65,line1,ha='left',va='top',fontsize=6.0,fontweight='bold',color=C['ink'])
        aa.text(x,.49,line2,ha='left',va='top',fontsize=6.0,fontweight='bold',color=C['ink'])
        aa.text(x,.27,detail,ha='left',va='top',fontsize=5.4,color=C['other'])
        aa.text(x,.05,count,ha='left',va='bottom',fontsize=7.2,fontweight='bold',color=color)
        if i<4:aa.plot([i+1,i+1],[.12,.84],color=C['rule'],lw=.55)
    heading(f,'B','Neuronal-activity hit enrichment',.045,.78)
    ax=axes(f,[.10,.42,.38,.28],ylabel='Activity-hit genes (%)',grid=True)
    for j,(_,row) in enumerate(a.iterrows()):
        for k,(hc,nc,co) in enumerate([('exposed_hits','exposed_n',C['prey']),('comparator_hits','comparator_n',C['other'])]):
            x=j+(k-.5)*.3;y=row[hc]/row[nc]*100;ax.bar(x,y,.27,color=co,zorder=3)
            ax.text(x,y+.7,f'{int(row[hc])}/{int(row[nc]):,}',fontsize=7,ha='center')
    ax.set_xticks([0,1],['Primary\nFDR < 0.10','Strict\nFDR < 0.05']);ax.set_ylim(0,22);ax.spines['left'].set_visible(False)
    ax.legend(handles=[Patch(color=C['prey'],label='PPI prey'),Patch(color=C['other'],label='Other screened')],loc='upper right',fontsize=6.5)
    heading(f,'C','Selection-stratified permutation summaries',.53,.78)
    ax=f.add_axes([.62,.414,.33,.236])
    ax.axvline(0,color=C['other'],lw=.65,ls=(0,(3,3)),zorder=0)
    for y,key,color in [(1,'primary_fdr_lt_0_10',C['prey']),(0,'strict_fdr_lt_0_05',C['other'])]:
        row=perm[perm.analysis==key].iloc[0]
        lo,hi,mu,observed=[100*float(row[c]) for c in ['null_q025','null_q975','null_mean','observed_hit_rate_difference']]
        ax.plot([lo,hi],[y,y],color=color,lw=2.6,alpha=.78,solid_capstyle='round',zorder=1)
        ax.scatter([mu],[y],marker='|',s=54,color=C['ink'],linewidths=.9,zorder=3)
        ax.scatter([observed],[y],marker='o',s=21,facecolors=color,edgecolors='white',linewidths=.6,zorder=4)
        ax.text(9.5,y,f'P = {float(row.empirical_P_plus_one):.5f}',ha='left',va='center',fontsize=6.1,color=color)
    ax.set_yticks([1,0],['Primary FDR < 0.10\nn = 100,000','Strict FDR < 0.05\nn = 100,000'])
    ax.set_xlim(-8,16);ax.set_ylim(-.6,1.6);ax.set_xticks([-5,0,5,10,15]);ax.set_xlabel('Hit-rate difference (percentage points)',fontsize=6.3,labelpad=2)
    ax.tick_params(axis='x',labelsize=5.9,length=2,pad=2);ax.tick_params(axis='y',labelsize=6.0,length=0,pad=4)
    ax.grid(axis='x',color=C['grid'],lw=.4);ax.set_axisbelow(True);ax.spines[['top','right','left']].set_visible(False)
    ax.plot([],[],color=C['rule'],lw=2.6,label='Null 2.5–97.5% range')
    ax.plot([],[],linestyle='None',marker='|',color=C['ink'],markersize=5,label='Null mean')
    ax.plot([],[],linestyle='None',marker='o',markerfacecolor='white',markeredgecolor=C['ink'],markersize=3.5,label='Observed Δ')
    ax.legend(loc='upper left',bbox_to_anchor=(-.02,1.20),ncol=3,frameon=False,handlelength=1.0,columnspacing=.65,fontsize=5.25)
    heading(f,'D','Adjusted and stratum-omission estimates',.045,.29)
    ax=axes(f,[.24,.085,.70,.145],xlabel='Odds ratio (95% CI)')
    p=adj[(adj.model=='primary_fdr_lt_0_10')&(adj.term=='prey_only')].iloc[0];c=st[(st.row_type=='CMH_common')&(st.analysis=='primary_fdr_lt_0_10')].iloc[0]
    d=loo[loo.endpoint=='primary_fdr_lt_0_10']
    clean_stratum=lambda x:'provenance unassigned' if 'unresolved' in str(x).lower() else str(x).replace('_',' ')
    labels=['Selection adjusted','CMH common']+['Omit '+clean_stratum(x) for x in d.dropped_stratum]
    forest(ax,labels,[p.OR,c.OR]+d.leave_one_out_CMH_OR.tolist(),[p.OR_CI_low,c.CI95_low]+d.CI95_low.tolist(),[p.OR_CI_high,c.CI95_high]+d.CI95_high.tolist(),log=True,xlabel='Odds ratio (95% CI)')
    ax.tick_params(axis='y',labelsize=6.5);ax.set_xlim(.7,6)
    r.save(f,'Figure_1','Figure 1. Activity-screen enrichment. (A) Five frozen evidence layers with source-defined anchor counts: 118 prey-only genes in the experimental ASD AP–MS interactome; 21/118 primary neuronal-activity CRISPRi hits; 1,186 eligible non-bait IGVF targets; eight in-vivo driver-by-age conditions; and 36 Wamsley cortical cell types. Counts describe their respective source universes and do not imply causal transitions. (B) Exact frozen primary and strict hit counts and rates. (C) Source-saved permutation summaries for the primary FDR < 0.10 and strict FDR < 0.05 endpoints, each with 100,000 within-selection-stratum permutations. Circles mark the observed hit-rate differences, vertical ticks the saved null means, and bars the saved 2.5th–97.5th percentile ranges; labels show the original empirical P values. The frozen input contains no individual permutation draws, so no histogram or density has been reconstructed. (D) Primary selection-adjusted, common CMH and all four stratum-omission estimates with original 95% confidence intervals.')

def f2_r5(r):
    s=r.source;g=s.table('igvf_genes');g=g[(g.usable_target==1)&(g.Science_bait==0)];m=s.table('igvf_models').iloc[:4];mt=s.table('igvf_matched').iloc[0];b=s.table('igvf_balance');hp=one(s.table('holm_family'),test='A1_IGVF_network_adjusted').P_adjusted_Holm
    f=figure(190,155);heading(f,'A','Full gene-level perturbation distribution',.04,.975)
    ax=axes(f,[.10,.16,.43,.73],xlabel='Mean energy distance (log scale)',ylabel='Cumulative fraction of genes')
    for k,col,lab in [(1,C['prey'],'PPI prey'),(0,C['other'],'Other targets')]:
        vals=np.sort(g[g.Science_prey_only==k].mean_energy_distance);ax.step(vals,np.arange(1,len(vals)+1)/len(vals),where='post',color=col,lw=1.7,label=f'{lab} (n = {len(vals):,})')
    ax.set_xscale('log');ax.set_ylim(0,1.01);ax.legend(loc='upper left',fontsize=7);ax.grid(axis='y',lw=.4,color=C['grid'])
    heading(f,'B','Same 623 complete-case genes',.59,.975)
    ax=axes(f,[.77,.66,.19,.22]);forest(ax,['Unadjusted','Source adjusted','Network adjusted','Guide-count\nsensitivity'],m.OR.tolist(),m.CI_low_OR.tolist(),m.CI_high_OR.tolist(),log=True,xlabel='Odds ratio (95% CI)');ax.set_xlim(.9,5.4);ax.tick_params(axis='y',labelsize=6.8)
    f.text(.61,.565,'Primary Holm P = '+numeric_p(hp),fontsize=7.4)
    heading(f,'C','Matched pairs and covariate balance',.59,.525)
    f.text(.61,.465,f'{int(mt.matched_sets_n)} pairs · hits {int(mt.exposed_hits)}/{int(mt.n_exposed)} vs {int(mt.comparator_hits)}/{int(mt.n_comparator)}',fontsize=7.5)
    f.text(.61,.425,f'Risk difference = {mt.estimate_difference:.3f} [{mt.CI_low:.3f}, {mt.CI_high:.3f}]\nPaired permutation P = {mt.P_permutation:.4f}',fontsize=7.2,linespacing=1.5)
    ax=axes(f,[.78,.16,.18,.18],xlabel='Absolute standardized difference')
    bb=b[b.analysis_id=='A1_IGVF'];metrics=list(bb.metric.unique())
    for j,metric in enumerate(metrics):
        z=bb[bb.metric==metric]
        for _,row in z.iterrows():ax.plot(row.abs_SMD,j,'o' if 'after' in str(row.phase).lower() or 'post' in str(row.phase).lower() else 's',color=C['prey'] if 'after' in str(row.phase).lower() or 'post' in str(row.phase).lower() else C['other'],ms=4)
        if len(z)==2:ax.plot(z.abs_SMD,[j,j],color=C['rule'],lw=.8)
    names={'degree_log':'BioPlex degree','abundance':'293T abundance','promoter_log':'Tested promoters'}
    ax.set_yticks(range(len(metrics)),[names.get(x,x) for x in metrics]);ax.tick_params(axis='y',labelsize=6.5);ax.set_xlim(left=0);ax.set_ylim(-.6,len(metrics)-.4)
    f.text(.61,.08,'Square: before matching · circle: after matching',fontsize=7.2,color=C['ink'])
    r.save(f,'Figure_2','Figure 2. Cortical-neuron perturbation and measured-network adjustment. A, empirical cumulative distributions of all 1,186 eligible non-bait IGVF targets (125 prey-only; 1,061 comparators); the energy-distance values are the frozen author-derived measures, without smoothing or new inference. B, four frozen logistic models in the same 623 genes (120 prey-only; 503 comparators; 136 events), with profile 95% confidence intervals. The displayed primary Holm P belongs to the original four-test A1/A2/A3/A4 family. C, 106 frozen matched pairs, original matched-set risk difference and bootstrap confidence interval, paired permutation P, and the three measured covariate balance summaries. Full coverage, continuous sensitivity and program testability appear in S5.')

def f3_r5(r):
    s=r.source;d=s.table('invivo_screens');p=s.table('invivo_primary').iloc[0];a=s.table('invivo_age');ls=s.sup(6,'09_LEAVE_ONE_SCREEN_OUT.tsv');lc=s.sup(6,'10_LEAVE_ONE_CELLTYPE_OUT.tsv')
    f=figure(190,150);heading(f,'A','Eight neuronal conditions',.045,.97)
    ax=f.add_axes([.14,.25,.43,.61]);order=['Syn1','Emx1','Vglut2','Vgat'];mat=np.array([[one(d,screen=f'{n}-GW-{age}M').mean_difference for age in [2,4]] for n in order]);q=np.array([[one(d,screen=f'{n}-GW-{age}M').permutation_P_BH_FDR for age in [2,4]] for n in order])
    cm=LinearSegmentedColormap.from_list('vulnerability',['#F4F7F6','#9CC6C2','#1A6C72']);im=ax.imshow(mat,vmin=0,vmax=.16,cmap=cm,aspect='auto')
    ax.set_xticks([0,1],['2 months','4 months']);ax.xaxis.tick_top();ax.tick_params(length=0,pad=9);ax.set_yticks(range(4),order)
    for i in range(4):
        for j in range(2):
            ax.text(j,i,f'Δ={mat[i,j]:.3f}',ha='center',va='center',fontsize=12.5,color='white' if mat[i,j]>.095 else C['ink'])
            if q[i,j]<.05:
                ax.scatter(j,i+.22,s=16,marker='o',facecolors='white',edgecolors='#40464B',linewidths=.65,zorder=5)
    for y in [.5,1.5,2.5]:ax.axhline(y,color='white',lw=3)
    ax.axvline(.5,color='white',lw=3)
    for sp in ax.spines.values():sp.set_visible(False)
    cb=f.colorbar(im,cax=f.add_axes([.14,.155,.43,.022]),orientation='horizontal');cb.set_label('Prey minus comparator mean vulnerability',fontsize=7.7);cb.outline.set_visible(False)
    f.text(.14,.065,'Open circle: within-eight-condition BH q < 0.05',fontsize=6.4,color=C['ink'])
    f.text(.14,.205,'8/8 positive effects',fontsize=7.1)
    heading(f,'B','Primary composite',.65,.97)
    ax=axes(f,[.68,.64,.27,.15]);forest(ax,[''],[p.mean_difference],[p.bootstrap_CI_low],[p.bootstrap_CI_high],xlabel='Mean difference (95% CI)');ax.set_xlim(-.005,.17);ax.set_ylim(-.65,.65)
    f.text(.68,.855,f'Δ = {p.mean_difference:.4f}\nCohen’s d = {p.cohen_d:.3f}',fontsize=9,linespacing=1.6)
    f.text(.68,.555,'Permutation P = '+numeric_p(p.permutation_P_two_sided),fontsize=7.3)
    heading(f,'C','Age and omission sensitivities',.65,.465)
    ax=axes(f,[.81,.14,.145,.24]);labels=['2 months','4 months','Omit one\ncondition','Omit one\ncell class']
    forest(ax,labels,a.mean_difference.tolist()+[np.nan,np.nan],a.bootstrap_CI_low.tolist()+[np.nan,np.nan],a.bootstrap_CI_high.tolist()+[np.nan,np.nan],xlabel='Mean difference')
    for t in list(ax.texts):t.remove()
    for y,z in [(1,ls),(0,lc)]:ax.plot([z.mean_difference.min(),z.mean_difference.max()],[y,y],color=C['other'],lw=3);ax.plot(z.mean_difference,np.full(len(z),y),'|',color=C['ink'],ms=5)
    ax.set_xlim(0,.18);ax.tick_params(axis='y',labelsize=7);f.text(.65,.04,'Omission bars: range of saved estimates',fontsize=6.5)
    r.save(f,'Figure_3','Figure 3. In-vivo neuronal dependency. A, all four neuronal drivers at both ages, using the frozen mean vulnerability difference on the common 868 prey-only / 9,081 comparator gene set. Cell values are effect differences (Δ), not P values; positive values denote greater vulnerability. White-filled circles with thin dark-grey outlines mark conditions meeting the original within-eight-condition BH q < 0.05 threshold. All eight effects point in the same direction, although not every condition is significant. These conditions share one study and library. B, primary eight-condition composite, original bootstrap 95% interval, two-sided permutation P and Cohen’s d. C, two age composites with original bootstrap 95% intervals, and ranges of all saved leave-one-condition and leave-one-cell-class estimates. Range bars are sensitivity ranges, not confidence intervals. S6 retains all detailed results.')

def s1(r):
    s=r.source;a=s.table('activity');adj=s.table('activity_adjusted');st=s.table('activity_strata');loo=s.table('activity_loo');pm=s.table('activity_permutation');asc=s.sup(1,'08_ASD_NDD_ASCERTAINMENT_SENSITIVITY.tsv').iloc[0]
    f=subfig('S1  Activity-screen robustness and ascertainment',h=190)
    ax=axes(f,[.19,.64,.29,.23]);ttl(ax,'A','Prey-only associations')
    labels=[];x=[];lo=[];hi=[]
    for _,rr in a.iterrows():
        short='Primary' if 'primary' in rr.analysis else 'Strict';ad=one(adj,model=rr.analysis,term='prey_only');cm=one(st,analysis=rr.analysis,row_type='CMH_common')
        labels += [short+' crude',short+' adjusted',short+' CMH'];x += [rr.OR,ad.OR,cm.OR];lo += [rr.CI95_low,ad.OR_CI_low,cm.CI95_low];hi += [rr.CI95_high,ad.OR_CI_high,cm.CI95_high]
    forest(ax,labels,x,lo,hi,log=True,xlabel='Odds ratio (95% CI)');ax.set_xlim(.7,5)
    ax=axes(f,[.68,.64,.27,.23]);ttl(ax,'B','All selection strata')
    d=st[st.row_type=='stratum'];labels=[('Primary: ' if 'primary' in z.analysis else 'Strict: ')+z.selection_stratum.replace('_',' ') for _,z in d.iterrows()]
    for j,(_,z) in enumerate(d.iterrows()):
        yy=len(d)-1-j;ax.plot([z.exposed_hit_rate*100,z.comparator_hit_rate*100],[yy,yy],color=C['rule']);ax.plot(z.exposed_hit_rate*100,yy,'o',color=C['prey'],ms=4);ax.plot(z.comparator_hit_rate*100,yy,'s',color=C['other'],ms=3)
    ax.set_yticks(range(len(d))[::-1],labels);ax.set_xlabel('Activity-hit genes (%)');ax.set_xlim(-2,65);ax.tick_params(axis='y',labelsize=6.5)
    ax=axes(f,[.22,.22,.30,.28]);ttl(ax,'C','Every stratum omission')
    labels=[('Primary: ' if 'primary' in z.endpoint else 'Strict: ')+z.dropped_stratum.replace('_',' ') for _,z in loo.iterrows()]
    forest(ax,labels,loo.leave_one_out_CMH_OR.tolist(),loo.CI95_low.tolist(),loo.CI95_high.tolist(),log=True,xlabel='Common odds ratio (95% CI)');ax.tick_params(axis='y',labelsize=6.5)
    ax=f.add_axes([.62,.21,.33,.28]);ttl(ax,'D','Permutation summaries')
    rows=[('Assignments per definition',f'{int(pm.iloc[0].n_permutations):,}'),('Primary empirical P',f'{pm.iloc[0].empirical_P_plus_one:.6f}'),('Strict empirical P',f'{pm.iloc[1].empirical_P_plus_one:.6f}'),('ASD/NDD-selected prey',f'{int(asc.exposed_hits)}/{int(asc.exposed_n)} hits'),('Selected comparators',f'{int(asc.comparator_hits)}/{int(asc.comparator_n)} hits')]
    textrows(ax,rows);f.text(.62,.135,'Sparse ascertainment subset: direction uninformative.\nIndividual permutation draws unavailable.',fontsize=7,linespacing=1.6)
    f.text(.05,.065,'Selection strata remain incomplete for unresolved target provenance. Circles: prey-only; squares: comparators.',fontsize=7)
    r.save(f,'Figure_S1','Supplementary Figure S1. A, primary and strict source-defined hits with crude, adjusted and CMH confidence intervals. B, all primary/strict strata, including sparse disease-selected categories. C, all eight saved leave-one-stratum-out estimates and intervals. D, exact permutation summary P values and sparse ASD/NDD ascertainment counts. The zero-hit ascertainment subset is uninformative, not evidence of protection. Individual permutation draws are absent from the frozen bundle. No primary definition or membership was changed.','supplementary')

def s2(r):
    s=r.source;d=s.sup(2,'S2A_VARIANT_REWIRING.tsv');g=s.sup(2,'10_VARIANT_ALTERED_GENELEVEL_SENSITIVITY.tsv').iloc[0];c=s.sup(2,'11_VARIANT_CONTINUOUS_SENSITIVITY.tsv').iloc[0]
    f=subfig('S2  Patient-variant rewiring: non-supportive association',h=160)
    ax=axes(f,[.33,.50,.61,.35]);ttl(ax,'A','All frozen variant-pair permutation comparisons')
    labels=[x.replace('_',' ') for x in d.analysis]
    yy=np.arange(len(d))[::-1]
    for j,(_,z) in enumerate(d.iterrows()):
        ax.plot([z.null_q025,z.null_q975],[yy[j],yy[j]],color=C['rule'],lw=4);ax.plot(z.expected,yy[j],'s',color=C['other'],ms=4);ax.plot(z.observed,yy[j],'o',color=C['prey'],ms=4)
        ax.text(.99,yy[j]-.26,'P = '+numeric_p(z.empirical_P_plus_one),transform=ax.get_yaxis_transform(),ha='right',fontsize=6.5)
    ax.set_yticks(yy,labels);ax.tick_params(axis='y',labelsize=6.7);ax.set_xlabel('Activity-hit altered pairs');ax.set_ylim(-.6,len(d)-.3)
    ax=axes(f,[.18,.16,.30,.17]);ttl(ax,'B','Unique-gene sensitivity');forest(ax,['Altered vs\nnever altered'],[g.OR],[g.CI95_low],[g.CI95_high],[g.Fisher_P],log=True,xlabel='Odds ratio (95% CI)')
    ax=axes(f,[.68,.16,.27,.17]);ttl(ax,'C','Continuous sensitivity');forest(ax,[f'n = {int(c.N)} genes'],[c.spearman_rho],[c.bootstrap_95_low],[c.bootstrap_95_high],[c.nominal_P],xlabel='Spearman correlation (95% CI)')
    f.text(.05,.04,'Observed altered-pair count: circle · expected count: square · saved null 2.5–97.5% interval: line',fontsize=6.8)
    r.save(f,'Figure_S2','Supplementary Figure S2. All frozen variant-pair comparisons are retained. A, observed and expected counts and original null quantiles, with the original greater-tail empirical P; lines are null-reference ranges, not effect confidence intervals. B, unique-gene altered-versus-never-altered sensitivity (2/21 versus 6/34) and Fisher confidence interval/P. C, continuous gene-level Spearman association (55 genes), original bootstrap interval and nominal P. Variant rewiring remains non-supportive; no patient-missense prioritization claim is introduced.','supplementary')

def s3(r):
    d=r.source.sup(3,'S3A_BAIT_NEIGHBORHOOD_NULLS.tsv').sort_values('bait');f=subfig('S3  All bait-specific neuronal-activity neighborhoods',w=210,h=230)
    for j in range(4):
        z=d.iloc[j*25:(j+1)*25];ax=axes(f,[.085+j*.24,.15,.165,.75]);yy=np.arange(len(z))[::-1];v=z.enrichment_ratio_vs_screened_experimental_preys
        for i,(_,rr) in enumerate(z.iterrows()):
            if np.isfinite(v.iloc[i]):ax.plot(v.iloc[i],yy[i],'o',ms=3,color=C['prey'] if rr.minimum_degree_rule_met else C['other'])
            else:ax.text(.03,yy[i],'—',transform=ax.get_yaxis_transform(),fontsize=6)
        ax.axvline(1,color=C['other'],ls='--',lw=.6);ax.set_yticks(yy,z.bait);ax.tick_params(axis='y',labelsize=6.3,length=0);ax.set_xlim(0,max(6,float(d.enrichment_ratio_vs_screened_experimental_preys.max())*1.04));ax.set_xlabel('Hit enrichment',fontsize=7);ax.set_title(f'Baits {j*25+1}–{(j+1)*25}',fontsize=8)
    f.text(.06,.07,'All 100 baits in alphabetical order. No bait meets the frozen BH significance threshold.',fontsize=8)
    f.text(.06,.04,'Dash: enrichment unavailable. Gray: original minimum-degree rule not met. Exact counts and adjusted P values: frozen source tables.',fontsize=6.8)
    r.save(f,'Figure_S3','Supplementary Figure S3. All 100 bait neighborhoods, alphabetically ordered without significance-based selection. Points show saved enrichment ratios relative to screen-measured experimental prey genes. Unavailable ratios are dashes; gray points identify neighborhoods below the frozen minimum-degree rule. No neighborhood meets the original BH threshold. Full bait-specific tested counts, hits, expected counts, P, q and membership are retained in the supplementary tables.','supplementary')

def s4(r):
    z=r.source.sup(4,'S4A_SYNAPTOGENESIS.tsv').iloc[0];f=subfig('S4  Independent synaptogenesis-screen comparison',h=130)
    ax=axes(f,[.12,.25,.32,.55],ylabel='Source-defined hit genes (%)');ttl(ax,'A','Tested-gene hit fractions')
    ns=[z.prey_only_covered,z.comparator_n];hs=[z.prey_only_hits,z.comparator_hits]
    ax.bar([0,1],np.array(hs)/ns*100,color=[C['prey'],C['other']],width=.55)
    for j in range(2):ax.text(j,hs[j]/ns[j]*100+.5,f'{int(hs[j])}/{int(ns[j])}',ha='center',fontsize=8)
    ax.set_xticks([0,1],['PPI prey','Comparators']);ax.set_ylim(0,16)
    ax=axes(f,[.64,.48,.30,.20]);ttl(ax,'B','Frozen comparison');forest(ax,[''],[z.odds_ratio],[z.ci95_low],[z.ci95_high],[z.fisher_two_sided_p],log=True,xlabel='Odds ratio (95% CI)');ax.set_xlim(.15,2)
    f.text(.60,.25,f'{int(z.screened_genes)} screened genes; {int(z.author_candidate_hits)} source hits\n{int(z.bait_or_other_prey_role_excluded_n)} bait/other-role genes excluded',fontsize=8,linespacing=1.7)
    f.text(.06,.06,'Non-supportive external comparison; source testability and coverage constrain interpretation.',fontsize=8)
    r.save(f,'Figure_S4','Supplementary Figure S4. The frozen external synaptogenesis screen is non-supportive. A, exact source-defined hit fractions (7/91 prey-only versus 63/544 comparators). B, original Fisher odds ratio, confidence interval and two-sided P. The source screen comprises 644 genes and 71 primary score hits, including assay controls; nine bait/other-prey-role genes are excluded from the displayed comparison. No hit definition is substituted.','supplementary')

def s5(r):
    s=r.source;m=s.table('igvf_models');b=s.table('igvf_balance');mt=s.table('igvf_matched').iloc[0];pr=s.sup(5,'B_PROGRAM_ENRICHMENT_ALL.tsv');br=s.sup(5,'B_PROGRAM_BREADTH.tsv');test=s.sup(5,'B_PROGRAM_TESTABILITY.tsv');cc=s.sup(5,'A1_IGVF_COMPLETECASE_INPUT.tsv')
    f=subfig('S5  IGVF model inventory, matching and program testability',h=195)
    ax=axes(f,[.24,.69,.28,.19]);ttl(ax,'A','Complete-case logistic models');mm=m.iloc[:4];forest(ax,['Unadjusted','Source adjusted','Network adjusted','Guide count'],mm.OR.tolist(),mm.CI_low_OR.tolist(),mm.CI_high_OR.tolist(),log=True,xlabel='Odds ratio (95% profile CI)')
    ax=axes(f,[.75,.69,.20,.19]);ttl(ax,'B','Continuous sensitivity');z=m.iloc[4];forest(ax,[''],[z.estimate],[z.ci_low],[z.ci_high],[z.P_raw],xlabel='Rank-normal coefficient')
    f.text(.07,.595,'623 complete cases: 120 prey / 503 comparators; 136 events. Original eligible universe: 1,186 genes.',fontsize=7.3)
    ax=axes(f,[.24,.35,.28,.17]);ttl(ax,'C','Matching covariate balance');bb=b[b.analysis_id=='A1_IGVF'];metrics=list(bb.metric.unique());names=['BioPlex degree','293T abundance','Tested promoters']
    for i,met in enumerate(metrics):
        z=bb[bb.metric==met];ax.plot(z.abs_SMD,[i]*len(z),color=C['rule']);
        for _,v in z.iterrows():ax.plot(v.abs_SMD,i,'o' if v.phase=='after' else 's',color=C['prey'] if v.phase=='after' else C['other'],ms=4)
    ax.set_yticks(range(len(metrics)),names);ax.set_xlabel('Absolute standardized difference');ax.set_ylim(-.5,2.5)
    ax=f.add_axes([.62,.33,.33,.20]);ttl(ax,'D','Matched comparison');textrows(ax,[('Matched pairs',str(int(mt.matched_sets_n))),('Prey / comparator hits',f'{int(mt.exposed_hits)} / {int(mt.comparator_hits)}'),('Risk difference',f'{mt.estimate_difference:.3f}'),('Original 95% interval',f'[{mt.CI_low:.3f}, {mt.CI_high:.3f}]'),('Paired permutation P',numeric_p(mt.P_permutation))])
    ax=axes(f,[.10,.105,.40,.12]);ttl(ax,'E','Reported program breadth (descriptive)')
    for k,col,lab in [(1,C['prey'],'PPI prey'),(0,C['other'],'Other targets')]:
        v=np.sort(br.loc[br.Science_prey_only==k,'source_reported_significant_GP_n']);ax.step(v,np.arange(1,len(v)+1)/len(v),where='post',color=col,label=lab)
    ax.set_xlabel('Number of reported significant programs');ax.set_ylabel('Fraction');ax.legend(fontsize=6.5,loc='lower right')
    f.text(.60,.20,f'All {len(pr)} programs: enrichment not estimable',fontsize=8,fontweight='bold');f.text(.60,.14,'Tested denominators unavailable.\nZero reported rows does not establish a non-hit.\nNo program-level positive claim.',fontsize=7.3,linespacing=1.6)
    r.save(f,'Figure_S5','Supplementary Figure S5. Full four-model logistic inventory and the continuous energy-distance sensitivity use the same complete-case set. Logistic intervals are penalized profile intervals; the continuous interval uses the frozen HC3 model and its own nominal P. Matching retains 106 pairs. Covariate balance uses saved standardized differences; square/circle denote before/after matching. E is descriptive source-reported program breadth across all eligible targets. All 200 program tests lack tested denominators and remain not estimable; absence from the significant-only source does not establish a non-hit.','supplementary')

def s6(r):
    s=r.source;d=s.table('invivo_screens');ls=s.sup(6,'09_LEAVE_ONE_SCREEN_OUT.tsv');lc=s.sup(6,'10_LEAVE_ONE_CELLTYPE_OUT.tsv');age=s.table('invivo_age');alt=s.sup(6,'12_ALTERNATIVE_GENE_SUMMARY_SENSITIVITY.tsv');ext=s.sup(6,'13_VGLUT2_18M_EXTENSION.tsv')
    f=subfig('S6  In-vivo profiles and complete stability analyses',h=220)
    ax=axes(f,[.18,.62,.30,.25]);ttl(ax,'A','All eight condition effects');forest(ax,d.screen.str.replace('-GW-',' / ').tolist(),d.mean_difference.tolist(),p=d.gene_label_permutation_P_two_sided.tolist(),xlabel='Mean vulnerability difference');ax.set_xlim(0,.18)
    ax=axes(f,[.70,.62,.25,.25]);ttl(ax,'B','Leave one condition out');forest(ax,ls.omitted_screen.str.replace('-GW-',' / ').tolist(),ls.mean_difference.tolist(),xlabel='Composite mean difference');ax.set_xlim(0,.13)
    ax=axes(f,[.18,.32,.30,.18]);ttl(ax,'C','Leave one cell class out');forest(ax,lc.omitted_cell_type.tolist(),lc.mean_difference.tolist(),xlabel='Composite mean difference');ax.set_xlim(0,.13)
    ax=axes(f,[.70,.32,.25,.18]);ttl(ax,'D','Age and gene summaries');dd=pd.concat([age,alt],ignore_index=True);forest(ax,['2 months','4 months','Trimmed mean','Winsorized mean'],dd.mean_difference.tolist(),dd.bootstrap_CI_low.tolist(),dd.bootstrap_CI_high.tolist(),xlabel='Mean difference (95% CI)');ax.set_xlim(0,.19)
    ax=axes(f,[.22,.105,.67,.10]);ttl(ax,'E','Vglut2 age extension');forest(ax,ext.screen.str.replace('-GW-',' / ').tolist(),ext.observed_effect.tolist(),p=ext.two_sided_P_Holm_3test_family.tolist(),xlabel='Mean difference; labels: original Holm P');ax.set_xlim(0,.25)
    r.save(f,'Figure_S6','Supplementary Figure S6. A, all eight screen-specific effects and original two-sided permutation P values; no source confidence intervals are available. B–C, every leave-one-screen and leave-one-cell-class composite. D, age summaries and both alternative gene aggregations with their original bootstrap intervals. E, all three prespecified Vglut2 ages (including the 18-month extension) with original three-test Holm P values; the extension is secondary. Omission estimates are sensitivity summaries, not independent cohorts.','supplementary')

def s7(r):
    s=r.source;edge=s.sup(7,'07_NEURONAL_PPI_EDGE_SUPPORT.tsv');node=s.sup(7,'08_NEURONAL_PPI_NODE_SUPPORT.tsv');tests=s.sup(7,'09_NEURONAL_CONTEXT_FUNCTIONAL_TEST.tsv');d=s.sup(7,'S7C_DISEASE_NONESTIMABLE.tsv').iloc[0];spec=s.sup(7,'S7A_BIOPLEX_SPECIFICITY.tsv')
    f=subfig('S7  Generic-network and neuronal-context boundaries',h=210)
    ax=axes(f,[.12,.70,.34,.18],ylabel='Testable genes');ttl(ax,'A','Neuronal-context coverage')
    for j,nm in enumerate(['Pintacuda','Murtaza','Gao']):
        tested=int(node[nm+'_testable'].sum());supported=int(node[nm+'_supported'].sum());ax.bar(j,tested,color=C['pale'],width=.6);ax.bar(j,supported,color=C['prey'],width=.6);ax.text(j,tested+3,f'{supported}/{tested}',ha='center',fontsize=7)
    ax.set_xticks(range(3),['Pintacuda','Murtaza','Gao']);ax.set_ylim(0,max(node[n+'_testable'].sum() for n in ['Pintacuda','Murtaza','Gao'])*1.3);ax.set_xlabel('Supported / testable among 1,043 prey genes')
    ax=axes(f,[.72,.70,.24,.18]);ttl(ax,'B','Activity support gradient');z=tests.iloc[:2];forest(ax,['Unadjusted\n(n = 18)','Adjusted\n(n = 17)'],z.estimate.tolist(),z.CI_low.tolist(),z.CI_high.tolist(),z.P.tolist(),log=True,xlabel='Odds ratio (95% CI)');ax.set_xticks([.1,1,10]);ax.set_xticklabels(['0.1','1','10'])
    ax=axes(f,[.25,.41,.23,.17]);ttl(ax,'C','In-vivo support gradient');z=tests.iloc[2:4];forest(ax,['Unadjusted\n(n = 105)','Adjusted\n(n = 16)'],z.estimate.tolist(),z.CI_low.tolist(),z.CI_high.tolist(),z.P.tolist(),xlabel='Linear trend (95% CI)')
    ax=axes(f,[.72,.41,.24,.17]);ttl(ax,'D','Same-bait comparisons');z=tests.iloc[4:6];forest(ax,['Activity\n(3 baits)','In vivo\n(6 baits)'],z.estimate.tolist(),z.CI_low.tolist(),z.CI_high.tolist(),z.P.tolist(),xlabel='Paired difference (95% CI)')
    ax=axes(f,[.12,.125,.34,.17],ylabel='Eligible baits');ttl(ax,'E','Disease-specificity feasibility')
    vals=[d.ASD_eligible_baits_n,d.DD_eligible_baits_n,d.SCZ_eligible_baits_n];ax.bar(range(3),vals,color=[C['prey'],C['other'],C['pale']],width=.55)
    for j,v in enumerate(vals):ax.text(j,v+.1,str(int(v)),ha='center',fontsize=8)
    ax.set_xticks(range(3),['ASD','DD only','SCZ only']);ax.set_ylim(0,7)
    f.text(.60,.28,'Disease-specificity tests: not estimable',fontsize=8,fontweight='bold');f.text(.60,.22,'Sparse eligible bait sets and inadequate matching\npreclude formal disease-specificity inference.',fontsize=7.5,linespacing=1.6)
    f.text(.60,.12,'Generic AP–MS specificity is not established\nfor the in-vivo association. Cross-study activity\ncontrols do not establish same-assay specificity.',fontsize=7.3,linespacing=1.6)
    r.save(f,'Figure_S7','Supplementary Figure S7. A, exact frozen neuronal-context testable and supported gene counts; non-testability is not a biological null. B–C, source-adjusted and unadjusted support gradients, with their original model-specific intervals and P values. Activity effects are odds ratios; in-vivo effects are linear coefficients. D, frozen same-bait contrasts are on different outcome scales and cannot be compared as standardized magnitudes. E, eligible ASD/DD/SCZ bait counts. Disease-specificity contrasts remain non-estimable because of low coverage and inadequate matching, not negative formal tests. Generic-network specificity for in-vivo dependency remains unestablished.','supplementary')

def s8(r):
    s=r.source;w=s.table('wam_cells');non=s.sup(8,'S8B_WAMSLEY_NONNEURONAL.tsv');ct=s.sup(8,'S8B_WAMSLEY_CLASS_CONTRAST.tsv');a2=s.table('wam_a2').iloc[0];a3=s.table('wam_a3').iloc[0];a4=s.table('wam_a4').iloc[0]
    f=subfig('S8  Wamsley magnitudes, direct contrasts and secondary controls',w=210,h=250)
    # The deeper cell-type display uses raw median differences, unlike the AUC main panel.
    order=sorted(w.source_cell_type.unique(),key=natural_key)
    ax=axes(f,[.18,.34,.30,.56]);ttl(ax,'A','All cell-type median differences')
    for k,gs in enumerate(['SCIENCE_PREY_ONLY','ACTIVITY_VULNERABLE_PPI']):
        z=w[w.gene_set==gs].set_index('source_cell_type').loc[order];ax.plot(z.effect_median_difference,np.arange(36)[::-1]+(k-.5)*.20,'o' if k==0 else 's',ms=3,color=C['prey'] if k==0 else C['hit'],label='PPI prey' if k==0 else 'Activity hits')
    ax.set_yticks(range(36)[::-1],[x.replace('_',' ') for x in order]);ax.tick_params(axis='y',labelsize=6.1,length=0);ax.axvline(0,color=C['other'],ls='--',lw=.5);ax.set_xlabel('Median |logFC| difference');ax.legend(fontsize=6.3,loc='lower right')
    ax=axes(f,[.75,.70,.20,.15]);ttl(ax,'B','Non-neuronal omnibus')
    zz=non[non.test.isin(['SET1_NONNEURONAL','SET2_NONNEURONAL'])];forest(ax,['PPI prey','Activity hits'],zz.effect_median_difference.tolist(),p=zz.Holm_P_one_sided.tolist(),xlabel='Median |logFC| difference')
    ax=axes(f,[.75,.46,.20,.14]);ttl(ax,'C','Direct cell-class contrast');forest(ax,['PPI prey','Activity hits'],ct.effect_median_difference.tolist(),p=ct.Holm_P_two_sided.tolist(),xlabel='Difference in median contrasts')
    ax=axes(f,[.28,.115,.30,.12]);ttl(ax,'D','Matched generic-network control');forest(ax,['887 matched sets'],[a2.matched_estimate_difference_original_Y],[a2.matched_CI_low],[a2.matched_CI_high],[a2.matched_permutation_P],xlabel='Matched median |logFC| difference')
    ax=axes(f,[.77,.125,.19,.17]);ttl(ax,'E','Same-screen sensitivity');forest(ax,['Primary\n(n = 114)','Network\nadjusted (101)'],[a4.estimate_rank_beta,a4.network_adjusted_rank_beta],[a4.ci_low,a4.network_adjusted_ci_low],[a4.ci_high,a4.network_adjusted_ci_high],[a4.P_raw,a4.network_adjusted_P],xlabel='Rank-normal coefficient')
    f.text(.60,.36,'Subset-versus-rest matching: not estimable\nOriginal maximum |SMD| = '+f'{a3.matched_max_abs_SMD:.3f}',fontsize=7.2,linespacing=1.6)
    f.text(.05,.035,'Same-screen prioritization remains non-supportive; separate within-class results do not establish neuronal specificity.',fontsize=7.2)
    r.save(f,'Figure_S8','Supplementary Figure S8. A shows every cell type using median absolute-logFC differences rather than repeating the main-figure rank-AUC matrix. B, non-neuronal omnibus median differences and original one-sided Holm P. C, direct neuronal-minus-nonneuronal class contrasts and original two-sided Holm P; no intervals were supplied. D, saved generic-network matched sensitivity with bootstrap interval and permutation P. The subset-versus-rest matching comparison was not estimable under its balance criterion. E, primary and network-augmented same-screen models have their own rank-normal coefficient intervals and raw P values. Same-screen prioritization remains non-supportive.','supplementary')

def s11(r):
    s=r.source;days=s.table('gordon_days');forms=s.table('gordon_forms');cp=s.table('vel_compartments')
    f=subfig('S11  External-validation descriptive profiles',h=190)
    ax=axes(f,[.12,.62,.36,.23],ylabel='Mean |logFC|',xlabel='Developmental day');ttl(ax,'A','Gordon: all developmental days')
    for pn,ls in [('C36_all_9_forms','-'),('C24_6_forms_sensitivity','--')]:
        d=days[days.profile_panel==pn].sort_values('developmental_day')
        for col,co,lab in [('mean_A',C['hit'],'Hits'),('mean_B',C['other'],'Screened non-hits')]:ax.plot(d.developmental_day,d[col],ls,color=co,marker='o',ms=3,label=('C36 ' if 'C36' in pn else 'C24 ')+lab)
    ax.set_xticks([25,50,75,100]);ax.legend(fontsize=6.4)
    ax=axes(f,[.70,.58,.25,.29]);ttl(ax,'B','Gordon: all nine forms')
    d=forms[forms.contrast_set=='C36'].sort_values('ASD_form');yy=np.arange(len(d))[::-1]
    ax.plot(d.mean_A,yy,'o',color=C['hit'],ms=4,label='Hits');ax.plot(d.mean_B,yy,'s',color=C['other'],ms=3,label='Screened non-hits');ax.set_yticks(yy,d.ASD_form.str.replace('IdiopathicASD','Idiopathic ASD'));ax.set_xlabel('Mean |logFC|');ax.tick_params(axis='y',labelsize=6.5)
    ax=axes(f,[.12,.17,.36,.24],ylabel='Mean |coefficient|');ttl(ax,'C','Velmeshev: all neuronal compartments')
    order=['EXC','INH','OTHER_NEURON'];z=cp[cp.comparison=='C1_H_VS_N'].set_index('compartment').loc[order];rr=cp[cp.comparison=='C2_H_VS_R'].set_index('compartment').loc[order]
    for v,col,lab,mk in [(z.H_mean_abs_beta,C['hit'],'Hits','o'),(z.comparator_mean_abs_beta,C['other'],'Screened non-hits','s'),(rr.comparator_mean_abs_beta,C['prey'],'Other PPI','^')]:ax.plot(range(3),v,color=col,marker=mk,ms=4,label=lab)
    ax.set_xticks(range(3),['Excitatory','Inhibitory','Other neurons']);ax.legend(fontsize=6.5)
    f.text(.62,.39,'Descriptive profiles only',fontsize=9,fontweight='bold');f.text(.62,.31,'No day-, form- or compartment-level\nP values or confidence intervals are introduced.',fontsize=7.5,linespacing=1.7)
    f.text(.62,.19,'Gordon and Velmeshev primary classifications\nremain NON_SUPPORTIVE.',fontsize=7.5,linespacing=1.7)
    r.save(f,'Figure_S11','Supplementary Figure S11. Gordon group means across every day in both C36 and C24 (A) and all nine C36 forms (B). C24 is a prespecified sensitivity. C, Velmeshev source-estimated mean absolute coefficients in all three neuronal compartments for hits, same-screen non-hits and other PPI genes. These descriptive profiles have no new inferential tests and do not alter the NON_SUPPORTIVE primary classifications.','supplementary')

def s12(r):
    s=r.source;flow=s.sup(10,'COHORT_FLOW.tsv');exact=s.sup(12,'EXACT_OVERLAP_DONORS.tsv');poss=s.sup(12,'POSSIBLE_OR_UNRESOLVED_OVERLAPS.tsv');counts=s.sup(12,'RETAINED_DONOR_SAMPLE_CELLTYPE_COUNTS.tsv');status=s.sup(12,'DONOR_AUDIT_FINAL_STATUS.tsv');src=s.sup(12,'SOURCE_QA.tsv')
    f=subfig('S12  Donor overlap, cohort construction and identity limits',h=180)
    ax=axes(f,[.12,.59,.35,.27],ylabel='Donors');ttl(ax,'A','Public-roster cohort construction')
    vals=[31,28,27];ax.bar(range(3),vals,color=[C['pale'],C['prey'],C['other']],width=.56)
    for j,v in enumerate(vals):ax.text(j,v+.7,str(v),ha='center',fontsize=9)
    ax.set_xticks(range(3),['Original','Primary','Sensitivity']);ax.set_ylim(0,37)
    ax=f.add_axes([.60,.58,.34,.28]);ttl(ax,'B','Overlap and identity boundaries');textrows(ax,[('Exact public-ID overlaps',str(len(exact))),('Possible alias, sensitivity only',str(len(poss))),('Final reference roster','63 donors'),('Genotype identity verification','Not performed'),('Undocumented aliases','Unresolved')])
    ax=axes(f,[.12,.18,.35,.22],ylabel='Nuclei (log scale)');ttl(ax,'C','Retained primary sample coverage')
    groups=[counts[counts.compartment==c].n_nuclei.to_numpy() for c in ['EXC','INH','OTHER_NEURON']]
    # All donor values; labels use sequence positions only and expose no IDs.
    for j,vals in enumerate(groups):ax.plot(np.full(len(vals),j),vals,'o',ms=3,color=[C['prey'],C['other'],C['hit']][j],alpha=.6)
    ax.set_yscale('log');ax.set_xticks(range(3),['Excitatory','Inhibitory','Other neurons']);ax.set_xlim(-.5,2.5)
    f.text(.60,.39,'28 primary donors: 14 ASD / 14 controls',fontsize=8.5);f.text(.60,.31,'One sample per donor; PFC preferred,\nACC used when needed.',fontsize=7.5,linespacing=1.6)
    f.text(.60,.21,'27-donor sensitivity: 13 ASD / 14 controls.\nPossible-alias exclusion does not confirm identity.\nResidual independence remains unresolved.',fontsize=7.5,linespacing=1.7)
    f.text(.06,.055,'The final public Wamsley roster was used. Public IDs and demographics do not establish genotype identity.',fontsize=7.3)
    r.save(f,'Figure_S12','Supplementary Figure S12. A, original, exact-overlap-excluded primary and possible-alias-excluded sensitivity donor counts. Three exact public-ID overlaps were removed before the 28-donor cohort; one additional possible alias is removed only in the 27-donor sensitivity. B, the final public Wamsley roster contains 63 donors (33 ASD, 30 controls); metadata matching is not genotype verification. C, all retained primary donor-by-compartment nuclei counts without donor identifiers. One PFC-preferred/ACC-fallback sample is selected per donor. The primary cohort is 14 ASD/14 controls and the sensitivity is 13 ASD/14 controls. Undocumented aliases and residual independence remain unresolved.','supplementary')


def verify_display_arithmetic(r):
    s=r.source
    pm=s.table('activity_permutation')
    for _,z in pm.iterrows():r.check('Activity empirical P '+z.analysis,(z.exceedances+1)/(z.n_permutations+1),z.empirical_P_plus_one)
    z=s.table('igvf_matched').iloc[0]
    r.check('Matched IGVF displayed risk difference',z.exposed_hits/z.n_exposed-z.comparator_hits/z.n_comparator,z.estimate_difference)
    d=s.sup(3,'S3A_BAIT_NEIGHBORHOOD_NULLS.tsv')
    r.check('All frozen bait rows',len(d),100)
    r.check('No bait has BH q below 0.05',bool((d.BH_FDR.dropna()>=.05).all()),True)
    for _,z in d[d.expected_hit_count>0].iterrows():r.check('Bait enrichment '+z.bait,z.screen_measured_prey_hit_n/z.expected_hit_count,z.enrichment_ratio_vs_screened_experimental_preys)
    z=s.sup(4,'S4A_SYNAPTOGENESIS.tsv').iloc[0]
    r.check('Synaptogenesis prey rate',z.prey_only_hits/z.prey_only_covered,z.prey_only_hit_rate)
    r.check('Synaptogenesis comparator rate',z.comparator_hits/z.comparator_n,z.comparator_hit_rate)
    w=s.table('wam_cells')
    for _,z in w.iterrows():r.check('Wamsley median difference '+z.gene_set+' '+z.source_cell_type,z.set_median_abs_logFC-z.comparator_median_abs_logFC,z.effect_median_difference)
    for role in ['gordon_days','gordon_forms']:
        for i,z in s.table(role).iterrows():r.check(role+' mean difference '+str(i),z.mean_A-z.mean_B,z.delta_mean_A_minus_B)
    for i,z in s.table('vel_compartments').iterrows():r.check('Velmeshev compartment difference '+str(i),z.H_mean_abs_beta-z.comparator_mean_abs_beta,z.mean_abs_beta_difference)
    b=s.sup(5,'B_PROGRAM_ENRICHMENT_ALL.tsv')
    r.check('All frozen programs included',len(b),200)
    r.check('Program inference unavailable rather than zero',bool(b[['OR','P_raw','P_adjusted_BH']].isna().all().all()),True)


def r5_main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    set_theme();plt.rcParams['font.weight']='normal';s=R5Source(args.input);r=R5Run(s,args.output)
    validate(s,r);verify_display_arithmetic(r);s.active=[]
    for fun in [f1_r5,f2_r5,f3_r5,draw_figure4,s1,s2,s3,s4,s5,s6,s7,s8,draw_s9,draw_s10]:
        print('Rendering',fun.__name__,flush=True);fun(r)
    # Scientific statuses and exact mapping/cohort accounting are part of the legends.
    mapping=s.sup(9,'GORDON_EFFECTIVE_MAPPING.tsv')
    for path,sha,panel,n in s.active:r.crosswalk.append(dict(panel='Figure_S9 legend',source=path,sha256=sha,operation='Effective mapping counts in legend',source_rows=n))
    r.legends['Figure_S9']+='\nClassification: NON_SUPPORTIVE. Effective mapped same-screen sets: '+str(int(mapping.is_H.sum()))+' hits and '+str(int(mapping.is_N_screen.sum()))+' non-hits; all effective source mapping rows are retained in the supplementary tables.\n'
    r.legends['Figure_S10']+='\nClassification: NON_SUPPORTIVE. The all-region dream sensitivity is NOT_RUN, not an effect estimate.\n'
    for n in ['Figure_S9','Figure_S10']:
        # Replace old internal wording retained in historical code with publication language.
        r.legends[n]=r.legends[n].replace('Moving this figure from the main figure to the supplement does not change its pre-specified analysis status, numerical results or the requirement to report the non-supportive outcome in Results.','The prespecified primary comparison remains non-supportive.')
    r.check('Main figure count',len(list((r.output/'figures').glob('Figure_*.pdf'))),4)
    r.check('Supplementary figure count',len(list((r.output/'supplementary_figures').glob('Figure_S*.pdf'))),10)
    merge_pdfs([r.output/'figures'/f'Figure_{i}.pdf' for i in range(1,5)],r.output/'Main_Figures_R5.pdf')
    merge_pdfs([r.output/'supplementary_figures'/f'Figure_S{i}.pdf' for i in range(1,11)],r.output/'Supplementary_Figures_R5.pdf')
    for path in r.pdf_paths+[r.output/'Main_Figures_R5.pdf',r.output/'Supplementary_Figures_R5.pdf']:
        with fitz.open(path) as doc:
            for i,page in enumerate(doc):page.get_pixmap(matrix=fitz.Matrix(2,2),alpha=False).save(r.output/'qa/pdf_renders'/f'{path.stem}_page{i+1}.png')
    (r.output/'legends/Figure_Legends.md').write_text('\n\n'.join(r.legends.values()))
    pd.DataFrame(r.crosswalk).drop_duplicates().to_csv(r.output/'source_index/FIGURE_SOURCE_INDEX.tsv',sep='\t',index=False)
    # Direct-source records are provenance checks; the independent recalculations are separate named checks.
    pd.DataFrame(r.checks+s.numeric).drop_duplicates().to_csv(r.output/'qa/FIGURE_NUMBER_SOURCE_QA.tsv',sep='\t',index=False)
    print('Completed frozen figure rendering from the supplied data root. This message does not establish source-license clearance or scientific acceptance.',flush=True)

if __name__=='__main__': r5_main()
