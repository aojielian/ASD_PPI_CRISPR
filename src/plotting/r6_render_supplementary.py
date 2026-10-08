#!/usr/bin/env python3
"""R6 supplementary page composition using unchanged frozen R5 plotting data.

No upstream analysis is run. The inherited functions read their original TSVs.
Only wording, type, axes placement and PDF page composition are changed.
"""
from pathlib import Path
import argparse, io, html, json, re, hashlib
import importlib.util
import fitz

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('r5_base',HERE/'r5_render_figures.py')
b=importlib.util.module_from_spec(spec); spec.loader.exec_module(b)
TITLES={
1:'Activity-screen robustness and ascertainment',
2:'Patient-variant rewiring comparisons',
3:'Bait-specific neuronal-activity neighborhoods',
4:'Independent synaptogenesis-screen comparison',
5:'IGVF model sensitivities, matching and program coverage',
6:'In-vivo condition profiles and sensitivity analyses',
7:'Neuronal-context support and comparison coverage',
8:'Wamsley cell-type magnitudes and secondary comparisons',
9:'Gordon cortical-organoid validation',
10:'Velmeshev patient-cortex validation and sensitivities',
11:'External-validation descriptive profiles',
12:'Donor overlap, cohort construction and identity limits'}
CAPTIONS={
1:'''(A) Primary (FDR < 0.10) and strict (FDR < 0.05) activity-hit definitions, with crude, selection-adjusted and Cochran-Mantel-Haenszel (CMH) odds ratios and their original 95% confidence intervals. (B) All primary and strict selection strata; circles denote prey-only genes and squares denote comparators. (C) All eight leave-one-stratum-out CMH estimates and intervals. (D) Summary results from 100,000 within-selection-stratum assignments per definition, plus ASD/NDD-selected subset counts. The zero-hit ascertainment subset is too sparse for stable directional interpretation. Individual permutation draws were not retained in the source package. Targets with incompletely documented selection histories could not be assigned to every stratum.''',
2:'''(A) All variant-pair comparisons: circles show observed counts, squares show expected counts and horizontal lines show the original null 2.5th-97.5th percentiles. Labels give greater-tail empirical P values; these lines are null-reference ranges, not confidence intervals for effects. (B) Unique-gene altered-versus-never-altered sensitivity (2/21 versus 6/34), with the original Fisher odds ratio, 95% confidence interval and P value. (C) Continuous gene-level Spearman association (55 genes), with the original bootstrap interval and nominal P value. These comparisons did not provide evidence that variant-altered interactions identify more neuronal-activity hits.''',
3:'''All 100 bait neighborhoods are shown in alphabetical order. Points show the original hit-enrichment ratios relative to screen-measured experimental prey genes. A dash indicates that a ratio could not be estimated; gray points mark neighborhoods below the prespecified minimum-degree rule. No neighborhood met the original BH q < 0.05 threshold. The supplementary tables contain bait-specific tested counts, hit counts, expected counts, P values, q values and gene memberships.''',
4:'''(A) Source-defined synaptogenesis-hit fractions: 7/91 prey-only genes and 63/544 comparators. (B) The original Fisher odds ratio, 95% confidence interval and two-sided P value. The screen comprised 644 genes and 71 primary-score hits, including assay controls; nine bait or other-prey-role genes were excluded from this comparison. There was no evidence of additional hit enrichment among prey-only genes. Assay coverage and the source-defined tested universe should be considered when interpreting this comparison.''',
5:'''(A) All four logistic models for the same 623 complete cases (120 prey, 503 comparators; 136 events), drawn from 1,186 eligible genes; intervals are penalized profile 95% intervals. (B) Continuous energy-distance sensitivity, with the original HC3 interval and nominal P value. (C) Covariate balance before (squares) and after (circles) matching. (D) The 106-pair matched comparison and original paired-permutation P value. (E) Descriptive cumulative distributions of the number of source-reported significant programs across eligible targets. Tested denominators were not provided for any of the 200 programs, so program-level enrichment could not be estimated. Absence from the significant-only results does not establish a non-hit.''',
6:'''(A) All eight condition-specific mean differences and original two-sided permutation P values; the source did not provide condition-specific confidence intervals. (B, C) Every leave-one-condition-out and leave-one-cell-class-out composite. (D) Age summaries and both alternative gene aggregations, with their original bootstrap 95% intervals. (E) All three prespecified Vglut2 ages, including the secondary 18-month extension; labels give the original three-test Holm-adjusted P values. Omission estimates are sensitivity summaries from overlapping data and do not constitute independent cohorts.''',
7:'''(A) Source-defined numbers of supported and testable prey genes for each neuronal-context dataset; limited testability does not indicate a biological null. (B, C) Unadjusted and source-adjusted support gradients, with original model-specific intervals and P values. Activity effects are odds ratios; in-vivo effects are linear coefficients. (D) Same-bait contrasts use different outcome scales and cannot be compared as standardized magnitudes. (E) Eligible ASD, DD-only and SCZ-only bait counts. Low coverage and inadequate matching prevented formal disease-specificity contrasts. The in-vivo results do not establish specificity relative to generic interaction networks, and cross-study activity comparisons do not establish same-assay specificity.''',
8:'''(A) All 36 cell types, shown as median absolute-logFC differences for PPI prey and activity-hit sets. (B) Non-neuronal omnibus median differences and original one-sided Holm P values. (C) Direct neuronal-minus-nonneuronal class contrasts and original two-sided Holm P values; the source did not provide intervals. (D) Generic-network matched sensitivity with the original bootstrap 95% interval and permutation P value. The subset-versus-rest matching comparison could not be estimated under its prespecified balance criterion. (E) Primary and network-augmented same-screen models, with their own rank-normal coefficient intervals and raw P values. The same-screen comparison did not support additional prioritization of activity hits. Separate within-class associations do not establish neuronal specificity.''',
9:'''(A) All 21 activity hits and 97 same-screen PPI non-hits. Each gene's score is its mean absolute author logFC across nine ASD forms and four developmental days (C36). Bars show the median, interquartile range and Tukey whiskers; all gene points are retained. The original two-sided selection-stratified competitive P value used 1,000,000 assignments. This gene-competitive comparison did not support an advantage for the activity-hit set. (B) All four days under C36 and the prespecified six-form C24 sensitivity. (C) All nine C36 forms in alphabetical order; gray squares identify low-donor forms (PCDH19, SHANK3 and TS). Day and form profiles are descriptive, with no new P values or intervals. (D) The primary, C24 and unstratified estimates and original P values. The source supplied no confidence intervals for these estimates. C24 is a sensitivity analysis; it does not replace the C36 primary endpoint.''',
10:'''The primary cohort contains 28 donors (14 ASD, 14 controls), after excluding three exact public-ID overlaps with Wamsley. One PFC-preferred or ACC sample was selected per donor. The 27-donor sensitivity additionally excludes one possible alias (13 ASD, 14 controls). (A) All expression-qualified hits (19), same-screen non-hits (89) and other PPI genes (888), with original mean absolute coefficients across three neuronal compartments. Bars show medians, interquartile ranges and Tukey whiskers, not confidence intervals. (B) Original ROMER Mixed rank-enrichment estimates and Holm P values from 49,999 rotations; the two comparisons form a correction family within each cohort. Neither primary comparison supported activity-hit prioritization. Rank differences use different comparator universes; no confidence intervals were supplied. (C) Descriptive absolute-coefficient differences in all three compartments. The voom sensitivity is retained in the supplementary tables. The all-region sensitivity analysis was not performed and is not represented as a zero effect. (D) Original, uniquely mapped and expression-qualified gene counts; bar heights are retained percentages. Public IDs and demographic matching do not establish genotype identity (S12).''',
11:'''(A) Gordon group means across all developmental days under C36 and the prespecified C24 sensitivity. (B) Group means for all nine C36 forms. (C) Velmeshev source-estimated mean absolute coefficients in each of the three neuronal compartments for activity hits, same-screen non-hits and other PPI genes. These profiles are descriptive: no day-, form- or compartment-specific P values or confidence intervals were added. The primary Gordon and Velmeshev comparisons did not support the proposed advantage for the activity-hit set.''',
12:'''(A) Donor counts before overlap exclusion, in the primary cohort and in the possible-alias sensitivity. Three exact public-ID overlaps were removed before the 28-donor primary cohort; one additional possible alias is removed only in the 27-donor sensitivity. (B) Public-roster overlap and identity information. The final Wamsley roster contains 63 donors (33 ASD, 30 controls). Public IDs and demographics were used for matching; genotype identity was not assessed. (C) Every retained primary donor-by-compartment nuclei count, with donor identifiers omitted. One PFC-preferred or ACC sample was selected per donor. The primary cohort contains 14 ASD and 14 controls; the sensitivity contains 13 ASD and 14 controls. Possible aliases do not establish identical donors, and undocumented aliases could affect independence.'''}

# Explicit display-only wording substitutions, never replacements in source tables.
REPLACE={
'All frozen variant-pair permutation comparisons':'All variant-pair comparisons',
'Frozen comparison':'Synaptogenesis association',
'Disease-specificity feasibility':'Disease-specificity coverage',
'Overlap and identity boundaries':'Public-roster identity information',
'Unresolved':'Not established',
'Undocumented aliases':'Additional alias identities',
'Original 95% interval':'95% interval',
'Original maximum |SMD|':'Maximum |SMD|',
}
REMOVE_PREFIXES=(
'Sparse ascertainment subset:', 'Selection strata remain incomplete',
'All 100 baits in alphabetical', 'Dash: enrichment unavailable',
'Non-supportive external comparison;',
'All 200 programs:', 'Tested denominators unavailable.',
'Disease-specificity tests:', 'Sparse eligible bait sets', 'Generic AP–MS specificity',
'Same-screen prioritization remains',
'Descriptive profiles only','No day-, form- or compartment-level',
'Gordon and Velmeshev primary classifications',
'27-donor sensitivity:','The final public Wamsley roster was used.',
)

def plot_payload(f):
    items=[]
    for ai,ax in enumerate(f.axes):
        for ln in ax.lines:items.append({'axis':ai,'kind':'line','x':b.np.asarray(ln.get_xdata()).tolist(),'y':b.np.asarray(ln.get_ydata()).tolist()})
        for c in ax.collections:
            if hasattr(c,'get_offsets'):items.append({'axis':ai,'kind':'offsets','xy':b.np.asarray(c.get_offsets()).tolist()})
        for p in ax.patches:
            if isinstance(p,b.matplotlib.patches.Rectangle):items.append({'axis':ai,'kind':'bar','x':p.get_x(),'y':p.get_y(),'w':p.get_width(),'h':p.get_height()})
    # JSON equality handles NaN consistently and gives an exact unchanged-data check.
    return json.dumps(items,sort_keys=True,default=str)

class R6Run(b.R5Run):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs); self.pages=[];self.records=[]
    def save(self,f,name,legend,category='supplementary'):
        n=int(name.split('S')[-1]); title=TITLES[n]
        before_payload = plot_payload(f)
        # Remove the old figure header; all pages receive one identical header.
        for t in list(f.texts):
            st=t.get_text()
            if re.match(r'^S\d+\s',st) or st.startswith(REMOVE_PREFIXES):t.remove()
        # Body caveats are conveyed in the bottom caption; numeric annotations remain.
        for t in f.findobj(match=b.Text):
            st=t.get_text()
            for old,new in REPLACE.items():st=st.replace(old,new)
            st=st.replace('unresolved','unspecified')
            t.set_text(st)
            if t.get_text().startswith('Observed altered-pair count:'):t.remove()
        # Use one 190 mm-wide body canvas, constant physical type, and panel area.
        f.set_size_inches(190/25.4,210/25.4,forward=True)
        for t in f.findobj(match=b.Text):
            if t.get_visible() and t.get_text():
                t.set_color(b.C['ink'])
                if t.get_fontsize()<7:t.set_fontsize(7)
        # Supplement-specific spacing adjustments; data and axis limits are retained.
        if n==1:
            for ax in f.axes:
                labels=[t.get_text().replace('unresolved','unspecified') for t in ax.get_yticklabels()]
                if labels:ax.set_yticks(ax.get_yticks(),labels)
        if n==3:
            for ax in f.axes:ax.tick_params(axis='y',labelsize=7)
        if n==5:
            f.axes[-1].set_position([.12,.10,.80,.13])
        if n==8:
            f.axes[0].set_position([.18,.36,.30,.54]);f.axes[0].tick_params(axis='y',labelsize=6.9)
            f.axes[1].set_title('B  Non-neuronal',loc='left',pad=12,fontsize=9,fontweight='bold')
            f.axes[2].set_title('C  Class contrast',loc='left',pad=12,fontsize=9,fontweight='bold')
            f.axes[3].set_title('D  Matched network control',loc='left',pad=12,fontsize=9,fontweight='bold')
            f.axes[4].set_title('E  Same-screen models',loc='left',pad=12,fontsize=9,fontweight='bold')
            f.axes[1].set_xlabel('Median |logFC| difference',fontsize=7)
            f.axes[2].set_xlabel('Difference in median contrasts',fontsize=7)
            f.axes[4].set_xlabel('Rank-normal coefficient',fontsize=7)
            for t in f.texts:
                if t.get_text().startswith('Subset-versus-rest matching:'):
                    t.set_text(t.get_text().replace('Subset-versus-rest matching: not estimable','Subset matching not estimated'))
        if n==12:
            f.axes[1].set_title('B  Public-roster identity information',loc='left',pad=12,fontsize=9,fontweight='bold')
            for t in f.axes[1].texts:t.set_fontsize(7.1)
            for t in f.texts:
                if t.get_text().startswith('One sample per donor;'):t.remove()
        f.canvas.draw()
        # Preserve all artist-data payloads for source comparison by the parent.
        coords=plot_payload(f)
        if coords != before_payload:raise RuntimeError(f'Data payload changed: {name}')
        raw=io.BytesIO();f.savefig(raw,format='pdf',facecolor='white',metadata={'Title':name,'Author':'','CreationDate':None,'ModDate':None})
        src=fitz.open(stream=raw.getvalue(),filetype='pdf');doc=fitz.open();page=doc.new_page(width=595.2756,height=841.8898)
        # A4: 10 mm margins; 210 mm panel body; 54 mm fixed caption region.
        page.insert_text((28.346,24),f'Figure S{n}  |  {title}',fontname='hebo',fontsize=11,color=(.14,.16,.18))
        page.show_pdf_page(fitz.Rect(28.346,38,566.929,633.276),src,0,keep_proportion=True)
        page.draw_line((28.346,645),(566.929,645),color=(.72,.75,.76),width=.5)
        caption=f'<b>Figure S{n}. {html.escape(title)}.</b> '+html.escape(CAPTIONS[n])
        css='* {font-family: sans-serif; font-size: 8.5pt; line-height: 1.20; color:#24282d;} p {text-align:justify; margin:0;}'
        spare,scale=page.insert_htmlbox(fitz.Rect(28.346,655,566.929,817),'<p>'+caption+'</p>',css=css,scale_low=1)
        if spare<0:raise RuntimeError(f'Caption does not fit at fixed size: {name}')
        base=self.output/'supplementary_figures'/name
        doc.set_metadata({'title':f'Figure S{n}. {title}','author':'','creator':'Scientific figure renderer'})
        doc.save(base.with_suffix('.pdf'),garbage=4,deflate=True)
        page.get_pixmap(dpi=self.dpi,alpha=False).save(base.with_suffix('.png'))
        self.pdf_paths.append(base.with_suffix('.pdf'));self.legends[name]=f'Figure S{n}. {title}.\n\n'+CAPTIONS[n]
        self.records.append({'figure':name,'page_points':[page.rect.width,page.rect.height],'caption_font_points':8.5,'caption_scale':scale,'caption_spare_points':spare,'pdf_sha256':hashlib.sha256(base.with_suffix('.pdf').read_bytes()).hexdigest(),'unchanged_plot_payload':True,'plot_payload_sha256':hashlib.sha256(json.dumps(coords,sort_keys=True,default=str).encode()).hexdigest()})
        for path,sha,panel,rows in self.source.active:self.crosswalk.append({'panel':name,'source':path,'sha256':sha,'source_rows':rows})
        self.source.active=[];b.plt.close(f);doc.close();src.close()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--dpi',type=int,default=200);a=ap.parse_args()
    b.set_theme();s=b.R5Source(a.input);r=R6Run(s,a.output,a.dpi)
    for fn in [b.s1,b.s2,b.s3,b.s4,b.s5,b.s6,b.s7,b.s8,b.draw_s9,b.draw_s10]:
        print('Rendering',fn.__name__,flush=True);fn(r)
    b.merge_pdfs(r.pdf_paths,a.output/'Supplementary_Figures_R6.pdf')
    (a.output/'legends'/'Supplementary_Figure_Legends.md').write_text('\n\n'.join(r.legends.values())+'\n',encoding='utf-8')
    # Structured execution evidence on stdout, so the parent retains ownership of QA.
    print(json.dumps({'page_records':r.records,'source_crosswalk':r.crosswalk,'source_numeric_checks':s.numeric},indent=2),flush=True)
if __name__=='__main__':main()
