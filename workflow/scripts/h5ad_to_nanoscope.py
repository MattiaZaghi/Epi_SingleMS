#!/usr/bin/env python
"""
Build nanoscope-style `cell_picking/metadata.csv` files from ONE published h5ad
(the paper's final object, all patients merged) plus ONE cellranger fragments
file PER PATIENT.

    inspect    look at the h5ad and a fragments file, work out the barcode
               format and which .obs column identifies the patient
    calibrate  decide how to scale all_unique_MB so it is comparable to your
               existing nanoscope samples
    build      write one metadata.csv per patient

WHY IT IS SHAPED THIS WAY
-------------------------
Raw 10x barcodes are only unique *within* a library, so the same 16-mer appears
in every patient's fragments file. A merged published object therefore
disambiguates cells somehow -- usually a prefix/suffix on the barcode
("P1_AAAC...-1", "AAAC...-1_3") and/or an .obs column. Matching must be on the
pair (patient, 16bp barcode); matching on barcode alone silently mixes patients
together. `inspect` exists to find that out before you commit to anything.

CELL CALLING
------------
The h5ad is the paper's final, already-filtered object, so its membership *is*
the published cell call. By default `passedMB` = "this barcode is in the h5ad
for this patient", which reproduces the paper's cell set -- what you want for a
comparison dataset. Metrics are still computed for every barcode in the
fragments file, so the full droplet cloud is available and plotPassedCells()
looks like it does for your own samples.

`--cell-call gmm` instead re-runs nanoscope's own criteria. Use it only if you
deliberately want your gate rather than the paper's, and note it will not
reproduce the paper's numbers.

WHAT CANNOT BE RECONSTRUCTED
----------------------------
metadata.csv normally merges cellranger's singlecell.csv (alignment QC) with
per-barcode read counts from possorted_bam.bam. From a fragments file you get
passed_filters, peak_region_fragments, mitochondrial, total, duplicate,
all_unique_MB, peak_MB and every derived ratio. You cannot get chimeric,
unmapped, lowmapq, nonprimary, peak_region_cutsites, or the
DNase/enhancer/on_target/blacklist counts; those are written NA rather than
invented. Nothing in plotPassedCells() or QC*.R reads them.

all_unique_MB IS NOT A FRAGMENT COUNT -- see `calibrate --help`.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

BARCODE_RE_DEFAULT = r"([ACGT]{16})"

CELLRANGER_COLS = [
    "total", "duplicate", "chimeric", "unmapped", "lowmapq", "mitochondrial",
    "nonprimary", "passed_filters", "is__cell_barcode", "excluded_reason",
    "TSS_fragments", "DNase_sensitive_region_fragments",
    "enhancer_region_fragments", "promoter_region_fragments",
    "on_target_fragments", "blacklist_region_fragments",
    "peak_region_fragments", "peak_region_cutsites",
]
DERIVED_COLS = [
    "logUMI", "promoter_ratio", "peak_region_ratio",
    "all_unique_MB", "peak_MB", "peak_ratio_MB", "sample",
]


def need(tool):
    if shutil.which(tool) is None:
        sys.exit(f"ERROR: '{tool}' not on PATH (conda env nanoscope_general).")


def extract_bc(series, regex):
    """Pull the 16bp barcode out of whatever decoration surrounds it."""
    ex = series.astype(str).str.extract(regex, expand=False)
    return ex


def count_per_barcode(fragments, bed=None, chrom=None, tmpdir=None):
    """Per-barcode fragment count and summed read-pair support (column 5)."""
    stream = f"gzip -dc {fragments} | grep -v '^#'"
    if chrom:
        stream += f" | awk -F'\\t' '$1==\"{chrom}\"'"
    if bed:
        need("bedtools")
        stream += f" | bedtools intersect -a stdin -b {bed} -u"
    awk = ("awk -F'\\t' '{n[$4]++; r[$4]+=($5==\"\"?1:$5)} "
           "END{for (b in n) printf \"%s\\t%d\\t%d\\n\", b, n[b], r[b]}'")
    cmd = f"set -o pipefail; {stream} | {awk}"
    with tempfile.NamedTemporaryFile("w+", suffix=".tsv", dir=tmpdir, delete=False) as fh:
        out = fh.name
    try:
        with open(out, "w") as fh:
            subprocess.run(["bash", "-c", cmd], stdout=fh, check=True)
        return pd.read_csv(out, sep="\t", header=None,
                           names=["raw_barcode", "n_frag", "n_reads"])
    finally:
        os.unlink(out)


def frag_barcodes_head(fragments, n=5):
    cmd = (f"set -o pipefail; gzip -dc {fragments} | grep -v '^#' | "
           f"cut -f4 | head -{n * 200} | sort -u | head -{n}")
    r = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    return [x for x in r.stdout.strip().split("\n") if x]


# ----------------------------------------------------------------- inspect ----

def cmd_inspect(a):
    try:
        import anndata as ad
    except ImportError:
        sys.exit("ERROR: pip install anndata")

    print(f"=== {a.h5ad} ===")
    adata = ad.read_h5ad(a.h5ad, backed="r")
    print(f"{adata.n_obs} cells x {adata.n_vars} features\n")

    idx = pd.Series(adata.obs_names.astype(str))
    print("obs_names (cell ids):")
    for v in idx.head(5):
        print("   ", v)
    bc = extract_bc(idx, a.barcode_regex)
    print(f"\n16bp barcode extracted with {a.barcode_regex!r}: "
          f"{bc.notna().sum()}/{len(bc)} matched")
    if bc.notna().any():
        print("    e.g.", list(bc.dropna().head(3)))
        print(f"    unique barcodes: {bc.nunique()}   "
              f"(if < n_cells, barcodes repeat across patients -> "
              f"you MUST match on patient+barcode)")
    else:
        print("    !! no match. Try --barcode-regex, e.g. r'([ACGT]{16})' or "
              r"r'([ACGTN]+)'")

    print("\n.obs columns (candidates for --h5ad-sample-col marked *):")
    for c in adata.obs.columns:
        s = adata.obs[c]
        nu = s.nunique(dropna=True)
        star = "*" if (1 < nu <= 64 and (s.dtype == object or
                                         str(s.dtype) == "category")) else " "
        vals = list(pd.Series(s.unique()).dropna().astype(str)[:6])
        print(f"  {star} {c:<32} {str(s.dtype):<10} nuniq={nu:<6} {vals}")

    if a.fragments:
        print(f"\n=== {a.fragments} ===")
        fb = frag_barcodes_head(a.fragments)
        print("fragments column 4:")
        for v in fb:
            print("   ", v)
        fbc = extract_bc(pd.Series(fb), a.barcode_regex)
        print(f"extracted: {list(fbc.dropna())}")
        if bc.notna().any() and fbc.notna().any():
            shared = set(bc.dropna()) & set(fbc.dropna())
            print(f"\n{len(shared)} of these {len(fbc.dropna())} example barcodes "
                  f"appear in the h5ad -> extraction is consistent"
                  if shared else
                  "\nNone of the example barcodes appear in the h5ad. Either this "
                  "patient's cells\nare a small subset, or the regex is wrong.")

    if getattr(a, "combo", None):
        cols = [c.strip() for c in a.combo.split(",") if c.strip()]
        miss = [c for c in cols if c not in adata.obs.columns]
        if miss:
            print(f"\n!! no such .obs column(s): {miss}")
        else:
            key = adata.obs[cols].astype(str).agg(a.h5ad_sample_sep.join, axis=1)
            vc = key.value_counts()
            print(f"\n{len(vc)} distinct {' + '.join(cols)} combinations "
                  f"(these are your samplesheet h5ad_value entries):")
            for k, n in vc.items():
                print(f"    {k:<48} {n:>7} cells")

    print("\nNext: write a samplesheet (sample,h5ad_value,fragments[,peaks]) and run "
          "`build`.\n"
          "      Use the starred .obs column above for --h5ad-sample-col, and put "
          "its\n      exact values in the h5ad_value column.")


# --------------------------------------------------------------- calibrate ----

def cmd_calibrate(a):
    """all_unique_MB is counted by nanoscope as
           samtools view -f2 possorted_bam.bam | count CB
       -f2 keeps properly-paired ALIGNMENTS and does not drop duplicates, so it
       is ~2x reads per barcode, not fragments per barcode. Counting a fragments
       file directly shifts logUMI down by ~0.3-0.6 and makes the GMM and hard
       cutoffs non-comparable with your existing samples. This picks the scale
       empirically against a sample you already processed."""
    cnt = count_per_barcode(a.fragments, tmpdir=a.tmpdir)
    cnt["barcode"] = extract_bc(cnt.raw_barcode, a.barcode_regex)
    cnt = cnt.dropna(subset=["barcode"]).set_index("barcode")

    ref = pd.read_csv(a.nanoscope_metadata)
    if "all_unique_MB" not in ref.columns:
        sys.exit("ERROR: reference metadata has no all_unique_MB column.")
    ref["barcode"] = extract_bc(ref["barcode"], a.barcode_regex)
    ref = ref.dropna(subset=["barcode"]).set_index("barcode")

    common = cnt.index.intersection(ref.index)
    if len(common) == 0:
        sys.exit("ERROR: no shared barcodes -- check --barcode-regex.")
    print(f"\n{len(common)} shared barcodes "
          f"({len(cnt)} in fragments, {len(ref)} in reference)\n")

    truth = ref.loc[common, "all_unique_MB"].astype(float)
    opts = {"fragments": cnt.loc[common, "n_frag"].astype(float),
            "reads":     cnt.loc[common, "n_reads"].astype(float) * 2,
            "x2":        cnt.loc[common, "n_frag"].astype(float) * 2}
    print(f"{'option':>10}  {'median ratio':>13}  {'log10 offset':>13}  {'spearman':>9}")
    print(f"{'reference':>10}  {1.0:13.3f}  {0.0:13.3f}  {1.0:9.3f}"
          f"   (median = {truth.median():.0f})")
    for nm, v in opts.items():
        ratio = (v / truth).replace([np.inf, -np.inf], np.nan).dropna()
        off = (np.log10(v.replace(0, np.nan)) -
               np.log10(truth.replace(0, np.nan))).dropna()
        print(f"{nm:>10}  {ratio.median():13.3f}  {off.median():13.3f}  "
              f"{v.corr(truth, method='spearman'):9.3f}")
    print("\nPass the option with median ratio nearest 1.000 to `build --read-scale`.")


# ------------------------------------------------------------------- build ----

def build_one(sample, h5ad_value, fragments, peaks, obs_sub, a):
    allc = count_per_barcode(fragments, tmpdir=a.tmpdir)
    allc["barcode"] = extract_bc(allc.raw_barcode, a.barcode_regex)
    allc = allc.dropna(subset=["barcode"])
    if allc.barcode.duplicated().any():
        allc = allc.groupby("barcode", as_index=False)[["n_frag", "n_reads"]].sum()
    allc = allc.set_index("barcode")

    pk = count_per_barcode(fragments, bed=peaks, tmpdir=a.tmpdir)
    pk["barcode"] = extract_bc(pk.raw_barcode, a.barcode_regex)
    pk = pk.dropna(subset=["barcode"]).groupby("barcode")["n_frag"].sum()

    mito = count_per_barcode(fragments, chrom=a.mito_chrom, tmpdir=a.tmpdir)
    mito["barcode"] = extract_bc(mito.raw_barcode, a.barcode_regex)
    mito = mito.dropna(subset=["barcode"]).groupby("barcode")["n_frag"].sum()

    prom = None
    if a.promoters:
        p = count_per_barcode(fragments, bed=a.promoters, tmpdir=a.tmpdir)
        p["barcode"] = extract_bc(p.raw_barcode, a.barcode_regex)
        prom = p.dropna(subset=["barcode"]).groupby("barcode")["n_frag"].sum()

    md = pd.DataFrame(index=allc.index)
    md.index.name = "barcode"
    for c in CELLRANGER_COLS:
        md[c] = np.nan

    n_frag = allc["n_frag"].astype(float)
    n_read = allc["n_reads"].astype(float)
    md["passed_filters"] = n_frag
    md["total"] = n_read
    md["duplicate"] = n_read - n_frag
    md["mitochondrial"] = mito.reindex(md.index).fillna(0)
    md["peak_region_fragments"] = pk.reindex(md.index).fillna(0)
    md["excluded_reason"] = 0
    if prom is not None:
        pv = prom.reindex(md.index).fillna(0)
        md["promoter_region_fragments"] = pv
        md["TSS_fragments"] = pv

    md["logUMI"] = np.log10(md["passed_filters"] + 1)
    md["promoter_ratio"] = (md["promoter_region_fragments"] + 1) / (md["passed_filters"] + 1)
    md["peak_region_ratio"] = (md["peak_region_fragments"] + 1) / (md["passed_filters"] + 1)

    scale = {"fragments": n_frag, "reads": n_read * 2, "x2": n_frag * 2}[a.read_scale]
    md["all_unique_MB"] = scale
    md["peak_MB"] = md["peak_region_fragments"] * (scale / n_frag.replace(0, np.nan))
    md["peak_ratio_MB"] = md["peak_MB"] / md["all_unique_MB"]
    md["sample"] = sample
    if a.modality:
        md["modality"] = a.modality
        md["mark"] = a.modality

    # ---- the paper's cell call ------------------------------------------------
    cells = set(obs_sub["_bc16"].dropna())
    md["is__cell_barcode"] = md.index.isin(cells).astype(int)
    n_match = int(md["is__cell_barcode"].sum())

    if a.cell_call == "h5ad":
        md["passedMB"] = md.index.isin(cells)
    else:
        md["passedMB"] = ((md.all_unique_MB > 10**a.min_reads) &
                          (md.all_unique_MB < 10**a.max_reads) &
                          (md.peak_ratio_MB > a.peak_fraction_min) &
                          (md.peak_ratio_MB < a.peak_fraction_max))
    md["passedMB_legacy"] = ((md.all_unique_MB > 10**a.min_reads) &
                             (md.all_unique_MB < 10**a.max_reads) &
                             (md.peak_ratio_MB > a.peak_fraction_min) &
                             (md.peak_ratio_MB < a.peak_fraction_max))
    md["class"] = np.where(md["passedMB"], 1, 2)

    # carry the paper's annotation onto the matched cells
    ann = obs_sub.dropna(subset=["_bc16"]).drop_duplicates("_bc16").set_index("_bc16")
    ann = ann.drop(columns=["_bc16", "_grp"], errors="ignore").add_prefix("h5ad_")
    md = md.join(ann, how="left")

    md = md.reset_index()
    front = ["barcode"] + CELLRANGER_COLS + DERIVED_COLS + \
            ["passedMB", "passedMB_legacy", "class"]
    md = md[[c for c in front if c in md.columns] +
            [c for c in md.columns if c not in front]]
    return md, n_match, len(obs_sub)


def cmd_build(a):
    try:
        import anndata as ad
    except ImportError:
        sys.exit("ERROR: pip install anndata")

    sheet = pd.read_csv(a.samplesheet)
    for c in ("sample", "fragments"):
        if c not in sheet.columns:
            sys.exit(f"ERROR: samplesheet needs a '{c}' column.")
    if "h5ad_value" not in sheet.columns:
        sheet["h5ad_value"] = sheet["sample"]

    adata = ad.read_h5ad(a.h5ad, backed="r")
    obs = adata.obs.copy()
    obs["_bc16"] = extract_bc(pd.Series(adata.obs_names.astype(str),
                                        index=obs.index), a.barcode_regex).values
    unmatched = obs["_bc16"].isna().sum()
    if unmatched:
        print(f"WARNING: {unmatched}/{len(obs)} h5ad cells gave no 16bp barcode "
              f"with {a.barcode_regex!r}", file=sys.stderr)

    if a.h5ad_sample_col:
        # comma-separated list -> composite key. Needed for designs like this
        # paper's (3 donors x 18 brain regions, hash-multiplexed), where a
        # single .obs column does not identify a library.
        cols = [c.strip() for c in a.h5ad_sample_col.split(",") if c.strip()]
        missing = [c for c in cols if c not in obs.columns]
        if missing:
            sys.exit(f"ERROR: .obs has no column(s) {missing}. "
                     f"Run `inspect` to list them.")
        obs["_grp"] = (obs[cols].astype(str)
                       .agg(a.h5ad_sample_sep.join, axis=1) if len(cols) > 1
                       else obs[cols[0]].astype(str))
        if len(cols) > 1:
            print(f"Composite patient key: {' + '.join(cols)} "
                  f"(joined with {a.h5ad_sample_sep!r})", file=sys.stderr)
    else:
        print("WARNING: no --h5ad-sample-col. Every patient will be matched "
              "against ALL h5ad cells, so barcodes shared between patients will "
              "be assigned to more than one. Only safe if barcodes are globally "
              "unique.", file=sys.stderr)
        obs["_grp"] = None

    rows = []
    for _, r in sheet.iterrows():
        sample = str(r["sample"])
        peaks = r.get("peaks") if isinstance(r.get("peaks"), str) and r.get("peaks") else a.peaks
        if not peaks:
            sys.exit(f"ERROR: no peaks for {sample} (samplesheet column or --peaks).")
        print(f"\n=== {sample} ===", file=sys.stderr)

        if a.h5ad_sample_col:
            sub = obs[obs["_grp"] == str(r["h5ad_value"])]
            if sub.empty:
                print(f"  WARNING: no h5ad cells with {a.h5ad_sample_col}="
                      f"{r['h5ad_value']!r}. Values present: "
                      f"{sorted(obs['_grp'].unique())[:10]}", file=sys.stderr)
        else:
            sub = obs

        md, n_match, n_h5ad = build_one(sample, r["h5ad_value"], r["fragments"],
                                        peaks, sub, a)
        out = os.path.join(a.outdir, sample, a.modality or "modality",
                           "cell_picking", "metadata.csv")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        md.to_csv(out, index=False)

        pct = 100 * n_match / n_h5ad if n_h5ad else float("nan")
        print(f"  {len(md)} barcodes in fragments; {n_h5ad} h5ad cells for this "
              f"patient; {n_match} matched ({pct:.1f}%)", file=sys.stderr)
        print(f"  median all_unique_MB={md.all_unique_MB.median():.0f}  "
              f"median peak_ratio_MB={md.peak_ratio_MB.median():.3f}", file=sys.stderr)
        print(f"  -> {out}", file=sys.stderr)
        if n_h5ad and pct < 80:
            print("  !! low match rate: check --barcode-regex and "
                  "--h5ad-sample-col mapping before using this.", file=sys.stderr)
        rows.append(dict(sample=sample, barcodes=len(md), h5ad_cells=n_h5ad,
                         matched=n_match, pct=round(pct, 1),
                         med_all_unique_MB=round(float(md.all_unique_MB.median()), 1),
                         med_peak_ratio_MB=round(float(md.peak_ratio_MB.median()), 4),
                         out=out))

    summ = pd.DataFrame(rows)
    print("\n" + summ.to_string(index=False), file=sys.stderr)
    tot_m, tot_h = summ.matched.sum(), summ.h5ad_cells.sum()
    print(f"\nmatched {tot_m} cells; h5ad has {adata.n_obs} in total", file=sys.stderr)
    if a.h5ad_sample_col and tot_m < 0.9 * adata.n_obs:
        print("NOTE: fewer cells matched than the h5ad contains -- either some "
              "patients\n      are missing from the samplesheet, or the "
              "patient-value mapping is off.", file=sys.stderr)
    if a.summary:
        summ.to_csv(a.summary, index=False)
        print(f"summary -> {a.summary}", file=sys.stderr)
    print("\nNOT reconstructed (BAM-only, NA): chimeric, unmapped, lowmapq, "
          "nonprimary,\n  peak_region_cutsites, DNase/enhancer/on_target/blacklist.",
          file=sys.stderr)


# -------------------------------------------------------------------- main ----

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("--barcode-regex", default=BARCODE_RE_DEFAULT,
                       help=f"capture group for the 16bp barcode "
                            f"(default {BARCODE_RE_DEFAULT!r})")
        p.add_argument("--tmpdir", default=None)

    p = sub.add_parser("inspect", help="show h5ad obs/barcode format")
    p.add_argument("--h5ad", required=True)
    p.add_argument("--fragments", help="one fragments file, to check barcodes match")
    p.add_argument("--combo", help="comma-separated .obs columns; print the "
                                   "distinct combinations (candidate libraries)")
    p.add_argument("--h5ad-sample-sep", default="|")
    common(p)

    p = sub.add_parser("calibrate", help="choose --read-scale against a known sample")
    p.add_argument("--fragments", required=True)
    p.add_argument("--nanoscope-metadata", required=True)
    common(p)

    p = sub.add_parser("build", help="write one metadata.csv per patient")
    p.add_argument("--h5ad", required=True)
    p.add_argument("--samplesheet", required=True,
                   help="CSV: sample,h5ad_value,fragments[,peaks]")
    p.add_argument("--h5ad-sample-col",
                   help=".obs column identifying the library/patient. Accepts a "
                        "comma-separated list for a composite key, e.g. "
                        "'donor,region' (see `inspect`).")
    p.add_argument("--h5ad-sample-sep", default="|",
                   help="separator joining a composite --h5ad-sample-col "
                        "(default '|'); h5ad_value in the samplesheet must use it")
    p.add_argument("--peaks", help="peak BED used for all samples unless "
                                   "overridden per row")
    p.add_argument("--promoters")
    p.add_argument("--modality", default="")
    p.add_argument("--outdir", default=".")
    p.add_argument("--summary", help="write a per-sample summary CSV here")
    p.add_argument("--read-scale", choices=["fragments", "reads", "x2"],
                   default="reads")
    p.add_argument("--cell-call", choices=["h5ad", "gmm"], default="h5ad",
                   help="h5ad (default): passedMB = present in the published "
                        "object. gmm: nanoscope's own cutoffs instead.")
    p.add_argument("--min-reads", type=float, default=3.0)
    p.add_argument("--max-reads", type=float, default=5.5)
    p.add_argument("--peak-fraction-min", type=float, default=0.2)
    p.add_argument("--peak-fraction-max", type=float, default=1.0)
    p.add_argument("--mito-chrom", default="chrM")
    common(p)

    a = ap.parse_args()
    {"inspect": cmd_inspect, "calibrate": cmd_calibrate, "build": cmd_build}[a.cmd](a)


if __name__ == "__main__":
    main()