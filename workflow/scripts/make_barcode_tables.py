#!/usr/bin/env python3
#!/usr/bin/env python3
"""Generate per-experiment, per-modality sinto barcode tables from a merged
cell-metadata Excel file that has one sheet per modality.

Each output <experiment>_<modality>.tsv is what `sinto filterbarcodes -c`
expects: tab-separated, no header, column 1 = cell barcode (matches the BAM
CB tag), column 2 = filesystem-safe cluster label.

Splitting is done PER (EXPERIMENT, MODALITY) because the same 10x barcode
recurs both across experiments and across modalities in a merged multimodal
object (each modality has its own BAM, e.g. {sample}/{modality}_{barcode}/
cellranger/outs/possorted_bam.bam) -- only (experiment, modality, barcode)
is unique.
"""
import argparse, os, re
import pandas as pd
import yaml


def sanitize(x):
    """Filesystem-safe cluster name: 'OPCs + COPs' -> 'OPCs_COPs'."""
    return re.sub(r"[^0-9A-Za-z]+", "_", str(x)).strip("_")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", default="MS_metadata.xlsx")
    ap.add_argument("--cluster-col", default="predicted.id_Human_MS")
    ap.add_argument("--barcode-col", default="barcode")
    ap.add_argument("--experiment-col", default="experiment")
    ap.add_argument("--filter-col", default="passedMB",
                    help="keep rows where this column is True; pass '' to skip")
    ap.add_argument("--modality-col", default="modality",
                    help="column name to store/read the per-sheet modality label")
    ap.add_argument("--outdir", default="cluster_tracks")
    a = ap.parse_args()

    if a.meta.lower().endswith((".xlsx", ".xls")):
        # sheet_name=None -> dict of {sheet_name: DataFrame} for ALL sheets.
        # A bare pd.read_excel(path) only reads the first sheet, which is why
        # a workbook with one sheet per modality was previously collapsing
        # down to just whichever modality happened to be listed first.
        sheets = pd.read_excel(a.meta, sheet_name=None)
        parts = []
        for sheet_name, sheet_df in sheets.items():
            sheet_df = sheet_df.copy()
            sheet_df[a.modality_col] = sheet_name
            parts.append(sheet_df)
        df = pd.concat(parts, ignore_index=True)
    else:
        df = pd.read_csv(a.meta)
        if a.modality_col not in df.columns:
            raise SystemExit(
                f"--meta is not an Excel file, so it can't be split by sheet; "
                f"expected a '{a.modality_col}' column to identify modality per row."
            )

    if a.filter_col and a.filter_col in df.columns:
        df = df[df[a.filter_col] == True]                      # noqa: E712
    df = df[df[a.cluster_col].notna()].copy()
    df["cluster_safe"] = df[a.cluster_col].map(sanitize)

    bcdir = os.path.join(a.outdir, "barcodes")
    os.makedirs(bcdir, exist_ok=True)

    (df[[a.cluster_col, "cluster_safe"]].drop_duplicates()
       .sort_values("cluster_safe")
       .to_csv(os.path.join(a.outdir, "cluster_name_map.tsv"), sep="\t", index=False))

    bams = {}
    for (exp, modality), sub in df.groupby([a.experiment_col, a.modality_col]):
        key = f"{exp}_{modality}"
        sub[[a.barcode_col, "cluster_safe"]].to_csv(
            os.path.join(bcdir, f"{key}.tsv"), sep="\t", header=False, index=False)
        # Matches the repo convention {sample}/{modality}_{barcode}/cellranger/outs/
        # possorted_bam.bam -- the antibody barcode suffix isn't in the metadata,
        # so it's still left for the user to fill in.
        bams[key] = f"FILL_ME/{exp}/{modality}_FILL_ME/cellranger/outs/possorted_bam.bam"
        print(f"{key:44s} {len(sub):6d} cells  {sub['cluster_safe'].nunique()} clusters")

    cfg = dict(outdir=os.path.join(a.outdir, "results"),
               barcodes_dir=bcdir, barcodetag="CB", threads=8,
               binsize=50, normalize="RPKM",
               effective_genome_size=2913022398, blacklist="",
               merge_across_experiments=False, bams=bams)
    with open(os.path.join(a.outdir, "config.yaml"), "w") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)
    print(f"\nWrote {len(bams)} tables to {bcdir} and a config skeleton (fill in BAM paths).")


if __name__ == "__main__":
    main()