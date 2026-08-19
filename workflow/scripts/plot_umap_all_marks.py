"""
UMAP plots for all 3 marks (ATAC, H3K27ac, H3K27me3), both the non-integrated
coordinates (PBMCs_UMAP_coordinates.xlsx, umap_1/umap_2) and the
Harmony-integrated coordinates (PBMCs_UMAP_coordinates_harmony.xlsx,
umapharmony_1/umapharmony_2). For each embedding:

  1. colored by sample
  2. colored by cluster        (seurat_clusters for normal / harmony_clusters for Harmony)
  3. colored by major cell type (T cells / B cells / NK cells / Monocytes /
     other), collapsed from the fine predicted.id_ATLAS_ATAC atlas labels.

The predicted.id_ATLAS_ATAC label is only stored on the PBMC|ATAC sheet, so it
is transferred onto the H3K27ac / H3K27me3 cells that share the same physical
nucleus (sample|patient|barcode), mirroring the Rmd's propagate_atlas_atac
chunk. Cells NOT present in the ATAC library are dropped from the H3K27ac /
H3K27me3 panels entirely (per request), so every plotted H3K27ac / H3K27me3
cell has a real ATAC-derived label.

Style: theme_minimal (no bounding box, no ticks, faint background grid).
"""
import os
import pandas as pd
import matplotlib.pyplot as plt

# ---- fine atlas label -> major lineage ----
MAJOR_MAP = {
    "T_naive": "T cells", "T_mem_CD4": "T cells", "T_mem_CD8": "T cells", "T_mix": "T cells",
    "B_naive": "B cells", "B_mem": "B cells", "Plasma": "B cells",
    "NK_CD16": "NK cells",
    "Mono_CD14": "Monocytes", "Mono_CD16": "Monocytes", "Mono_other": "Monocytes",
}
MAJOR_ORDER = ["T cells", "B cells", "NK cells", "Monocytes", "other"]

# Explicit, distinct, colorblind-safe colors for the major lineages
MAJOR_COLORS = {
    "T cells":   "#0072B2",   # blue
    "B cells":   "#E69F00",   # orange
    "NK cells":  "#009E73",   # bluish green
    "Monocytes": "#D55E00",   # vermillion
    "other":     "#999999",   # grey
}

# A larger, high-contrast qualitative palette for many-category groupings
# (e.g. 16-23 Seurat / Harmony clusters). Built from tab20 + tab20b + tab20c.
def _extended_palette(n):
    base = []
    for name in ("tab20", "tab20b", "tab20c"):
        cmap = plt.get_cmap(name)
        base.extend([cmap(i) for i in range(cmap.N)])
    if n > len(base):
        base = (base * (n // len(base) + 1))
    return base[:n]


def qualitative_colors(categories):
    cats = list(categories)
    pal = _extended_palette(len(cats))
    return {c: pal[i] for i, c in enumerate(cats)}


def save_publication_figure(fig, path, formats=("pdf", "png"), dpi=300):
    for fmt in formats:
        fig.savefig(f"{path}.{fmt}", dpi=dpi, bbox_inches="tight")


COORD_XLSX = "/project/PBMCs_UMAP_coordinates.xlsx"
HARMONY_XLSX = "/project/PBMCs_UMAP_coordinates_harmony.xlsx"
META_XLSX = "/project/PBMCs_metadata.xlsx"
OUT_DIR = "/project/figures/umap_all_marks"
MARKS = ["ATAC", "H3K27ac", "H3K27me3"]
SHEET = {"ATAC": "PBMC|ATAC", "H3K27ac": "PBMC|H3K27ac", "H3K27me3": "PBMC|H3K27me3"}


def apply_theme_minimal(ax):
    for spine in ["top", "right", "left", "bottom"]:
        ax.spines[spine].set_visible(False)
    ax.tick_params(length=0, labelbottom=False, labelleft=False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(True, color="#EBEBEB", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)


def build_match_key(df):
    return df["sample"].astype(str) + "|" + df["patient"].astype(str) + "|" + df["barcode"].astype(str)


def get_atac_label():
    """predicted.id_ATLAS_ATAC keyed by sample|patient|barcode, from the ATAC sheet."""
    md = pd.read_excel(META_XLSX, sheet_name=SHEET["ATAC"],
                       usecols=["sample", "patient", "barcode", "predicted.id_ATLAS_ATAC"])
    md["match_key"] = build_match_key(md)
    assert md["match_key"].is_unique, "duplicate ATAC match keys"
    return md.set_index("match_key")["predicted.id_ATLAS_ATAC"]


def load_mark(mark, coord_xlsx, x_col, y_col, cluster_col, atac_label):
    coords = pd.read_excel(coord_xlsx, sheet_name=SHEET[mark])
    coords = coords.rename(columns={x_col: "UMAP_1", y_col: "UMAP_2"})
    meta = pd.read_excel(META_XLSX, sheet_name=SHEET[mark],
                          usecols=["cell", "sample", "patient", "barcode", cluster_col])
    df = coords.merge(meta, on="cell", how="left")
    df["match_key"] = build_match_key(df)
    df["predicted.id_ATLAS_ATAC"] = df["match_key"].map(atac_label)
    # drop cells not present in ATAC (H3K27ac / H3K27me3); ATAC keeps everything
    if mark != "ATAC":
        df = df[df["predicted.id_ATLAS_ATAC"].notna()].copy()
    df["major_celltype"] = df["predicted.id_ATLAS_ATAC"].map(MAJOR_MAP).fillna("other")
    return df


# ---------- grouped scatter (categorical color) ----------

def plot_grouped(ax, df, group_col, order, colors):
    present = [g for g in order if g in df[group_col].unique()]
    for g in present:
        sub = df[df[group_col] == g]
        ax.scatter(sub["UMAP_1"], sub["UMAP_2"], s=3, alpha=0.6, linewidths=0,
                   color=colors[g], label=f"{g} (n={len(sub)})", zorder=2)


def make_grouped_figures(dfs, group_col, order, colors, tag, legend_title, embedding_tag):
    fig, axes = plt.subplots(1, len(MARKS), figsize=(4.5 * len(MARKS), 4.2))
    for ax, mark in zip(axes, MARKS):
        df = dfs[mark].sample(frac=1.0, random_state=0).reset_index(drop=True)
        plot_grouped(ax, df, group_col, order, colors)
        ax.set_title(f"{mark} (n={len(df)})")
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        apply_theme_minimal(ax)
    # pool legend entries across all panels (clusters differ per modality), one
    # swatch per category in `order`, colored consistently
    seen = set()
    handles, labels = [], []
    for ax in axes:
        for h, l in zip(*ax.get_legend_handles_labels()):
            key = l.split(" (n=")[0]
            if key not in seen:
                seen.add(key)
                handles.append(h)
                labels.append(key)
    # reorder to match `order`
    order_idx = {str(g): i for i, g in enumerate(order)}
    paired = sorted(zip(labels, handles), key=lambda t: order_idx.get(t[0], len(order_idx)))
    labels, handles = [p[0] for p in paired], [p[1] for p in paired]
    fig.legend(handles, labels, title=legend_title, loc="center left",
               bbox_to_anchor=(1.0, 0.5), frameon=False, markerscale=4, ncol=1)
    fig.suptitle(f"{embedding_tag} UMAP colored by {legend_title}", y=1.02)
    fig.tight_layout()
    save_publication_figure(fig, f"{OUT_DIR}/umap_{tag}_combined", formats=["pdf", "png"], dpi=300)
    plt.close(fig)

    for mark in MARKS:
        df = dfs[mark].sample(frac=1.0, random_state=0).reset_index(drop=True)
        fig, ax = plt.subplots(figsize=(6, 5))
        plot_grouped(ax, df, group_col, order, colors)
        ax.set_title(f"{mark} (n={len(df)})")
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        apply_theme_minimal(ax)
        ax.legend(title=legend_title, loc="center left", bbox_to_anchor=(1.0, 0.5),
                  frameon=False, markerscale=4)
        fig.tight_layout()
        save_publication_figure(fig, f"{OUT_DIR}/umap_{mark}_{tag}", formats=["pdf", "png"], dpi=300)
        plt.close(fig)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    atac_label = get_atac_label()

    embeddings = [
        ("Normal", "normal", COORD_XLSX, "umap_1", "umap_2", "seurat_clusters"),
        ("Harmony", "harmony", HARMONY_XLSX, "umapharmony_1", "umapharmony_2", "harmony_clusters"),
    ]

    for embedding_tag, tag_prefix, coord_xlsx, x_col, y_col, cluster_col in embeddings:
        dfs = {mark: load_mark(mark, coord_xlsx, x_col, y_col, cluster_col, atac_label)
               for mark in MARKS}

        sample_order = sorted(pd.concat([dfs[m]["sample"] for m in MARKS]).dropna().unique())
        sample_colors = qualitative_colors(sample_order)

        cluster_order = sorted(pd.concat([dfs[m][cluster_col] for m in MARKS]).dropna().unique())
        cluster_colors = qualitative_colors(cluster_order)

        major_order = [c for c in MAJOR_ORDER
                       if c in pd.concat([dfs[m]["major_celltype"] for m in MARKS]).unique()]

        make_grouped_figures(dfs, "sample", sample_order, sample_colors,
                              f"{tag_prefix}_by_sample", "Sample", embedding_tag)
        make_grouped_figures(dfs, cluster_col, cluster_order, cluster_colors,
                              f"{tag_prefix}_by_cluster", "Cluster", embedding_tag)
        make_grouped_figures(dfs, "major_celltype", major_order, MAJOR_COLORS,
                              f"{tag_prefix}_by_major_celltype", "Major cell type", embedding_tag)

        for mark, df in dfs.items():
            print(f"[{embedding_tag}] {mark}: {len(df)} cells | "
                  + ", ".join(f"{k}={v}" for k, v in df['major_celltype'].value_counts().items()))


if __name__ == "__main__":
    main()