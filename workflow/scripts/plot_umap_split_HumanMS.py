"""
Two standalone (split) UMAP plots -- one for H3K27ac, one for H3K27me3 --
colored by the Human_MS reference label-transfer annotation
('predicted.id_Human_MS'), using the same color palette and theme_minimal
styling as the combined/linked figure (umap_H3K27ac_H3K27me3_linked_HumanMS).

H3K27me3 only has this annotation via label transfer from its paired
H3K27ac cell (same sample + patient + 10x barcode); cells without a
matching H3K27ac barcode are shown in grey as "unpaired".
"""
import os
import pandas as pd
import matplotlib.pyplot as plt


def save_publication_figure(fig, path, formats=("pdf", "png"), dpi=300):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    for fmt in formats:
        fig.savefig(f"{path}.{fmt}", dpi=dpi, bbox_inches="tight")

COORD_XLSX = "/date/gcb/gcb_wq/nanoCTAR_pipeline_barcode_correction/first_eval/MS_UMAP_coordinates_harmony.xlsx"
META_XLSX = "/date/gcb/gcb_wq/nanoCTAR_pipeline_barcode_correction/first_eval/MS_metadata.xlsx"
OUT_DIR = "/date/gcb/gcb_wq/nanoCTAR_pipeline_barcode_correction/predicted.id_Human_MS_harmony/"
RNG_SEED = 0

# Same explicit color mapping as the linked figure (see plot_umap_link_barcodes.py):
# Okabe-Ito based, with "OPCs + COPs" (blue) and "Oligodendrocytes" (vermillion)
# deliberately placed at opposite ends of the color wheel for clear separation.
CELLTYPE_COLORS = {
    "Astrocytes":          "#CC79A7",  # reddish purple
    "B cells":             "#332288",  # indigo
    "Endo + Peri":         "#000000",  # black
    "Excitatory neurons":  "#009E73",  # bluish green
    "Inhibitory neurons":  "#999933",  # olive
    "Microglia":           "#56B4E9",  # sky blue
    "OPCs + COPs":         "#0072B2",  # blue
    "Oligodendrocytes":    "#D55E00",  # vermillion
    "T cells":             "#F0E442",  # yellow
}


def apply_theme_minimal(ax):
    for spine in ["top", "right", "left", "bottom"]:
        ax.spines[spine].set_visible(False)
    ax.tick_params(length=0, labelbottom=False, labelleft=False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(True, color="#EBEBEB", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)


def load_mark(mark, meta_cols):
    umap = pd.read_excel(COORD_XLSX, sheet_name=mark)
    meta = pd.read_excel(META_XLSX, sheet_name=mark, usecols=meta_cols)
    return umap.merge(meta, on="cell", how="left")


def plot_split(df, mark, categories, colors):
    # drop cells with no Human_MS label (H3K27me3 cells with no matching
    # H3K27ac barcode, i.e. "unpaired") -- not shown in this plot
    df_plot = df[df["predicted.id_Human_MS"].notna()]
    df_plot = df_plot.sample(frac=1.0, random_state=RNG_SEED).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(6, 5))

    for ct in categories:
        sub = df_plot[df_plot["predicted.id_Human_MS"] == ct]
        ax.scatter(sub["umapharmony_1"], sub["umapharmony_2"], s=3, alpha=0.75, linewidths=0,
                   color=colors[ct], label=f"{ct} (n={len(sub)})", zorder=2)

    ax.set_title(mark)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    apply_theme_minimal(ax)
    ax.legend(title="predicted.id_Human_MS", loc="center left", bbox_to_anchor=(1.0, 0.5),
              frameon=False, markerscale=4)
    fig.tight_layout()

    save_publication_figure(fig, f"{OUT_DIR}/umap_{mark}_HumanMS_split", formats=["pdf", "png"], dpi=300)
    plt.close(fig)


def main():
    ac = load_mark("H3K27ac", ["cell", "sample", "patient", "barcode", "predicted.id_Human_MS"])
    me = load_mark("H3K27me3", ["cell", "sample", "patient", "barcode"])

    ac["match_key"] = ac["sample"] + "|" + ac["patient"] + "|" + ac["barcode"].astype(str)
    me["match_key"] = me["sample"] + "|" + me["patient"] + "|" + me["barcode"].astype(str)

    label_map = ac.set_index("match_key")["predicted.id_Human_MS"]
    me["predicted.id_Human_MS"] = me["match_key"].map(label_map)

    # same category set / color mapping as the combined linked figure
    categories = sorted(ac["predicted.id_Human_MS"].dropna().unique())
    colors = {ct: CELLTYPE_COLORS[ct] for ct in categories}

    plot_split(ac, "H3K27ac", categories, colors)
    plot_split(me, "H3K27me3", categories, colors)

    print("H3K27ac:", ac["predicted.id_Human_MS"].value_counts().to_dict())
    print("H3K27me3 (transferred):", me["predicted.id_Human_MS"].value_counts(dropna=False).to_dict())


if __name__ == "__main__":
    main()