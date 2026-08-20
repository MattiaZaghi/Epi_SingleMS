"""
Two-panel UMAP (H3K27ac | H3K27me3), colored by the Human_MS reference
label-transfer annotation (metadata column 'predicted.id_Human_MS', only
available for H3K27ac). Since H3K27ac and H3K27me3 cells sharing the same
sample + patient + 10x barcode are the same physical nucleus profiled for
both marks, that Human_MS label is transferred onto the matching H3K27me3
cell for coloring, and a thin line connects each matched pair across panels.

Style: theme_minimal (no bounding box, no ticks, faint background grid).
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch


def save_publication_figure(fig, path, formats=("pdf", "png"), dpi=300):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    for fmt in formats:
        fig.savefig(f"{path}.{fmt}", dpi=dpi, bbox_inches="tight")


COORD_XLSX = "/date/gcb/gcb_wq/nanoCTAR_pipeline_barcode_correction/first_eval/MS_UMAP_coordinates_harmony.xlsx"
META_XLSX = "/date/gcb/gcb_wq/nanoCTAR_pipeline_barcode_correction/first_eval/MS_metadata.xlsx"
OUT_DIR = "/date/gcb/gcb_wq/nanoCTAR_pipeline_barcode_correction/figures_manual/umap_by_patient"

N_LINKS = 1500          # subsample of connecting lines (avoid a hairball)
RNG_SEED = 0

LINK_COLOR = "#666666"

# Explicit, colorblind-safe cell-type -> color mapping (Okabe-Ito based, with
# two extra hues for the 9th/10th category). Chosen so that "OPCs + COPs" and
# "Oligodendrocytes" -- both oligodendrocyte-lineage populations that were
# nearly indistinguishable muted red/wine tones in the previous palette -- are
# now maximally separated (blue vs. vermillion, opposite ends of the wheel).
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


SHEET = {"H3K27ac": "H3K27ac", "H3K27me3": "H3K27me3"}


def load_mark(mark, meta_cols):
    umap = pd.read_excel(COORD_XLSX, sheet_name=SHEET[mark])
    all_meta_cols = ["cell", "sample", "patient"] + meta_cols
    meta = pd.read_excel(META_XLSX, sheet_name=SHEET[mark], usecols=all_meta_cols)
    df = umap.merge(meta, on="cell", how="left")
    return df


def main():
    ac = load_mark("H3K27ac", ["barcode", "predicted.id_Human_MS"])
    me = load_mark("H3K27me3", ["barcode"])

    # Key identifying the same physical nucleus across both marks
    ac["match_key"] = ac["sample"] + "|" + ac["patient"] + "|" + ac["barcode"].astype(str)
    me["match_key"] = me["sample"] + "|" + me["patient"] + "|" + me["barcode"].astype(str)

    assert ac["match_key"].is_unique, "duplicate AC match keys"
    assert me["match_key"].is_unique, "duplicate ME match keys"

    # Transfer the Human_MS label from H3K27ac onto its paired H3K27me3 cell
    label_map = ac.set_index("match_key")["predicted.id_Human_MS"]
    me["predicted.id_Human_MS"] = me["match_key"].map(label_map)
    me["paired"] = me["predicted.id_Human_MS"].notna()

    n_paired = me["paired"].sum()
    print(f"H3K27ac cells: {len(ac)}  |  H3K27me3 cells: {len(me)}  |  paired (shared barcode): {n_paired}")

    categories = sorted(ac["predicted.id_Human_MS"].dropna().unique())
    colors = {ct: CELLTYPE_COLORS[ct] for ct in categories}

    # shuffle plotting order within each panel (avoid systematic overplotting)
    ac_plot = ac.sample(frac=1.0, random_state=RNG_SEED).reset_index(drop=True)
    me_plot = me.sample(frac=1.0, random_state=RNG_SEED).reset_index(drop=True)

    fig, (ax_ac, ax_me) = plt.subplots(1, 2, figsize=(11, 5))

    # --- left panel: H3K27ac, colored by Human_MS label ---
    for ct in categories:
        sub = ac_plot[ac_plot["predicted.id_Human_MS"] == ct]
        ax_ac.scatter(sub["umapharmony_1"], sub["umapharmony_2"], s=3, alpha=0.75, linewidths=0,
                      color=colors[ct], label=f"{ct} (n={len(sub)})", zorder=2)
    ax_ac.set_title("H3K27ac")
    ax_ac.set_xlabel("UMAP 1")
    ax_ac.set_ylabel("UMAP 2")
    apply_theme_minimal(ax_ac)

    # --- right panel: H3K27me3, colored by transferred Human_MS label ---
    # (unpaired cells, i.e. no matching H3K27ac barcode, are dropped from the plot)
    me_plot = me_plot[me_plot["paired"]]
    for ct in categories:
        sub = me_plot[me_plot["predicted.id_Human_MS"] == ct]
        ax_me.scatter(sub["umapharmony_1"], sub["umapharmony_2"], s=3, alpha=0.75, linewidths=0,
                      color=colors[ct], label=f"{ct} (n={len(sub)})", zorder=2)
    ax_me.set_title("H3K27me3")
    ax_me.set_xlabel("UMAP 1")
    ax_me.set_ylabel("UMAP 2")
    apply_theme_minimal(ax_me)

    # shared legend (categories only)
    handles, labels = ax_ac.get_legend_handles_labels()
    fig.legend(handles, labels, title="predicted.id_Human_MS", loc="center left",
               bbox_to_anchor=(1.0, 0.5), frameon=False, markerscale=4)

    # --- connecting lines between paired cells (subsampled) ---
    paired = me[me["paired"]][["match_key", "umapharmony_1", "umapharmony_2"]].rename(
        columns={"umapharmony_1": "umapharmony_1_me", "umapharmony_2": "umapharmony_2_me"}
    )
    ac_xy = ac.set_index("match_key")[["umapharmony_1", "umapharmony_2"]].rename(
        columns={"umapharmony_1": "umapharmony_1_ac", "umapharmony_2": "umapharmony_2_ac"}
    )
    pairs = paired.join(ac_xy, on="match_key").dropna()

    rng = np.random.default_rng(RNG_SEED)
    n_show = min(N_LINKS, len(pairs))
    sample_idx = rng.choice(len(pairs), size=n_show, replace=False)
    pairs_sub = pairs.iloc[sample_idx]

    for _, row in pairs_sub.iterrows():
        con = ConnectionPatch(
            xyA=(row["umapharmony_1_ac"], row["umapharmony_2_ac"]), coordsA=ax_ac.transData,
            xyB=(row["umapharmony_1_me"], row["umapharmony_2_me"]), coordsB=ax_me.transData,
            color=LINK_COLOR, linewidth=0.25, alpha=0.15, zorder=1,
        )
        fig.add_artist(con)

    fig.suptitle(
        f"Paired scCUT&Tag UMAPs linked by shared barcode "
        f"(showing {n_show:,} of {len(pairs):,} matched cells)",
        y=1.03,
    )
    fig.tight_layout()

    save_publication_figure(fig, f"{OUT_DIR}/umap_H3K27ac_H3K27me3_linked_HumanMS",
                             formats=["pdf", "png"], dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()