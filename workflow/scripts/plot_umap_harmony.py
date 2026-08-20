"""
UMAP plots using the Harmony-integrated coordinates (MS_UMAP_coordinates_harmony.xlsx,
columns 'umapharmony_1' / 'umapharmony_2'), replacing the earlier non-integrated
UMAP_1/UMAP_2 coordinates. For each mark (H3K27ac, H3K27me3):

  1. colored by patient       (combined 2-panel + individual + grey-background highlight)
  2. colored by protocol      (fresh / frozen)
  3. colored by technology    (GEM-X / Next-GEM)

Metadata (sample, patient, technology, protocol) comes from MS_metadata.xlsx,
sheets 'MS|H3K27ac' / 'MS|H3K27me3', merged onto the harmony coordinates by 'cell'.

Style: theme_minimal (no bounding box, no ticks, faint background grid).
"""
import os
import pandas as pd
import matplotlib.pyplot as plt

OKABE_ITO_LIST = [
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#999999",  # grey
]


def save_publication_figure(fig, path, formats=("pdf", "png"), dpi=300):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    for fmt in formats:
        fig.savefig(f"{path}.{fmt}", dpi=dpi, bbox_inches="tight")

HARMONY_XLSX = "/date/gcb/gcb_wq/nanoCTAR_pipeline_barcode_correction/first_eval/MS_UMAP_coordinates.xlsx"
META_XLSX = "/date/gcb/gcb_wq/nanoCTAR_pipeline_barcode_correction/first_eval/MS_metadata.xlsx"
OUT_DIR = "/date/gcb/gcb_wq/nanoCTAR_pipeline_barcode_correction/figures_manual/umap"
MARKS = ["H3K27ac", "H3K27me3"]
SHEET = {"H3K27ac": "H3K27ac", "H3K27me3": "H3K27me3"}

PATIENT_ORDER = ["MS058CT1","MS058CT2", "MS381CT1", "MS381CT2","MS549CT","MS549_CT"]
PROTOCOL_ORDER = ["fresh", "frozen"]
TECHNOLOGY_ORDER = ["GEM-X", "Next-GEM"]

PATIENT_COLORS = dict(zip(PATIENT_ORDER, OKABE_ITO_LIST))
PROTOCOL_COLORS = dict(zip(PROTOCOL_ORDER, OKABE_ITO_LIST))
TECHNOLOGY_COLORS = dict(zip(TECHNOLOGY_ORDER, OKABE_ITO_LIST))

GREY = "#D9D9D9"


def apply_theme_minimal(ax):
    for spine in ["top", "right", "left", "bottom"]:
        ax.spines[spine].set_visible(False)
    ax.tick_params(length=0, labelbottom=False, labelleft=False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(True, color="#EBEBEB", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)


def load_mark(mark):
    coords = pd.read_excel(HARMONY_XLSX, sheet_name=SHEET[mark])
    meta = pd.read_excel(META_XLSX, sheet_name=SHEET[mark],
                          usecols=["cell", "sample", "patient", "technology", "protocol"])
    df = coords.merge(meta, on="cell", how="left")
    df = df.rename(columns={"umap_1": "UMAP_1", "umap_2": "UMAP_2"})
    assert df[["sample", "patient", "technology", "protocol"]].isna().sum().sum() == 0, \
        f"unmatched cells after merge for {mark}"
    return df


# ---------- grouped scatter (categorical color) ----------

def plot_grouped(ax, df, group_col, order, colors):
    present = [g for g in order if g in df[group_col].unique()]
    for g in present:
        sub = df[df[group_col] == g]
        ax.scatter(sub["UMAP_1"], sub["UMAP_2"], s=3, alpha=0.6, linewidths=0,
                   color=colors[g], label=f"{g} (n={len(sub)})", zorder=2)


def make_grouped_figures(dfs, group_col, order, colors, tag, legend_title):
    # combined 2-panel
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2))
    for ax, mark in zip(axes, MARKS):
        df = dfs[mark].sample(frac=1.0, random_state=0).reset_index(drop=True)
        plot_grouped(ax, df, group_col, order, colors)
        ax.set_title(mark)
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        apply_theme_minimal(ax)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, title=legend_title, loc="center left",
               bbox_to_anchor=(1.0, 0.5), frameon=False, markerscale=4)
    fig.suptitle(f"Harmony UMAP colored by {legend_title}", y=1.02)
    fig.tight_layout()
    save_publication_figure(fig, f"{OUT_DIR}/umap_harmony_by_{tag}_combined", formats=["pdf", "png"], dpi=300)
    plt.close(fig)

    # individual per-mark
    for mark in MARKS:
        df = dfs[mark].sample(frac=1.0, random_state=0).reset_index(drop=True)
        fig, ax = plt.subplots(figsize=(6, 5))
        plot_grouped(ax, df, group_col, order, colors)
        ax.set_title(mark)
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        apply_theme_minimal(ax)
        ax.legend(title=legend_title, loc="center left", bbox_to_anchor=(1.0, 0.5),
                  frameon=False, markerscale=4)
        fig.tight_layout()
        save_publication_figure(fig, f"{OUT_DIR}/umap_{mark}_by_{tag}",
                                 formats=["pdf", "png"], dpi=300)
        plt.close(fig)


# ---------- patient highlight (grey background) ----------

def highlight_panel(ax, df, patient):
    background = df[df["patient"] != patient]
    foreground = df[df["patient"] == patient]
    ax.scatter(background["UMAP_1"], background["UMAP_2"], s=3, c=GREY, alpha=0.5,
               linewidths=0, zorder=1)
    ax.scatter(foreground["UMAP_1"], foreground["UMAP_2"], s=3, c=PATIENT_COLORS[patient],
               alpha=0.85, linewidths=0, zorder=2)
    ax.set_title(f"{patient} (n={len(foreground)})")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    apply_theme_minimal(ax)


def make_highlight_figures(dfs):
    for mark in MARKS:
        df = dfs[mark]
        fig, axes = plt.subplots(1, len(PATIENT_ORDER), figsize=(4.2 * len(PATIENT_ORDER), 4.2))
        for ax, patient in zip(axes, PATIENT_ORDER):
            highlight_panel(ax, df, patient)
        fig.suptitle(f"{mark} — Harmony UMAP highlighted by patient", y=1.03)
        fig.tight_layout()
        save_publication_figure(fig, f"{OUT_DIR}/umap_{mark}_highlight_by_patient",
                                 formats=["pdf", "png"], dpi=300)
        plt.close(fig)

        for patient in PATIENT_ORDER:
            fig_i, ax_i = plt.subplots(figsize=(4.5, 4.2))
            highlight_panel(ax_i, df, patient)
            fig_i.tight_layout()
            save_publication_figure(fig_i, f"{OUT_DIR}/umap_{mark}_highlight_{patient}",
                                     formats=["pdf", "png"], dpi=300)
            plt.close(fig_i)


def main():
    dfs = {mark: load_mark(mark) for mark in MARKS}

    make_grouped_figures(dfs, "patient", PATIENT_ORDER, PATIENT_COLORS, "patient", "Patient")
    make_highlight_figures(dfs)
    make_grouped_figures(dfs, "protocol", PROTOCOL_ORDER, PROTOCOL_COLORS, "protocol", "Protocol")
    make_grouped_figures(dfs, "technology", TECHNOLOGY_ORDER, TECHNOLOGY_COLORS, "technology", "Technology")

    for mark, df in dfs.items():
        print(f"--- {mark} ---")
        print("patient:", df["patient"].value_counts().to_dict())
        print("protocol:", df["protocol"].value_counts().to_dict())
        print("technology:", df["technology"].value_counts().to_dict())


if __name__ == "__main__":
    main()