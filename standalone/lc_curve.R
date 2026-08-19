library(ggplot2)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(readr)
library(stringr)
library(purrr)

# 1. List all your parent directories
base_dirs <- c(
  "/date/gcb/gcb_MZ/nanoCTAR_YD/MS/"
  # Add more base paths as needed
)

# 2. Find all yield.txt files from all base directories
file_paths <- unlist(lapply(base_dirs, function(dir) {
  list.files(
    path = dir,
    pattern = "yield\\.txt$",
    recursive = TRUE,
    full.names = TRUE
  )
}))

cat("Found file paths:\n")
print(head(file_paths, 4))

# 3. Process and extract sample and assay info (adapt "FC_Droplet_Paired-Tag" logic if needed)
combined_data <- lapply(file_paths, function(f) {
  df <- read_delim(f, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  path_parts <- str_split(f, "/")[[1]]
  path_parts <- path_parts[path_parts != ""]
  # Find parent directory (adapt this name if you're searching others)
  parent_matches <- which(path_parts %in% basename(base_dirs))
  if (length(parent_matches) > 0) {
    match <- parent_matches[1]
    sample_folder <- path_parts[match + 1]
    assay_folder <- path_parts[match + 2]
    assay_type <- str_extract(assay_folder, "^[^_]+")
    df$Sample <- paste(sample_folder, assay_type, sep = "_")
    df$SampleFolder <- sample_folder
    df$AssayFolder <- assay_folder
    df$AssayType <- assay_type
  } else {
    df$Sample <- "UNKNOWN"
    df$SampleFolder <- "UNKNOWN"
    df$AssayFolder <- "UNKNOWN"
    df$AssayType <- "UNKNOWN"
  }
  df
}) %>% bind_rows()

# 4. Filter and prepare data
filtered_data <- combined_data %>%
  dplyr::filter(TOTAL_READS <= 500000000, Sample != "UNKNOWN") %>%
  mutate(
    TOTAL_READS_M = TOTAL_READS / 1e6,
    EXPECTED_DISTINCT_M = EXPECTED_DISTINCT / 1e6
  ) %>%
  arrange(Sample, TOTAL_READS_M)

filtered_data$Sample<- gsub("_yield\\.txt$", "", filtered_data$Sample)

# Let's say you want to plot only a specific subset of samples:
selected_samples <- c( #"FC_Droplet_Paired-Tag_rep1_H3K27ac",
                       #"FC_Droplet_Paired-Tag_rep1_H3K27me3",
                       #"FC_Droplet_Paired-Tag_rep2_H3K27ac","FC_Droplet_Paired-Tag_rep2_H3K27me3",
                       #"FC_Droplet_Paired-Tag_rep3_H3K27ac","FC_Droplet_Paired-Tag_rep3_H3K27me3",
                      #"nanoCTAR_e13_1_ATAC","nanoCTAR_e13_1_H3K27ac","nanoCTAR_e13_1_H3K27me3",
                      #"nanoCTAR_e13_2_ATAC","nanoCTAR_e13_2_H3K27ac","nanoCTAR_e13_2_H3K27me3")
                      #"nanoCTAR_e13_K4me3_1_ATAC","nanoCTAR_e13_K4me3_1_H3K4me3","nanoCTAR_e13_K4me3_1_H3K27me3",
                      #"nanoCTAR_e13_K4me3_2_ATAC","nanoCTAR_e13_K4me3_2_H3K4me3","nanoCTAR_e13_K4me3_2_H3K27me3",
                      "MS_1_H3K4me3","MS_1_H3K27me3",
                      "MS_2_H3K4me3","MS_2_H3K27me3")        
# Replace these names with exactly those you want to plot, matching entries in filtered_data$Sample


# Filter your data for those samples only
filtered_selected <- filtered_data %>% 
  dplyr::filter(Sample %in% selected_samples)

# Suppose your selected samples are as before
selected_samples <- unique(filtered_selected$Sample)

# Enter your actual observed total reads (raw counts or in millions, adjust as needed)
# EXAMPLE: These are in millions; adjust names and values to match your use case
actual_reads_input <- data.frame(
  Sample = selected_samples,
  ACTUAL_TOTAL_READS_M = c(#159.1,167.8,119.6,152.5,183.8,174.9,
                           15.3,5.9,88.7,31.3))
                           #4.8,16.8,2.6,5.5,16.3,1.7))
                           #36.7,91.3,14.1,29.9,79.6,15.3))
# If you have raw read counts (not millions), just keep as 'ACTUAL_TOTAL_READS' and don't divide by 1e6 below

# For each sample, interpolate (from its LC curve) the expected unique reads at your actual sequencing depth
library(purrr)

actual_points <- actual_reads_input %>%
  mutate(
    EXPECTED_DISTINCT_M = purrr::map2_dbl(
      Sample, ACTUAL_TOTAL_READS_M,
      function(s, x) {
        # Get the LC curve points for this sample
        df_sample <- filtered_selected %>% dplyr::filter(Sample == s)
        # Interpolate (linear) the expected unique reads at the given total reads
        approx(
          x = df_sample$TOTAL_READS_M,
          y = df_sample$EXPECTED_DISTINCT_M,
          xout = x, 
          rule = 2
        )$y
      }
    )
  )

# Use a more discreet marker, e.g., shape 21 (filled circle) or shape 17 (triangle)
actual_marker_shape <- 21  # Try 17 for triangle, or 21 for a nice filled circle

nr_samples <- 4
library(RColorBrewer)
palette <- colorRampPalette(brewer.pal(8, "Set2"))(nr_samples)


p <- ggplot(filtered_selected, aes(
  x = TOTAL_READS_M,
  y = EXPECTED_DISTINCT_M,
  color = Sample,
  group = Sample
)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 1.5) +
  # Discreet actual points: small, filled, subtle outline
  geom_point(
    data = actual_points,
    aes(x = ACTUAL_TOTAL_READS_M, y = EXPECTED_DISTINCT_M, fill = Sample),
    color = "white",         # White outline for subtle touch
    shape = actual_marker_shape,
    size = 3,
    stroke = 1,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(
    title = "LC curve: nanoCTAR vs scCut&Tag pro",
    x = "Total Reads (millions)",
    y = "Expected Distinct Reads (millions)"
  ) +
  theme_classic()

print(p)

# SAVE THE PLOT as a PNG
ggsave("/proj/user/mattia/Embryo_2/plots/LC_Curve_nanoCTAR_e13_modalities_scCut_Tag_pro.png", plot = p, width = 10, height = 7, dpi = 350)



# Update the palette for the number of selected samples
library(RColorBrewer)
nr_samples <- length(unique(filtered_selected$Sample))
palette <- colorRampPalette(brewer.pal(8, "Set2"))(nr_samples)



# Now make the plot with only the selected samples
p <- ggplot(filtered_selected, aes(
  x = TOTAL_READS_M,
  y = EXPECTED_DISTINCT_M,
  color = Sample,
  group = Sample
)) +
  geom_line(size = 1.2) +
  geom_point(size = 1.5) +
  scale_color_manual(values = palette) +
  labs(
    title = "Comparative LC Curves: Selected Samples",
    x = "Total Reads (millions)",
    y = "Expected Distinct Reads (millions)"
  ) +
  theme_classic()

print(p)

