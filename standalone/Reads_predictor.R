library(ggplot2)
library(dplyr)

# Function to predict unique reads for a single modality
predict_unique_reads <- function(sample_data, target_reads) {
  # Use linear interpolation for prediction
  approx(x = sample_data$TOTAL_READS, 
         y = sample_data$EXPECTED_DISTINCT, 
         xout = target_reads, 
         rule = 2)$y  # rule=2 allows extrapolation
}

# Function to predict combined unique reads for multiple modalities
predict_multimodal_unique <- function(filtered_data, sample_prefix, 
                                      total_target_reads, proportions = NULL) {
  
  # Extract all modalities for this sample
  sample_modalities <- filtered_data %>%
    filter(str_detect(Sample, paste0("^", sample_prefix)))
  
  unique_modalities <- unique(sample_modalities$Sample)
  
  # If no proportions provided, use equal distribution
  if (is.null(proportions)) {
    proportions <- rep(1/length(unique_modalities), length(unique_modalities))
  }
  
  # Ensure proportions sum to 1
  proportions <- proportions / sum(proportions)
  
  # Calculate target reads per modality
  target_reads_per_modality <- total_target_reads * proportions
  
  # Predict unique reads for each modality
  predictions <- data.frame(
    Modality = unique_modalities,
    Target_Reads = target_reads_per_modality,
    Predicted_Unique = NA
  )
  
  for (i in seq_along(unique_modalities)) {
    modality_data <- sample_modalities %>%
      filter(Sample == unique_modalities[i])
    
    predictions$Predicted_Unique[i] <- predict_unique_reads(
      modality_data, target_reads_per_modality[i]
    )
  }
  
  # Calculate total predicted unique reads
  # Note: This assumes independence between modalities
  total_unique <- sum(predictions$Predicted_Unique, na.rm = TRUE)
  
  return(list(
    modality_predictions = predictions,
    total_unique_reads = total_unique,
    total_target_reads = total_target_reads,
    efficiency = total_unique / total_target_reads
  ))
}

# Example usage:
# Predict for nanoCTAR_e13_1 sample with 50M total reads
sample_prediction <- predict_multimodal_unique(
  filtered_data = filtered_data,
  sample_prefix = "nanoCTAR_e13_2",
  total_target_reads = 200e6,
  proportions = c(0.193, 0.706, 0.101)  # Custom proportions for 3 modalities
)

print(sample_prediction)
