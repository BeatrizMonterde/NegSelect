# NegSelect

Statistical method to prove negative selection in the context of cancer and make comparisons between groups (e.g. primary vs metastases). 

Required information:

- Size of the immunopeptidome (bp), generated with SOPRANO (https://github.com/instituteofcancerresearch/SOPRANO)
- Length of genome covered by the sequencing reaction (bp)
- VCF annotated files: total number of somatic mutations and mutations affecting the immunopeptidome 

library(dplyr)
library(stats)

## Function to set working directory, read files, and perform calculations
process_immunopeptidome <- function(immunopeptidome_file, primary_file, metastasis_file, annotated_file) {
  
  # Read the immunopeptidome bed file
  immunopeptidome <- read.delim(immunopeptidome_file, sep = "\t", header = FALSE, col.names = c("Transcript_ID", "Start", "End"))
  
  # Calculate the number of immunogenic amino acids
  Im <- sum(immunopeptidome$End - immunopeptidome$Start)
  
  # Calculate the number of immunogenic base pairs
  Im_bp <- Im * 3
  
  # Calculate the mutation rate WGS
  p <- Im_bp / (3 * 10^9)
  
  # Read the filtered data for primary and metastasis samples
  primary <- read.csv(primary_file, sep = "\t", header = TRUE)
  metastasis <- read.csv(metastasis_file, sep = "\t", header = TRUE)
  
  # Calculate the number of rows for primary and metastasis samples
  nrows_primary <- nrow(primary)
  nrows_metastasis <- nrow(metastasis)
  
  # Calculate the expected values for primary and metastasis samples
  exp_primary <- p * nrows_primary
  exp_metastasis <- p * nrows_metastasis
  
  # Read the SOPRANO results
  annot <- read.delim(annotated_file, sep = "\t", header = FALSE)
  
  # Extract the observed values for primary and metastasis samples
  obs_primary <- annot[annot$V1 == "Primary", "V2"]
  obs_metastasis <- annot[annot$V1 == "Metastasis", "V2"]
  
  # Calculate the ppois value for primary and metastasis samples
  ppois_primary <- ppois(obs_primary, lambda = exp_primary)
  ppois_metastasis <- ppois(obs_metastasis, lambda = exp_metastasis)
  
  # Create a final dataframe to store the results
  final_df <- data.frame(
    Sample = c("Primary", "Metastasis"),
    Size_Im = rep(Im_bp, 2),  # Add Im_bp column
    p_Im = rep(p, 2),      # Add p column
    Observed = c(obs_primary, obs_metastasis),
    Expected = c(exp_primary, exp_metastasis),
    ppois = c(ppois_primary, ppois_metastasis)
  )
  
  return(final_df)
}

# Ideally, you should construct a final_df with this structure 
> head(final_df)
                    Sample   Im_bp        p_Im Total_muts Observed Expected     ppois
1        ES0001_Metastasis 1606785 0.000535595       3538        2 1.894935 0.4351901
2           ES0001_Primary 1606785 0.000535595       2608        1 1.396832 0.5929270
3 ES0002_Metastasis_Lung_1 1562901 0.000520967       2649        1 1.380042 0.5987425
4 ES0002_Metastasis_Lung_2 1562901 0.000520967       3129        0 1.630106 0.9169827
5 ES0002_Metastasis_Lung_3 1562901 0.000520967       2748        1 1.431617 0.5809674
6 ES0002_Metastasis_Lung_4 1562901 0.000520967       2677        1 1.394629 0.2479251

# Combinatorial analysis and Poisson test
library(purrr)

# Function to process data for primary or metastasis samples
process_samples <- function(final_results_file, sample_type, output_file) {
  
  # Read the data
  final_results <- read.delim(final_results_file, sep = "\t")
  
  # Filter based on sample type (Primary or Metastasis)
  samples_data <- final_results %>% filter(grepl(sample_type, Sample))
  
  mean(samples_data$Observed, na.rm = TRUE)
  mean(samples_data$Expected, na.rm = TRUE)
  
  # Function to calculate ppois for a given combination of samples
  calculate_ppois <- function(samples, data) {
    sample_data <- data %>% filter(Sample %in% samples)
    
    avg_Im <- mean(sample_data$Im_bp)
    total_muts <- sum(sample_data$Total_muts)
    p_Im_combined <- avg_Im / (3 * 10^9)
    expected_muts <- p_Im_combined * total_muts
    observed_muts <- sum(sample_data$Observed)
    
    ppois_value <- ppois(observed_muts, lambda = expected_muts, lower.tail = TRUE)
    
    return(data.frame(Cohort_Size = length(samples), Samples = paste(samples, collapse = ", "), 
                      Avg_Im = avg_Im, Total_Muts = total_muts, Expected_Muts = expected_muts,
                      Observed_Muts = observed_muts, ppois = ppois_value))
  }
  
  # Initialize a list to store results
  results <- list()
  
  # Determine the number of samples
  num_samples <- nrow(samples_data)
  
  # Iterate over all cohort sizes from n=1 to the number of samples
  for (n in 1:num_samples) {
    combinations <- combn(samples_data$Sample, n, simplify = FALSE)
    cohort_results <- map_df(combinations, calculate_ppois, data = samples_data)
    results <- append(results, list(cohort_results))
  }
  
  # Combine all results into a single dataframe
  results_df <- bind_rows(results)
  
  mean(results_df$Observed_Muts, na.rm = TRUE)
  mean(results_df$Expected_Muts, na.rm = TRUE)
  
  # Save results to a text file
  write.table(results_df, output_file, row.names = FALSE, quote = FALSE, sep = "\t")
}

# SISMO: generate somatic mutations based on the trinucleotide context 
https://github.com/luisgls/SISMO

# Generate ~ 100 simulations per sample 

# Poisson graph comparing 2 groups 
# Create a final_df containing the information from primary, metastasis and their null distributions

final_df <- bind_rows(primary, metastasis, null_P, null_M)

# Transform Y values to -log(Y)
final_df$neglog_Y <- -log(final_df2$Y)

# Create the ggplot
ggplot(final_df_v2, aes(x = X, y = neglog_Y, color = factor(Site))) +
  geom_hline(yintercept = -log(0.05), linetype = "dashed", color = "grey24") +
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun = mean, geom = "line", aes(group = factor(Site))) +
  labs(x = "Cohort size", y = "-log(ppois)", title = "Negative selection", color = "Site") +
  +
  scale_color_manual(values = c("Primary" = "color1", "Metastasis" = "color2", "Null_P" = "color3", "Null_M" = "color4")) 

