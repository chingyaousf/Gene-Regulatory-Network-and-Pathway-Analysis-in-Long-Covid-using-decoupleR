BiocManager::install('saezlab/decoupleR')

if (!requireNamespace("OmnipathR", quietly = TRUE)) {
    BiocManager::install("OmnipathR")
}


# Load required libraries
library(Seurat)
library(dplyr)
library(decoupleR)

## for specific cell type
# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Inspect the Seurat object metadata
head(seurat_object@meta.data)
unique(seurat_object@meta.data$disease)
unique(seurat_object@meta.data$CT_Sub_combined_02)

# Subset to one cell type (e.g., CD4.Naive)
cell_type <- "CD4.Tfh"
current_cell_type <- subset(seurat_object, CT_Sub_combined_02 == cell_type)

# Further subset to Long Covid samples
longcovid_cells <- subset(current_cell_type, disease == "LongCovid")

# Extract normalized log-transformed gene expression matrix
mat <- as.matrix(longcovid_cells@assays$RNA@data)

# Load CollecTRI regulons (human)
net <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)

# Run Univariate Linear Model (ULM) method for Long Covid
ulm_results <- decoupleR::run_ulm(
  mat = mat,
  net = net,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 5
)

# View a preview of the results
print(head(ulm_results))

# Save results for future analysis
saveRDS(
  ulm_results, 
  file = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/TF_longcovid_", cell_type, ".rds")
)

message(paste("TF activity analysis completed for Long Covid in cell type:", cell_type))



## for all cell types
# Load required libraries
library(Seurat)
library(dplyr)
library(decoupleR)

# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Inspect the Seurat object metadata
head(seurat_object@meta.data)
unique(seurat_object@meta.data$disease)
unique(seurat_object@meta.data$CT_Sub_combined_02)


# Define the cell types to analyze
cell_types_to_analyze <- c("B_naive", "B_exhausted", "B_immature", "B_malignant", 
                           "B_memory", "CD16_mono", "CD4.Naive", "CD4.CM", "CD4.EM", 
                           "CD4.IL22", "CD4.Tfh", "CD8.Naive", "CD8.EM", "CD8.TE", 
                           "CD14_mono", "DC_cells", "HSC", "Lymphocytes", "MAIT", 
                           "NKT", "NK", "Plasma", "Platelets", "RBC", "Treg", "gdT", "pDC")

# Load CollecTRI regulons (human)
net <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)

# Define groups (LongCovid and Control)
groups <- c("LongCovid", "Control")

# Loop over each cell type and group (LongCovid and Control)
for (cell_type in cell_types_to_analyze) {
  for (group in groups) {
    # Subset data for the current cell type and group
    current_cell_type <- subset(seurat_object, CT_Sub_combined_02 == cell_type)
    current_group_cells <- subset(current_cell_type, disease == group)
    
    # Check if there are cells in this subset
    if (ncol(current_group_cells) == 0) {
      message(paste("No cells for", cell_type, "in", group))
      next
    }
    
    # Extract normalized log-transformed gene expression matrix
    mat <- as.matrix(current_group_cells@assays$RNA@data)
    
    # Run Univariate Linear Model (ULM) method
    ulm_results <- decoupleR::run_ulm(
      mat = mat,
      net = net,
      .source = "source",
      .target = "target",
      .mor = "mor",
      minsize = 5
    )
    
    # Save results for future analysis
    save_path <- paste0(
      "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/TF_", 
      tolower(group), "_", cell_type, ".rds"
    )
    saveRDS(ulm_results, file = save_path)
    
    message(paste("TF activity analysis completed for", cell_type, "in", group, "and saved to", save_path))
  }
}



# Load required libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)

# Define the file path for RDS files
path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/"

# List all RDS files
files <- list.files(path, pattern = "\\.rds$", full.names = TRUE)

# Initialize an empty list to store TF activities
tf_activity_list <- list()

# Loop over RDS files to extract and combine results
for (file in files) {
  # Load the RDS file
  ulm_results <- readRDS(file)
  
  # Pivot to wide format for each TF
  ulm_wide <- ulm_results %>%
    tidyr::pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
    tibble::column_to_rownames('source')
  
  # Store in the list (keyed by file name)
  tf_activity_list[[basename(file)]] <- ulm_wide
}

# Combine all results into a single matrix
tf_combined_matrix <- do.call(cbind, tf_activity_list)

# Check dimensions of the combined matrix
print(dim(tf_combined_matrix))

# Save the combined matrix as an RDS file
output_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/tf_combined_matrix.rds"
saveRDS(tf_combined_matrix, file = output_file)

# Confirm save location
message("Combined matrix saved at: ", output_file)



# Check the first few rows of the combined matrix
head(tf_combined_matrix)

# Optionally, check column names (conditions)
colnames(tf_combined_matrix)

# Optionally, check row names (TF sources)
rownames(tf_combined_matrix)

# Display the first 5 rows and the first 5 columns
tf_combined_matrix[1:5, 1:5]

# Use the colnames() function to check the actual column names:
colnames(tf_combined_matrix)[1:5]




# Load required libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(pheatmap)

# Define the file path for the saved combined matrix
matrix_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/tf_combined_matrix.rds"

# Load the combined matrix from the saved RDS file
tf_combined_matrix <- readRDS(matrix_file)

# Calculate mean TF activity per cell type and condition
mean_activity <- tf_combined_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "TF") %>%
  pivot_longer(cols = -TF, names_to = "sample", values_to = "activity") %>%
  mutate(
    cell_type = sub("TF_(longcovid|control)_(.+)\\.rds", "\\2", sample),
    group = sub("TF_(longcovid|control)_(.+)\\.rds", "\\1", sample)
  ) %>%
  group_by(TF, cell_type, group) %>%
  summarise(mean_activity = mean(activity, na.rm = TRUE), .groups = "drop")


# Save the mean_activity data frame as an RDS file
output_mean_activity <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/mean_activity.rds"
saveRDS(mean_activity, file = output_mean_activity)

# Confirm save location
message("Mean activity saved at: ", output_mean_activity)

head(mean_activity)
str(mean_activity)
summary(mean_activity)

# Show the full cell_type for the first few rows
mean_activity %>%
  slice(1:6) %>%
  select(cell_type)

# Clean the cell_type and group columns
mean_activity_cleaned <- mean_activity %>%
  mutate(
    #cell_type = sub("\\..*", "", cell_type), # Retain only the part before the first "."
    # Remove everything after the second dot (the part that contains metadata, not cell type)
    cell_type = sub("\\.[^.]*$", "", cell_type), # Keeps valid cell types with dots intact
    group = sub("\\..*", "", group)         # Retain only the part before the first "."
  )

# Check the unique cell types
unique(mean_activity_cleaned$cell_type)

head(mean_activity_cleaned)

# Define the output path for cleaned data
output_cleaned_mean_activity <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/cleaned_mean_activity.rds"

# Save the cleaned mean_activity data frame as an RDS file
saveRDS(mean_activity_cleaned, file = output_cleaned_mean_activity)

# Confirm save location
message("Cleaned mean activity saved at: ", output_cleaned_mean_activity)

# Inspect the cleaned data
head(mean_activity_cleaned)
str(mean_activity_cleaned)
summary(mean_activity_cleaned)


## top 20 TFs heatmap across all cell type; the top 20 TFs globally and then filters data for Long Covid and Contro; identify globally variable TFs and compare how their activity differs between Long Covid and Control
# Load required libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(tibble)


# Define the file path for the cleaned mean activity data
cleaned_mean_activity_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/cleaned_mean_activity.rds"

# Load the cleaned mean activity data
mean_activity_cleaned <- readRDS(cleaned_mean_activity_file)

# Confirm the structure of the cleaned data
head(mean_activity_cleaned)

# Check for TFs with zero or NA activity in any cell type
unexpressed_tfs <- mean_activity_cleaned %>%
  group_by(TF, cell_type) %>%
  summarise(mean_activity = mean(mean_activity, na.rm = TRUE), .groups = "drop") %>%
  filter(mean_activity == 0 | is.na(mean_activity)) %>%
  pull(TF)

# List unique TFs not expressed in some cell types
unique(unexpressed_tfs)

# Select the top 20 most variable TFs
top_tfs <- mean_activity_cleaned %>%
  group_by(TF) %>%
  summarise(sd_activity = sd(mean_activity, na.rm = TRUE)) %>%
  arrange(desc(sd_activity)) %>%
  slice_head(n = 20) %>%
  pull(TF)

# Define heatmap colors
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

# Define the output file paths for the heatmaps
longcovid_heatmap_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/top20_TF_activities_heatmap_longcovid.png"
control_heatmap_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/top20_TF_activities_heatmap_control.png"

# Summarize duplicates before reshaping
heatmap_data_longcovid <- mean_activity_cleaned %>%
  filter(TF %in% top_tfs, group == "longcovid") %>%
  group_by(cell_type, TF) %>%  # Group by cell type and TF to handle duplicates
  summarise(mean_activity = mean(mean_activity, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = TF, values_from = mean_activity) %>%
  column_to_rownames("cell_type") %>%
  as.matrix()


# Plot the heatmap for LongCovid
pheatmap::pheatmap(
  mat = heatmap_data_longcovid,
  color = colors.use,
  border_color = "white",
  cellwidth = 15,
  cellheight = 15,
  treeheight_row = 20,
  treeheight_col = 20,
  main = "Top 20 Most Variable TF Activities (LongCovid)",
  filename = longcovid_heatmap_file, # Save the LongCovid plot
  width = 10,
  height = 8
)

# Confirm save location for LongCovid
message("LongCovid heatmap saved at: ", longcovid_heatmap_file)

# Repeat for Control group
heatmap_data_control <- mean_activity_cleaned %>%
  filter(TF %in% top_tfs, group == "control") %>%
  group_by(cell_type, TF) %>%  # Group by cell type and TF to handle duplicates
  summarise(mean_activity = mean(mean_activity, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = TF, values_from = mean_activity) %>%
  column_to_rownames("cell_type") %>%
  as.matrix()

# Plot the heatmap for Control
pheatmap::pheatmap(
  mat = heatmap_data_control,
  color = colors.use,
  border_color = "white",
  cellwidth = 15,
  cellheight = 15,
  treeheight_row = 20,
  treeheight_col = 20,
  main = "Top 20 Most Variable TF Activities (Control)",
  filename = control_heatmap_file, # Save the Control plot
  width = 10,
  height = 8
)

# Confirm save location for Control
message("Control heatmap saved at: ", control_heatmap_file)


##  top 20 TFs heatmap across all cell type; top 20 TFs separately for longcovid and control; to highlight TFs that are most variable within each group (Long Covid and Control)
# Load required libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(tibble)

# Define the file path for the cleaned mean activity data
cleaned_mean_activity_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/cleaned_mean_activity.rds"

# Load the cleaned mean activity data
mean_activity_cleaned <- readRDS(cleaned_mean_activity_file)

# Confirm the structure of the cleaned data
head(mean_activity_cleaned)

# Define heatmap colors
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

# Define the output file paths for the heatmaps
longcovid_heatmap_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/top20_TF_activities_heatmap_longcovid_02.png"
control_heatmap_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/top20_TF_activities_heatmap_control_02.png"

# Step 1: Select the top 20 most variable TFs for Long Covid and Control
top_tfs_longcovid <- mean_activity_cleaned %>%
  filter(group == "longcovid") %>%
  group_by(TF) %>%
  summarise(sd_activity = sd(mean_activity, na.rm = TRUE)) %>%
  arrange(desc(sd_activity)) %>%
  slice_head(n = 20) %>%
  pull(TF)

top_tfs_control <- mean_activity_cleaned %>%
  filter(group == "control") %>%
  group_by(TF) %>%
  summarise(sd_activity = sd(mean_activity, na.rm = TRUE)) %>%
  arrange(desc(sd_activity)) %>%
  slice_head(n = 20) %>%
  pull(TF)

# Step 2: Prepare heatmap data for Long Covid
heatmap_data_longcovid <- mean_activity_cleaned %>%
  filter(TF %in% top_tfs_longcovid, group == "longcovid") %>%
  group_by(cell_type, TF) %>%
  summarise(mean_activity = mean(mean_activity, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = TF, values_from = mean_activity) %>%
  column_to_rownames("cell_type") %>%
  as.matrix()

# Plot and save the heatmap for Long Covid
pheatmap::pheatmap(
  mat = heatmap_data_longcovid,
  color = colors.use,
  border_color = "white",
  cellwidth = 15,
  cellheight = 15,
  treeheight_row = 20,
  treeheight_col = 20,
  main = "Top 20 Most Variable TF Activities (LongCovid)",
  filename = longcovid_heatmap_file, # Save the LongCovid plot
  width = 10,
  height = 8
)

# Confirm save location for LongCovid
message("LongCovid heatmap saved at: ", longcovid_heatmap_file)

# Step 3: Prepare heatmap data for Control
heatmap_data_control <- mean_activity_cleaned %>%
  filter(TF %in% top_tfs_control, group == "control") %>%
  group_by(cell_type, TF) %>%
  summarise(mean_activity = mean(mean_activity, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = TF, values_from = mean_activity) %>%
  column_to_rownames("cell_type") %>%
  as.matrix()

# Plot and save the heatmap for Control
pheatmap::pheatmap(
  mat = heatmap_data_control,
  color = colors.use,
  border_color = "white",
  cellwidth = 15,
  cellheight = 15,
  treeheight_row = 20,
  treeheight_col = 20,
  main = "Top 20 Most Variable TF Activities (Control)",
  filename = control_heatmap_file, # Save the Control plot
  width = 10,
  height = 8
)

# Confirm save location for Control
message("Control heatmap saved at: ", control_heatmap_file)



## ## top20_TF_activities_difference; identify globally variable TFs and compare how their activity differs between Long Covid and Control
# Load required libraries
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(tibble)

# Define the file path for the cleaned mean activity data
cleaned_mean_activity_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/cleaned_mean_activity.rds"

# Load the cleaned mean activity data
mean_activity_cleaned <- readRDS(cleaned_mean_activity_file)

# Confirm the structure of the cleaned data
head(mean_activity_cleaned)

# Define heatmap colors
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

# Calculate mean activity per cell type and TF for both groups
mean_activity_diff <- mean_activity_cleaned %>%
  group_by(cell_type, TF, group) %>%
  summarise(mean_activity = mean(mean_activity, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean_activity) %>%
  mutate(difference = longcovid - control) %>%
  select(TF, cell_type, difference)

# Select top 20 most variable TFs based on the difference
top_tfs_diff <- mean_activity_diff %>%
  group_by(TF) %>%
  summarise(sd_diff = sd(difference, na.rm = TRUE)) %>%
  arrange(desc(sd_diff)) %>%
  slice_head(n = 20) %>%
  pull(TF)

# Prepare heatmap data for differences
heatmap_data_diff <- mean_activity_diff %>%
  filter(TF %in% top_tfs_diff) %>%
  group_by(cell_type, TF) %>%  # Handle potential duplicates
  summarise(difference = mean(difference, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = TF, values_from = difference) %>%
  column_to_rownames("cell_type") %>%
  as.matrix()

# Define output file for difference heatmap
heatmap_output_file_diff <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/top20_TF_activities_difference_heatmap.png"

# Plot and save the difference heatmap
pheatmap::pheatmap(
  mat = heatmap_data_diff,
  color = colors.use,
  border_color = "white",
  cellwidth = 15,
  cellheight = 15,
  treeheight_row = 20,
  treeheight_col = 20,
  main = "Top 20 Most Variable TF Activity Differences (LongCovid - Control)",
  filename = heatmap_output_file_diff,
  width = 10,
  height = 8
)

# Confirm save location
message("Difference heatmap saved at: ", heatmap_output_file_diff)



## annotate top 20 TF with target genes in Longcovid and control 
# Load required libraries
library(Seurat)
library(dplyr)
library(decoupleR)
library(pheatmap)
library(RColorBrewer)
library(tibble)
library(openxlsx)

# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Load the cleaned mean activity data
cleaned_mean_activity_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/cleaned_mean_activity.rds"
mean_activity_cleaned <- readRDS(cleaned_mean_activity_file)

head(mean_activity_cleaned)


# Step 1: Extract genes from the Seurat object
dataset_genes <- rownames(seurat_object@assays$RNA@data)

# Extract unique TFs from mean_activity_cleaned
unique_tfs <- unique(mean_activity_cleaned$TF)

# Step 2: Select the top 20 most variable TFs for Long Covid and Control
top_tfs_longcovid <- mean_activity_cleaned %>%
  filter(group == "longcovid") %>%
  group_by(TF) %>%
  summarise(sd_activity = sd(mean_activity, na.rm = TRUE)) %>%
  filter(TF %in% unique_tfs) %>%  # Ensure TFs are from mean_activity_cleaned
  arrange(desc(sd_activity)) %>%
  slice_head(n = 20) %>%
  pull(TF)

top_tfs_control <- mean_activity_cleaned %>%
  filter(group == "control") %>%
  group_by(TF) %>%
  summarise(sd_activity = sd(mean_activity, na.rm = TRUE)) %>%
  filter(TF %in% unique_tfs) %>%  # Ensure TFs are from mean_activity_cleaned
  arrange(desc(sd_activity)) %>%
  slice_head(n = 20) %>%
  pull(TF)

# Step 3: Load CollecTRI regulons (human)
net <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)

# Step 4: Filter net for Long Covid and Control TFs, and ensure targets are in dataset_genes
filtered_net_longcovid <- net %>%
  filter(source %in% top_tfs_longcovid & target %in% dataset_genes)

filtered_net_control <- net %>%
  filter(source %in% top_tfs_control & target %in% dataset_genes)

# Step 5: Annotate target genes for Long Covid and Control
annotated_targets_longcovid <- filtered_net_longcovid %>%
  select(source, target)

annotated_targets_control <- filtered_net_control %>%
  select(source, target)

# Step 6: Organize target genes by TF
organized_targets_longcovid <- annotated_targets_longcovid %>%
  group_by(source) %>%
  summarise(target_genes = paste(target, collapse = ", "), .groups = "drop")

organized_targets_control <- annotated_targets_control %>%
  group_by(source) %>%
  summarise(target_genes = paste(target, collapse = ", "), .groups = "drop")

# Step 7: Save to Excel
output_excel_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/organized_top20_TF_targets_filtered.xlsx"

# Create an Excel workbook
wb <- createWorkbook()

# Add worksheets for Long Covid and Control
addWorksheet(wb, "LongCovid Top 20 TFs")
writeData(wb, "LongCovid Top 20 TFs", organized_targets_longcovid)

addWorksheet(wb, "Control Top 20 TFs")
writeData(wb, "Control Top 20 TFs", organized_targets_control)

# Save the workbook
saveWorkbook(wb, file = output_excel_file, overwrite = TRUE)

message("Filtered and organized target genes saved at: ", output_excel_file)



## annotate top 20 TF difference between Longcovid and control with target genes  
# Load required libraries
library(Seurat)
library(dplyr)
library(decoupleR)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(tibble)
library(openxlsx)

# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Load the cleaned mean activity data
cleaned_mean_activity_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/cleaned_mean_activity.rds"
mean_activity_cleaned <- readRDS(cleaned_mean_activity_file)

head(mean_activity_cleaned)

# Step 1: Extract genes from the Seurat object
dataset_genes <- rownames(seurat_object@assays$RNA@data)

# Step 2: Calculate mean activity differences between Long Covid and Control
mean_activity_diff <- mean_activity_cleaned %>%
  group_by(cell_type, TF, group) %>%
  summarise(mean_activity = mean(mean_activity, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean_activity) %>%
  mutate(difference = longcovid - control) %>%
  select(TF, cell_type, difference)

# Step 3: Select the top 20 most variable TFs based on the differences
top_tfs_diff <- mean_activity_diff %>%
  group_by(TF) %>%
  summarise(sd_diff = sd(difference, na.rm = TRUE)) %>%
  arrange(desc(sd_diff)) %>%
  slice_head(n = 20) %>%
  pull(TF)

# Step 4: Load CollecTRI regulons (human)
net <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)

# Step 5: Filter net for the selected top 20 TFs and ensure targets are in the dataset genes
filtered_net_diff <- net %>%
  filter(source %in% top_tfs_diff & target %in% dataset_genes)

# Step 6: Annotate target genes for the selected top 20 TFs
annotated_targets_diff <- filtered_net_diff %>%
  select(source, target)

# Step 7: Organize target genes by TF
organized_targets_diff <- annotated_targets_diff %>%
  group_by(source) %>%
  summarise(target_genes = paste(target, collapse = ", "), .groups = "drop")

# Step 8: Save to Excel
output_excel_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/organized_top20_TF_targets_diff.xlsx"

# Create an Excel workbook
wb <- createWorkbook()

# Add a worksheet for the annotated TF targets
addWorksheet(wb, "Top 20 TF Differences")
writeData(wb, "Top 20 TF Differences", organized_targets_diff)

# Save the workbook
saveWorkbook(wb, file = output_excel_file, overwrite = TRUE)

message("Filtered and organized target genes based on TF differences saved at: ", output_excel_file)





library(Seurat)
library(dplyr)
library(decoupleR)
library(pheatmap)
library(RColorBrewer)
library(tibble)

# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Load an example `ulm_results` file
ulm_results <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/TF_control_B_exhausted.rds")

head(ulm_results)
str(ulm_results)

# Load CollecTRI regulons (human)
net <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)

# Filter `net` to include only TFs present in `ulm_results`
filtered_net <- net %>% filter(source %in% unique(ulm_results$source))

# Annotate `ulm_results` with target genes
annotated_results <- ulm_results %>%
  left_join(filtered_net, by = c("source" = "source"))

# Check the structure
head(annotated_results)

# Step 1: Extract genes from the Seurat object
dataset_genes <- rownames(seurat_object@assays$RNA@data)

# Step 2: Filter targets based on genes in the dataset
filtered_annotated_results <- annotated_results %>%
  filter(target %in% dataset_genes)

# Step 3: Summarize target gene activities
target_gene_scores <- filtered_annotated_results %>%
  mutate(condition = paste(condition)) %>%  # Use the condition column
  group_by(target, condition) %>%
  summarise(mean_activity = mean(score, na.rm = TRUE), .groups = "drop")

# Step 4: Pivot to wide format for visualization
target_gene_matrix <- target_gene_scores %>%
  pivot_wider(names_from = condition, values_from = mean_activity) %>%
  column_to_rownames("target")

# Check the resulting matrix
print(dim(target_gene_matrix))
target_gene_matrix[1:3, 1:3]




