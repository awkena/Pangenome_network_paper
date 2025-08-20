
# # Set working directory to save figures
# setwd("~/Desktop/lgs1_figures")
# list.files()

# Please install panGB from GitHub to get updated scripts
# Install panGenomeBreedr
if (!require("devtools")) install.packages("devtools")
devtools::install_github("awkena/panGenomeBreedr")
library(panGenomeBreedr)

# Load Agriplex data for Mota Maradi NILs from GitHub repo
path1 <- "https://raw.githubusercontent.com/awkena/Pangenome_network_paper/main/mota_maradi_agriplex_marker_data.RData"

download.file(path1, destfile = file.path(tempdir(), basename(path1)))

load(file.path(tempdir(), basename(path1)))

# Subset data for lgs1 NILs and assign row names
lgs1 <- geno[, -c(1:8)]
rownames(lgs1) <- geno$Geno_Ord[1:10]  

# Before generating heatmaps, we must clean the data to remove the following:
# 1. Markers with at least one missing parent genotype call
# 2. Markers where the parents are monomorphic
# 3. Markers where at least one parent is heterozygous

# There are custom functions in panGB for all of these data cleaning steps

# Find loci with at least one missing parent genotype
par_miss <- parent_missing(x = lgs1,
                           rp_row = 1,
                           dp_row = 2,
                           na_code = "FAIL")

View(par_miss$par_missing)

# Extract loci with all parents genotype calls present
lgs1 <- par_miss$par_present

# Find polymorphic loci between MM and SRN39; better than using the rm_mono()
lgs1 <- parent_poly(x = lgs1,
                        rp_row = 1,
                        dp_row = 2,
                        sep = ' / ')

# Find loci with at least one heterozygous parent genotype
par_het <- parent_het(x = lgs1,
                      rp_row = 1,
                      dp_row = 2,
                      sep = ' / ')

View(par_het$par_het)

# Extract markers with parents being homozygous
lgs1 <- par_het$par_hom


# Find loci with potential genotype call errors:
# Check for genotype call error for each marker
geno_mat <- geno_error(x = lgs1,
                       rp_row = 1,
                       dp_row = 2,
                       sep = ' / ',
                       data_type = 'agriplex')

# No genotyping error
View(geno_mat$geno_err)

# Parse snp ids to generate a map file
snps <- colnames(lgs1) # Get snp ids
map_file <- parse_marker_ns(x = snps, sep = '_', prefix = 'S')

# order markers in map file
map_file <- order_markers(x = map_file)

# Process genotype data to re-order SNPs based on chromosome and positions
lgs1 <- proc_kasp(x = lgs1,
                  kasp_map = map_file,
                  map_snp_id = "snpid",
                  sample_id = "GENOTYPE",
                  marker_start = 1,
                  chr = 'chr',
                  chr_pos = 'pos')

# Convert to numeric format for plotting
num_geno <- kasp_numeric(x = lgs1,
                         rp_row = 1,
                         dp_row = 2,
                         sep = ' / ',
                         data_type = 'agriplex')

# Check for unique numeric codes in num_geno
unique(as.vector(num_geno))


# # Find loci with missing data
# cols_miss <- colnames(num_geno)[sapply(as.data.frame(num_geno), function(x) any(x == -5.0, na.rm = TRUE))]
# lgs1[,cols_miss]
# num_geno[,cols_miss]

# Create a heatmap that compares the parents to progenies
# Using the color codes given by Carl
# Annotate agronomic loci with vertical lines
trait_pos <- list(Stg5 = c(chr = 1, pos = 1e6),
                 Stg1 = c(chr = 3, pos = 66e6),
                 LGS1 = c(chr = 5, pos = 70e6),
                 Tan1 = c(chr = 4, pos = 62e6),
                 Tan2 = c(chr = 2, pos = 7.9e6),
                 Stg3a = c(chr = 2, pos = 56e6),
                 Stg3b = c(chr = 2, pos = 70.5e6),
                 Stg4 = c(chr = 5, pos = 13e6),
                 Dw1 = c(chr = 9, pos = 57e6),
                 Dw2 = c(chr = 6, pos = 43.8e6),
                 Dw3 = c(chr = 7, pos = 59e6),
                 Ma1 = c(chr = 6, pos = 40e6),
                 Ma2 = c(chr = 2, pos = 66e6),
                 Ma3 = c(chr = 1, pos = 68e6),
                 Ma6 = c(chr = 6, pos = 1e6),
                 Bmr6 = c(chr = 4, pos = 5.8e6),
                 Bmr12 = c(chr = 7, pos = 4.7e6))

cross_qc_heatmap2(x = num_geno,
                  map_file = map_file,
                  snp_ids = 'snpid',
                  chr = 'chr',
                  chr_pos = 'pos',
                  parents = c("Mota Maradi (S)", "SRN39 (R)"),
                  col_mapping=c('-1' = "grey",  '-2'= "pink", '0' = "#F54646",
                                '0.5' = "yellow", '1' = "#5249f5", '-5' = "purple"),
                  pdf = FALSE,
                  trait_pos = trait_pos,
                  filename = 'background_heatmap',
                  text_scale_fct = 0.4,
                  legend_title = '',
                  alpha = 0.9,
                  text_size = 10)

###########################################################################
# # Zoom in to annotate lgs1 position on Chr 5
# lgs1_ch5 <- num_geno[, map_file$chr == 5 & map_file$pos >= 50e6
#                      & map_file$pos <= 71726595]
# rownames(lgs1_ch5)
# 
# # Get the map file for subset data
# lgs1_ch5_map <- parse_marker_ns(colnames(lgs1_ch5))
# 
# 
# cross_qc_annotate(x = lgs1_ch5,
#                   map_file = lgs1_ch5_map,
#                   snp_ids = 'snpid',
#                   chr = 'chr',
#                   chr_pos = 'pos',
#                   parents = c("Mota Maradi (S)", "SRN39 (R)"),
#                   col_mapping=c('-1' = "grey",  '-2'= "pink", '0' = "#F54646",
#                                 '0.5' = "yellow", '1' = "#5249f5", '-5' = "purple"),
#                   trait_pos = list(lgs1 = c(chr = 5, start = 69.8e6, end = 70.02e6)),
#                   text_scale_fct = 0.3,
#                   pdf = FALSE,
#                   legend_title = 'Lgs1_NILs',
#                   alpha = 0.9,
#                   text_size = 12)
# 
# lgs1_R <- lgs1_ch5[c("Mota Maradi (S)",
#                      "SRN39 (R)",
#                      "Line 17 (R-MM)",
#                      "Line 19 (R-MM)",
#                      "Line 13 (R-MM)",
#                      "Line 15 (R-MM)"),]
# 
# cross_qc_annotate(x = lgs1_R,
#                   map_file = lgs1_ch5_map,
#                   snp_ids = 'snpid',
#                   chr = 'chr',
#                   chr_pos = 'pos',
#                   parents = c("Mota Maradi (S)", "SRN39 (R)"),
#                   col_mapping=c('-1' = "grey",  '-2'= "pink", '0' = "#F54646",
#                                 '0.5' = "yellow", '1' = "#5249f5", '-5' = "purple"),
#                   trait_pos = list(lgs1 = c(chr = 5, start = 69.8e6, end = 70.02e6)),
#                   text_scale_fct = 0.3,
#                   pdf = FALSE,
#                   legend_title = 'Lgs1_R_NILs',
#                   alpha = 0.9,
#                   text_size = 12)
# 
# write.csv(lgs1_R, file = 'lgs1_linkage_drag_in_R_NILs.csv')

# zoomed-in plots for all 16 genes
# Define start and end positions for zoomed-in plots for all 16 genes
trait_pos <- list(Stg5 = c(chr = 1, start = 0.5e6, end = 6.1e6),
                  Stg1 = c(chr = 3, start = 65.5e6, end = 66.5e6),
                  LGS1 = c(chr = 5, start = 68.8e6, end = 71.02e6),
                  Tan1 = c(chr = 4, start = 60.3e6, end = 65.5e6),
                  Tan2 = c(chr = 2, start = 5.1e6, end = 10.5e6),
                  Stg3a = c(chr = 2, start = 55e6, end = 60e6),
                  Stg3b = c(chr = 2, start = 67e6, end = 71e6),
                  Stg4 = c(chr = 5, start = 12e6, end = 14e6),
                  Dw1 = c(chr = 9, start = 56.3e6, end = 57e6),
                  Dw2 = c(chr = 6, start = 44.3e6, end = 46e6),
                  Dw3 = c(chr = 7, start = 57.7e6, end = 60e6),
                  Ma1 = c(chr = 6, start = 35e6, end = 47e6),
                  Ma2 = c(chr = 2, start = 60e6, end = 67e6),
                  Ma3 = c(chr = 1, start = 67e6, end = 69e6),
                  Ma6 = c(chr = 6, start = 0.8e6, 2.5e6),
                  Bmr6 = c(chr = 4, start = 3e6, end = 6e6),
                  Bmr12 = c(chr = 7, start = 3.5e6, end = 5.7e6))

# Create output directory
output_dir <- "zoomed_in_plots_for_MM"
if (!dir.exists(output_dir)) dir.create(output_dir)

for (gene in seq_len(length(trait_pos))) {

  # Get numeric coded data for the zoomed in region
  num_dat <- num_geno[, map_file$chr == trait_pos[[gene]][1] &
                        map_file$pos >= trait_pos[[gene]][2] &
                        map_file$pos <= trait_pos[[gene]][3]]

  # Get the map file for subset data
  gene_map <- parse_marker_ns(colnames(num_dat))

  # Generate zoomed-in plot
  cross_qc_annotate(x = num_dat,
                    map_file = gene_map,
                    snp_ids = 'snpid',
                    chr = 'chr',
                    chr_pos = 'pos',
                    parents = c("Mota Maradi (S)", "SRN39 (R)"),
                    col_mapping = c('-1' = "grey",  '-2'= "pink", '0' = "#F54646",
                                    '0.5' = "yellow", '1' = "#5249f5", '-5' = "purple"),
                    trait_pos = trait_pos[gene],
                    text_scale_fct = 0.4,
                    pdf = TRUE,
                    filename = file.path(getwd(), output_dir, paste0(paste0(names(trait_pos[gene]), '_NILs_for_MM'))),
                    legend_title = paste0(names(trait_pos[gene]), '_NILs'),
                    alpha = 0.9,
                    text_size = 12)

}


# Calculate weighted RPP
rpp <- calc_rpp_bc(x = num_geno,
                   map_file = map_file,
                   map_chr = 'chr',
                   map_pos = 'pos',
                   map_snp_ids = 'snpid',
                   rp = "Mota Maradi (S)",
                   rp_num_code = 1,
                   na_code =-5,
                   weighted = TRUE,
                   het_code = 0.5)

write.csv(rpp, file = 'recovery_percentage_MM.csv')

# Generate bar plot for RPP values
rpp_barplot(rpp_df = rpp,
            rpp_threshold = 0.93,
            text_size = 12,
            text_scale_fct = 0.2,
            alpha = 0.9,
            bar_width = 0.9,
            aspect_ratio = 0.8,
            pdf = FALSE)


