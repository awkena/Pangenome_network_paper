
# Install panGenomeBreedr from GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("awkena/panGenomeBreedr")

# Import raw KASP genotyping file (.csv) using the read_kasp_csv() function
library(panGenomeBreedr)

# Set path to the directory where your data is located
# path1 <-  "inst/extdata/Genotyping_141.010_01.csv"
path1 <-  system.file("extdata", "Genotyping_141.010_01.csv",
                      package = "panGenomeBreedr",
                      mustWork = TRUE)

# Import raw data file
file1 <- read_kasp_csv(file = path1,
                       row_tags = c("Statistics", "DNA", "SNPs", "Scaling", "Data"),
                       data_type = 'raw')

# Get KASP genotyping data for plotting
kasp_dat <- file1$Data


# Assign KASP fluorescence colors using the kasp_color() function
# Create a subet variable called plates: masterplate x snpid
kasp_dat$plates <- paste0(kasp_dat$MasterPlate, '_',
                          kasp_dat$SNPID)

dat1 <- kasp_color(x = kasp_dat,
                   subset = 'plates',
                   sep = ':',
                   geno_call = 'Call',
                   uncallable = 'Uncallable',
                   unused = '?',
                   blank = 'NTC',
                   assign_cols = c(FAM = "cornflowerblue", HEX = "firebrick2" ,
                                   het = "forestgreen"))


# KASP QC plot for Plate 5
kasp_qc_ggplot2(x = dat1[5],
                pdf = TRUE,
                filename = 'lgs1_paper_cluster_plot2',
                Group_id = NULL,
                scale = TRUE,
                expand_axis = 0.6,
                alpha = 0.9,
                legend.pos.x = 0.6,
                legend.pos.y = 0.8)


# KASP QC plot for Plate 5
kasp_qc_ggplot2(x = dat1[5],
                pdf = TRUE,
                Group_id = 'Group',
                filename = 'lgs1_paper_cluster_plot_with_control',
                Group_unknown = '?',
                scale = TRUE,
                pred_cols = c('Blank' = 'black', 'False' = 'firebrick3',
                              'True' = 'blue', 'Unverified' = 'beige'),
                expand_axis = 0.6,
                alpha = 0.9,
                legend.pos.x = 0.6,
                legend.pos.y = 0.8)

# Get prediction summary for all plates
my_sum <- pred_summary(x = dat1,
                       snp_id = 'SNPID',
                       Group_id = 'Group',
                       Group_unknown = '?',
                       geno_call = 'Call',
                       rate_out = TRUE)

View(my_sum$summ)

# Get prediction summary for snp:snpSB00804
my_sum <- my_sum$summ
my_sum <- my_sum[my_sum$snp_id == 'snpSB00804',]

pred_summary_plot(x = my_sum,
                  pdf = TRUE,
                  filename = 'prediction_summary_plot',
                  pred_cols = c('false' = 'firebrick3', 'true' = 'blue',
                                'unverified' = 'beige'),
                  alpha = 1,
                  text_size = 12,
                  width = 6,
                  height = 6,
                  angle = 45)
