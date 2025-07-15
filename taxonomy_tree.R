# ╔═══════════════════════════════════════════════════════════════════╗
# ║                      taxonomy_tree.R                              ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : antarctica                                       ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-07-10                                       ║
# ║ Last Modified  : 2025-07-15                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/antarctica        ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(metacoder)
library(qiime2R)


## Prepare data

set.seed(1234)
local_metadata <- "antartida-16S"
amplicon <- "16S" # ITS or 16S

metadata_file_path <- file.path("/home/sergio/scratch",
                                local_metadata,
                                "metadata.tsv")

sample <- read.csv(metadata_file_path,
                   header = TRUE,
                   sep = "\t")
sample$dummy_col <- rep("yes", nrow(sample))

setwd(file.path("~/projects/antarctica", amplicon))

feature_table <- qiime2R::read_qza("filtered_table.qza")$data

taxonomy <- qiime2R::read_qza("taxonomy.qza")$data # |> qiime2R::parse_taxonomy()

asv <- data.frame(
  ASV = rownames(feature_table),
  lineage = sub(";sh__.*$", "", taxonomy$Taxon[match(rownames(feature_table), taxonomy$Feature.ID)]),
  feature_table,
  row.names = NULL
)

obj <- parse_tax_data(asv,
                      class_cols = "lineage",
                      class_sep = ";",
                      class_regex = "^(.+)__(.+)$",
                      class_key = c(tax_rank = "info",
                                    tax_name = "taxon_name"))

obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data",
                                       cols = sample$ID)

obj$data$tax_occ <- calc_n_samples(obj, "tax_abund",
                                   groups = sample$dummy_col,
                                   cols = sample$ID)

obj$data$diff_table <- compare_groups(obj,
                                      dataset = "tax_abund",
                                      cols = sample$ID,
                                      groups = sample$Depth)
print(obj$data$diff_table)

## Plots

### General heat tree

heat_tree(obj, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = yes, 
          node_size_axis_label = "ASV count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel",
          initial_layout = "reingold-tilford")

### Pairwise comparison

heat_tree(obj, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = log2_median_ratio,
          node_color_interval = c(-2, 2),
          node_color_range = c("cyan", "gray", "tan"),
          node_size_axis_label = "ASV count",
          node_color_axis_label = "Log 2 ratio of median proportions",
          layout = "davidson-harel",
          initial_layout = "reingold-tilford")

obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value,
                                               method = "fdr")
range(obj$data$diff_table$wilcox_p_value, finite = TRUE) 

### All-vs-all comparison
