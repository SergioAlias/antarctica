# ╔═══════════════════════════════════════════════════════════════════╗
# ║                      taxonomy_tree.R                              ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : antarctica                                       ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-07-10                                       ║
# ║ Last Modified  : 2025-07-10                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/antarctica        ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(metacoder)
library(qiime2R)


## Prepare data

setwd("projects/antarctica/")

sample <- data.frame(
  sample_id = c("AN1_1__ITS", "AN10_1__ITS", "AN10_2__ITS"),
  sample_name = c("AN1__ITS", "AN10__ITS", "AN10__ITS"),
  depth = c("338", "337", "337") 
)

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
                                       cols = sample$sample_id)

obj$data$diff_table <- compare_groups(obj,
                                      dataset = "tax_abund",
                                      cols = sample$sample_id,
                                      groups = sample$depth)
print(obj$data$diff_table)

## Plot

set.seed(1234)

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
