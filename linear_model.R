# ╔═══════════════════════════════════════════════════════════════════╗
# ║                         linear_model.R                            ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : antarctica                                       ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-07-22                                       ║
# ║ Last Modified  : 2025-07-22                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/antarctica        ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(magrittr, include.only = "%<>%")
library(tidyverse)
library(ggpmisc)
library(qiime2R)


## Import QIIME 2 files

project_name <- "antartida_16S"
local_metadata <- "antartida-16S"
out <- "antartida-16S"
remove_blanks <- TRUE
alpha_mode <- FALSE
asv_id <- "9c8fe4f590ed08b997a8ca49bb3bb647"

readRenviron("/home/sergio/Renvs/.RenvBrigit")
brigit_IP <- Sys.getenv("IP_ADDRESS")
cluster_path <- paste0("/run/user/1001/gvfs/sftp:host=",
                       brigit_IP,
                       ",user=salias/mnt/lustre")
project_dir <- file.path(cluster_path,
                         "scratch/salias/projects",
                         project_name)
outdir <- file.path("/home/sergio/scratch",
                    out,
                    "lm")

metadata <- read.csv(file.path("/home/sergio/scratch",
                               local_metadata,
                               "metadata.tsv"),
                     sep = "\t")

shannon_file_path <- file.path(project_dir,
                               "qiime2/diversity/shannon_vector.qza")
dada2_file_path <- file.path(project_dir,
                             "qiime2/feature_tables/filtered_table.qza")

shannon <- read_qza(shannon_file_path)
shannon <- shannon$data %>% rownames_to_column("SampleID")
dada2 <- read_qza(dada2_file_path)$data
dada2 %<>% apply(2, function(x) x / sum(x))
asv <- dada2[asv_id,, drop = FALSE] %>% t() %>% as.data.frame()
asv$SampleID <- rownames(asv)
colnames(asv) <- c("asv", "SampleID")
rownames(asv) <- NULL

## Merge by sample ID

metadata %<>% rename(SampleID = ID) %>% left_join(shannon) %>% left_join(asv)

## Remove blanks (optional)

if (remove_blanks == TRUE){
  metadata %<>% filter(Depth != 0.0)
}


## Scatterplot with linear fit

if (alpha_mode == TRUE){
  metadata$response <- metadata$shannon_entropy
  plot_title <- "Shannon Diversity vs. Depth"
  plot_y_axis <- "Shannon Index"
} else {
  metadata$response <- metadata$asv
  plot_title <- paste0(asv_id, " vs. Depth")
  plot_y_axis <- paste0(substr(asv_id, 1, 7), " relative abundance (%)")
}

ggplot(metadata, aes(x = Depth, y = response)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue", formula = y ~ x) +
  stat_poly_eq(
    formula = y ~ x,
    aes(label = paste(..eq.label.., ..p.value.label.., sep = "~~~")),
    parse = TRUE,
    label.x = "right",
    label.y = "top"
  ) +
  theme_minimal() +
  labs(
    title = plot_title,
    x = "Depth (m)",
    y = plot_y_axis
  )


## Fit linear model

model <- lm(shannon_entropy ~ Depth, data = metadata)
summary(model)
