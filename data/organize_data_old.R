library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

count_cluster <- read_csv("RAW/counts_clusters.csv")
count_genera <- read_csv("RAW/OLD/genera_counts.csv")
metadata <- read_csv("RAW/metadata_fixed.csv") %>% select(-1) 


count_cluster <- count_cluster %>%
  mutate(sample_name = if_else(sample_name == "DTU_2017_314_1_MG_CI_YA", 
                               "DTU_2017_1011594_1_MG_CI_YA", 
                               sample_name)) %>%
  mutate(sample_id = str_split_fixed(sample_name, "_", 4)[,3],
         replicate_id = str_split_fixed(sample_name, "_", 5)[,4],
         .before = 1) %>%
  filter(sample_id != "")

count_genera <- count_genera %>%
  mutate(run_accession = if_else(run_accession == "DTU_2017_314_1_MG_CI_YA", 
                                 "DTU_2017_1011594_1_MG_CI_YA", 
                                 run_accession)) %>%
  mutate(sample_id = str_split_fixed(run_accession, "_", 4)[,3],
         replicate_id = str_split_fixed(run_accession, "_", 5)[,4],
         .before = 1) %>%
  filter(sample_id != "")

metadata$bindCol[metadata$bindCol == "DTU_2017_314_1_MG_CI_YA"] <- "DTU_2017_1011594_1_MG_CI_YA"
metadata$genepid[metadata$bindCol == "DTU_2017_1011594_1_MG_CI_YA"] <- "1011594"
all(count_genera$sample_id %in% metadata$genepid)
all(count_cluster$sample_id %in% metadata$genepid)


count_genera_noRep_max <- count_genera %>%
  mutate(sample_counts = rowSums(select(., Abiotrophia:last_col())),
         .after = run_accession) %>%
  group_by(sample_id) %>%
  arrange(desc(sample_counts), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  filter(run_accession != "DTU_2016_1012307_1_MG_USA_74c") %>%
  unite("sample_replicate_id", sample_id:replicate_id, remove = FALSE) %>%
  relocate(1)

# count_genera_noRep_sum <-  count_genera %>%
#   select(-run_accession) %>%
#   pivot_longer(
#     cols = Abiotrophia:last_col(),
#     names_to = "genera",
#     values_to = "counts"
#   ) %>%
#   reframe(counts = sum(counts), .by = c(sample_id, genera)) %>%
#   pivot_wider(names_from = genera, values_from = counts, values_fill = 0)

# AMR !!!
#------------------------------------------------------------------------------#
count_cluster_noRep_max <- count_cluster %>%
  unite("sample_replicate_id", sample_id:replicate_id, remove = FALSE) %>%
  filter(sample_replicate_id %in% count_genera_noRep_max$sample_replicate_id) %>%
  mutate(sample_depth = sum(fragmentCountAln_adj), 
         .by = c(sample_id, replicate_id), .after = sample_name) %>%
  group_by(sample_id) %>%
  filter(sample_depth == max(sample_depth)) %>%
  ungroup() %>%
  filter(sample_name != "DTU_2016_1012307_1_MG_USA_74c")

abun_genera <- count_genera_noRep_max %>%
  select(-sample_id, -replicate_id, -run_accession, -sample_counts) %>%
  column_to_rownames("sample_replicate_id") %>%
  as.matrix()

abun_cluster <- count_cluster_noRep_max %>%
  select(sample_replicate_id, cluster_representative_98, fragmentCountAln_adj) %>%
  pivot_wider(names_from = cluster_representative_98, values_from = fragmentCountAln_adj, 
              values_fill = 0) %>%
  column_to_rownames("sample_replicate_id") %>%
  as.matrix()

all(sort(rownames(abun_cluster)) == sort(rownames(abun_genera)))
abun_cluster <- abun_cluster[rownames(abun_genera),]

info_sample <- metadata %>%
  select(bindCol, genepid, country_alt, city, year, lat, lon, country_alpha3, Region) %>%
  mutate(sample_id = str_split_fixed(bindCol, "_", 4)[,3],
         replicate_id = str_split_fixed(bindCol, "_", 5)[,4],
         .before = 1) %>%
  unite("sample_replicate_id", sample_id:replicate_id, remove = FALSE) %>%
  filter(sample_replicate_id %in% rownames(abun_genera)) %>%
  filter(bindCol != "DTU_2016_1012307_1_MG_USA_74c") %>%
  select(-bindCol) %>%
  distinct() %>%
  group_by(sample_replicate_id) %>%
  arrange(year, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  column_to_rownames("sample_replicate_id")
info_sample <- info_sample[rownames(abun_genera), ]


info_taxa_AMR <- count_cluster %>% select(cluster_representative_98:group) %>%
  distinct() %>%
  column_to_rownames("cluster_representative_98")

all(colnames(abun_cluster) %in% rownames(info_taxa_AMR))
info_taxa_AMR <- info_taxa_AMR[colnames(abun_cluster), ]


## LAST MINOR FIXES
info_sample <- info_sample %>%
  rename("sample_code" = "sample_id")

if(all(info_sample$sample_code == info_sample$genepid)){
  info_sample <- info_sample %>%
    select(-genepid)
}

## PanRes



library(mgnet)
global_bact_genera <- mgnet(abundance = abun_genera,
                            info_sample = info_sample)

global_AMR <- mgnet(abundance = abun_cluster,
                    info_sample = info_sample,
                    info_taxa = info_taxa_AMR)

saveRDS(global_bact_genera, "R/global_sewage_bacteria_genera.rds")
saveRDS(global_AMR, "R/global_sewage_AMR_cluster98.rds")
