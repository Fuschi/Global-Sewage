library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

count_species <- read_csv("RAW/motus_counts/motus_agg_files/mOTU_gs3_counts.csv", name_repair = "minimal")
#count_genera <- read_csv("RAW/motus_counts/motus_agg_files/genus_gs3_counts.csv", name_repair = "minimal")
#taxa_genera_raw <- read_csv("RAW/motus_counts/motus_agg_files/genus_gs3_taxa_paths.csv", name_repair = "minimal")
taxa_species_raw <- read_csv("RAW/motus_counts/motus_agg_files/mOTU_gs3_taxa_paths.csv", name_repair = "minimal")

df_amr <- read_csv("RAW/counts_clusters.csv")
#taxa_amr <- readxl::read_excel("RAW/panres_annotations.xlsx", sheet = "wide_format")

info_sample <- read_csv("RAW/metadata_fixed.csv") %>% select(-1) 

# SAMPLE METADATA
#------------------------------------------------------------------------------#
info_sample <- info_sample %>%
  select(genepid, Region, country_alt, city, year, lat, lon) %>%
  distinct() %>%
  filter(year != 2023 & !is.na(genepid)) %>%
  column_to_rownames("genepid") 

# AMR
#------------------------------------------------------------------------------#
df_amr <- df_amr %>%
  mutate(sample_name = if_else(sample_name == "DTU_2017_314_1_MG_CI_YA", 
                               "DTU_2017_1011594_1_MG_CI_YA", 
                               sample_name)) %>%
  filter(!str_starts(sample_name, "ERR")) %>%
  mutate(genepid = str_split(sample_name, "_") %>% map_chr(3),
         replicate_id = str_split(sample_name, "_") %>% map_chr(4), .after = genepid) %>%
  unite("sample_id", genepid:replicate_id, sep = "_", remove = F) %>%
  unite("sample_long", sample_name:sample_name2, sep = "_", remove = F) %>%
  mutate(sample_sum = sum(fragmentCountAln_adj), .by = "sample_long") %>%
  ungroup() %>%
  arrange(desc(sample_sum)) %>%
  group_by(genepid, cluster_representative_98) %>%
  slice(1) %>%
  ungroup()

count_amr <- df_amr %>%
  select(genepid, cluster_representative_98, fragmentCountAln_adj) %>%
  pivot_wider(names_from = "cluster_representative_98", 
              values_from = "fragmentCountAln_adj", values_fill = 0) %>%
  column_to_rownames("genepid") %>%
  as.matrix()

info_taxa_amr <- df_amr %>%
  select(cluster_representative_98, group) %>%
  distinct() %>%
  rename(cluster_id = cluster_representative_98)
  
info_taxa_amr <- readxl::read_xlsx("RAW/panres_annotations.xlsx", sheet = "wide_format") %>%
  filter(!is.na(cluster_representative)) %>%
  select(-c(fa_name,gene_length,database)) %>%
  distinct() %>%
  separate(cluster_representative, sep = "_", remove = FALSE,
           into = c("temp1", "temp2", "cluster_version")) %>%
  unite("cluster_id", temp1, temp2, sep = "_") %>%
  group_by(cluster_id) %>%
  filter(!is.na(class)) %>%
  filter(nchar(class) == max(nchar(class))) %>%
  slice(1) %>%
  ungroup() %>%
  right_join(info_taxa_amr, by = "cluster_id") %>%
  relocate(group) %>%
  column_to_rownames("cluster_id")

  



# BACTERIA mOTU SPECIES
#------------------------------------------------------------------------------#
count_species <- count_species %>%
  column_to_rownames("genepid") %>%
  as.matrix()
count_species <- count_species[, -which(duplicated(colnames(count_species)))]

taxa_species <- taxa_species_raw %>%
  separate(path_taxs,
           into = c("kingdom_id","phylum_id","class_id","order_id","family_id","genus_id","mOTU_id"),
           remove = FALSE) %>%
  mutate(kingdom = taxizedb::taxid2name(kingdom_id),
         phylum = taxizedb::taxid2name(phylum_id),
         class = taxizedb::taxid2name(class_id),
         order = taxizedb::taxid2name(order_id),
         family = taxizedb::taxid2name(family_id),
         genus = taxizedb::taxid2name(genus_id)) %>%
  group_by(name) %>%
  slice(1) %>%
  column_to_rownames("name") %>%
  select(-c(ends_with("_id"), path))

#count_species <- count_species[, rownames(taxa_species)]
#count_species <- count_species[rownames(info_sample)]
  
# BACTERIA mOTU GENERA
#------------------------------------------------------------------------------#
# count_genera <- count_genera %>%
#   column_to_rownames("genepid") %>%
#   as.matrix()
# 
# taxa_genera <- taxa_genera_raw %>%
#   separate(path_taxs, 
#            into = c("kingdom_id","phylum_id","class_id","order_id","family_id","genus_id"),
#            remove = FALSE) %>%
#   mutate(kingdom = taxizedb::taxid2name(kingdom_id),
#          phylum = taxizedb::taxid2name(phylum_id),
#          class = taxizedb::taxid2name(class_id),
#          order = taxizedb::taxid2name(order_id),
#          family = taxizedb::taxid2name(family_id)) %>%
#   group_by(name) %>%
#   slice_max(nchar(path_taxs), with_ties = FALSE) %>%
#   column_to_rownames("name")



# SAMPLE INTERSECTION
#------------------------------------------------------------------------------#
sample_intersection <- purrr::reduce(
  .x = list(rownames(count_amr), rownames(info_sample), rownames(count_species)),
  .f = intersect
)

# SELECT THE INTERSECTION
# From 1240 to 1147 (fair enough)
#------------------------------------------------------------------------------#
count_amr <- count_amr[sample_intersection, ]
count_species <- count_species[sample_intersection, ]
#count_genera <- count_genera[sample_intersection, ]
info_sample <- info_sample[sample_intersection,]

# MERGE TAXA AND COUNT INFORMATION IN AMR
#------------------------------------------------------------------------------#
all(colnames(count_amr) %in% rownames(info_taxa_amr))
info_taxa_amr <- info_taxa_amr[colnames(count_amr), , drop = F]
info_taxa_amr <- info_taxa_amr %>% mutate(source = "amr", .before = 1) %>%
  rename(class_amr = class, group_amr = group, biocide_class_amr = biocide_class,
         metal_class_amr = metal_class) 

# MERGE TAXA AND COUNT INFORMATION IN GENERA
#------------------------------------------------------------------------------#
all(rownames(taxa_species) %in% colnames(count_species))
taxa_species <- taxa_species[colnames(count_species), ]
taxa_species <- taxa_species %>% rename(source = level)

# CREATE THE MGNET OBJECTS
#------------------------------------------------------------------------------#
library(mgnet)

gs_mOTU <- mgnet(
  abun = count_species,
  sample = info_sample,
  taxa = taxa_species
)

gs_amr <- mgnet(
  abun = as.matrix(count_amr),
  sample = info_sample,
  taxa = info_taxa_amr
)

gs <- mgnet(
  abun = cbind(as.matrix(count_species), as.matrix(count_amr)),
  sample = info_sample,
  taxa = bind_rows(taxa_species, info_taxa_amr)
)

gs@taxa <- taxa(gs) %>% mutate(type = case_when(
  source == "mOTU" ~ "mOTU",
  group_amr == "ResFinder" ~ "ResFinder",
  group_amr == "Functional" ~ "Functional",
), .after = 1) 

# SAVE DATA IN MGNET FORMATS
#------------------------------------------------------------------------------#
saveRDS(gs, "R/global-sewage.rds")

# Merge in the same mgnet
#------------------------------------------------------------------------------#
count_merge <- as_tibble(count_genera,rownames = "sample_id") %>%
  left_join(as_tibble(count_amr,rownames = "sample_id"),
            by = "sample_id") %>%
  column_to_rownames("sample_id") %>% as.matrix()

taxa_merge <- bind_rows(cbind(source = "Bacteria", taxa_genera),
                        cbind(source = "AMR", taxa_amr)) %>%
  arrange(match(colnames(count_merge), rownames(count_merge))) %>%
  mutate(source_group = if_else(source == "Bacteria", source, group)) %>%
  relocate(c(group, source_group), .after = source) %>% as.data.frame()

# Save Merged Information
#------------------------------------------------------------------------------#
gs <- mgnet(
  abundance = count_merge,
  info_sample = info_sample,
  info_taxa = taxa_merge
)

saveRDS(gs, "../data/R/gs_genera_amr.rds")
