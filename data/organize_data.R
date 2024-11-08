library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

count_species <- read_csv("RAW/mOTU_gs3_counts.csv", name_repair = "minimal")
count_genera <- read_csv("RAW/genus_gs3_counts.csv", name_repair = "minimal")
taxa_genera_raw <- read_csv("RAW/genus_gs3_taxa_paths.csv", name_repair = "minimal")
taxa_species_raw <- read_csv("RAW/mOTU_gs3_taxa_paths.csv", name_repair = "minimal")

df_amr <- read_csv("RAW/counts_clusters.csv")
#taxa_amr <- readxl::read_excel("RAW/panres_annotations.xlsx", sheet = "wide_format")
#taxa_amr <- read_tsv("RAW/panres_annotations.tsv")
info_sample <- read_csv("RAW/metadata_fixed.csv") %>% select(-1) 
taxa_amr <- read_csv("RAW/chosen_gene_metadata.csv")

dict_ncbi_gtbtk <- read_tsv("RAW/mOTUs_ncbi_gtdb_species.tsv")
uhgg <- read_csv("RAW/uhgg_gtdb_r220.csv")

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

# taxa_amr <- df_amr %>%
#   select(cluster_representative_98, group) %>%
#   distinct() %>%
#   rename(cluster_id = cluster_representative_98)
# 
# taxa_amr <- readxl::read_xlsx("RAW/panres_annotations.xlsx", sheet = "wide_format") %>%
#   filter(!is.na(cluster_representative)) %>%
#   select(-c(fa_name,gene_length,database)) %>%
#   distinct() %>%
#   separate(cluster_representative, sep = "_", remove = FALSE,
#            into = c("temp1", "temp2", "cluster_version")) %>%
#   unite("cluster_id", temp1, temp2, sep = "_") %>%
#   group_by(cluster_id) %>%
#   filter(!is.na(class)) %>%
#   filter(nchar(class) == max(nchar(class))) %>%
#   slice(1) %>%
#   ungroup() %>%
#   right_join(taxa_amr, by = "cluster_id") %>%
#   relocate(group) %>%
#   column_to_rownames("cluster_id")

# tmp <- taxa_amr %>%
#   filter(cluster_representative_98 %in% colnames(count_amr)) %>%
#   # group_by(cluster_representative_98) %>%
#   # filter(database == "resfinder") %>%
#   # ungroup() %>%
#   select(cluster_representative_98, group, class, database, resistance_type) %>%
#   distinct() %>%
#   group_by(cluster_representative_98) %>%
#   filter((cluster_representative_98 == "pan_188" & class == "['polymyxin']" & 
#             resistance_type == "antimicrobial" & database == "resfinder") | 
#            # Preserve all other rows that are not `pan_188`
#            cluster_representative_98 != "pan_188") %>%
#   filter((cluster_representative_98 == "pan_4862" & class == "['quinolone']") | 
#            # Preserve all other rows that are not `pan_4862`
#            cluster_representative_98 != "pan_4862") %>%
#   filter((cluster_representative_98 == "pan_1655" & class == "['quinolone']") | 
#            # Preserve all other rows that are not `pan_1655`
#            cluster_representative_98 != "pan_1655")
  
# taxa_amr <- tibble(taxa_id = colnames(count_amr)) %>%
#   left_join(taxa_amr, by = join_by(taxa_id == cluster_representative_98)) %>%
#   column_to_rownames("taxa_id")

taxa_amr_clean <- taxa_amr %>%
  select(-ResFinder, -ResFinderFG, -csabapal, -fa_name) %>%
  mutate(database_priority = factor(database, levels = c("resfinder", "functional_amr", "csabapal"), ordered = TRUE)) %>%
  group_by(cluster_representative_98) %>%
  filter(database_priority == max(database_priority)) %>%
  ungroup() %>%
  select(-database_priority) %>%
  distinct() %>%
  select(-group) %>%
  distinct() %>%
  select(-gene_length) %>%
  distinct() %>%
  select(-cluster_representative_90) %>%
  distinct() %>%
  mutate(class_length = nchar(gsub("\\W", "", class))) %>%
  group_by(cluster_representative_98) %>%
  slice_max(class_length) %>%
  select(-class_length) %>%
  mutate(
    class = sapply(class, function(x) {
      # Extract elements between brackets, split by comma, trim spaces, sort, and reconstruct
      elements <- str_extract(x, "\\[.*\\]") %>%
        gsub("\\[|\\]|'", "", .) %>%
        str_split(",") %>%
        unlist() %>%
        trimws() %>%
        sort()
      sprintf("['%s']", paste(elements, collapse = "', '"))
    })) %>%
  distinct() %>%
  select(-resistance_type) %>%
  distinct() %>%
  mutate(group = if_else(database == "resfinder", "Resfinder", "Functional")) %>%
  column_to_rownames("cluster_representative_98") 

colnames(taxa_amr_clean) <- paste("amr", colnames(taxa_amr_clean), sep = "_")

taxa_amr_clean <- taxa_amr_clean %>%
  mutate(source = "amr", .before = 1)


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

# Example usage where you iterate over each element of species_gtdb and compare it to all of uhgg$SPECIES
dict_ncbi_gtbtk <- dict_ncbi_gtbtk %>%
  mutate(
    human_gut = map_chr(species_gtdb, ~ {
      matching_species <- uhgg$SPECIES[str_detect(.x, fixed(uhgg$SPECIES))]
      if (length(matching_species) > 0) {
        return(matching_species[1])  # Return the first match if found
      } else {
        return(NA_character_)  # Return NA if no match
      }
    })
  )

taxa_species <- taxa_species %>%
  rownames_to_column("taxa_id") %>%
  left_join(dict_ncbi_gtbtk, by = join_by(taxa_id == mOTU_ncbi_name)) %>%
  mutate(is_from_gut = if_else(is.na(human_gut), FALSE, TRUE)) %>%
  group_by(taxa_id) %>%
  slice(1) %>%
  ungroup() %>%
  column_to_rownames("taxa_id")

#count_species <- count_species[, rownames(taxa_species)]
#count_species <- count_species[rownames(info_sample)]

# Example usage where you iterate over each element of species_gtdb and compare it to all of uhgg$SPECIES
# dict_ncbi_gtbtk <- dict_ncbi_gtbtk %>%
#   mutate(
#     human_gut = map_chr(species_gtdb, ~ {
#       matching_species <- uhgg$SPECIES[str_detect(.x, fixed(uhgg$SPECIES))]
#       if (length(matching_species) > 0) {
#         return(matching_species[1])  # Return the first match if found
#       } else {
#         return(NA_character_)  # Return NA if no match
#       }
#     })
#   )

# taxa_species <- taxa_species %>%
#   rownames_to_column("taxa_id") %>%
#   left_join(dict_ncbi_gtbtk, by = join_by(taxa_id == mOTU_ncbi_name)) %>%
#   mutate(is_from_gut = if_else(is.na(human_gut), FALSE, TRUE)) %>%
#   group_by(taxa_id) %>%
#   slice(1) %>%
#   ungroup() %>%
#   column_to_rownames("taxa_id")
  
# BACTERIA mOTU GENERA
#------------------------------------------------------------------------------#
count_genera <- count_genera %>%
  column_to_rownames("genepid") %>%
  as.matrix()

taxa_genera <- taxa_genera_raw %>%
  separate(path_taxs,
           into = c("kingdom_id","phylum_id","class_id","order_id","family_id","genus_id"),
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



# SAMPLE INTERSECTION
#------------------------------------------------------------------------------#
sample_intersection <- purrr::reduce(
  .x = list(rownames(count_amr), rownames(info_sample), rownames(count_species)),
  .f = intersect
)

sample_intersection_genus <- purrr::reduce(
  .x = list(rownames(count_amr), rownames(info_sample), rownames(count_genera)),
  .f = intersect
)

all(sample_intersection == sample_intersection_genus)

# SELECT THE INTERSECTION
# From 1240 to 1147 (fair enough)
#------------------------------------------------------------------------------#
count_amr <- count_amr[sample_intersection, ]
count_species <- count_species[sample_intersection, ]
count_genera <- count_genera[sample_intersection, ]
info_sample <- info_sample[sample_intersection,]

# MERGE TAXA AND COUNT INFORMATION IN AMR
#------------------------------------------------------------------------------#
all(colnames(count_amr) %in% rownames(taxa_amr_clean))
taxa_amr_clean <- taxa_amr_clean[colnames(count_amr), , drop = F]
#taxa_amr <- taxa_amr %>% mutate(source = "amr", .before = 1) %>%
#  rename(class_amr = class, group_amr = group, resistance_type_amr = resistance_type) 
# taxa_amr <- taxa_amr %>% mutate(source = "amr", .before = 1) %>%
#   rename(class_amr = class, group_amr = group, biocide_class_amr = biocide_class,
#          metal_class_amr = metal_class) 

# MERGE TAXA AND COUNT INFORMATION IN GENERA
#------------------------------------------------------------------------------#
all(rownames(taxa_species) %in% colnames(count_species))
taxa_species <- taxa_species[colnames(count_species), ]
taxa_species <- taxa_species %>% rename(source = level)

all(rownames(taxa_genera) %in% colnames(count_genera))
taxa_genera <- taxa_genera[colnames(count_genera), ]
taxa_genera <- taxa_genera %>% rename(source = level)

# CREATE THE MGNET OBJECTS
#------------------------------------------------------------------------------#
library(mgnet)

gs_mOTU <- mgnet(
  abun = count_species,
  meta = info_sample,
  taxa = taxa_species
)

gs_genera <- mgnet(
  abun = count_genera,
  meta = info_sample,
  taxa = taxa_genera
)

gs_amr <- mgnet(
  abun = as.matrix(count_amr),
  meta = info_sample,
  taxa = taxa_amr_clean
)

gs_mOTUs_amr <- mgnet(
  abun = cbind(as.matrix(count_species), as.matrix(count_amr)),
  meta = info_sample,
  taxa = bind_rows(taxa_species, taxa_amr_clean)
)

gs_genera_amr <- mgnet(
  abun = cbind(as.matrix(count_genera), as.matrix(count_amr)),
  meta = info_sample,
  taxa = bind_rows(taxa_genera, taxa_amr_clean)
)

gs_mOTUs_amr@taxa <- taxa(gs_mOTUs_amr) %>% mutate(type = case_when(
  source == "mOTU" ~ "mOTU",
  amr_group == "Resfinder" ~ "Resfinder",
  amr_group == "Functional" ~ "Functional",
), .after = 1) 

gs_genera_amr@taxa <- taxa(gs_genera_amr) %>% mutate(type = case_when(
  source == "genus" ~ "genus",
  amr_group == "Resfinder" ~ "Resfinder",
  amr_group == "Functional" ~ "Functional",
), .after = 1) 

# SAVE DATA IN MGNET FORMATS
#------------------------------------------------------------------------------#
saveRDS(gs_mOTUs_amr, "R/global-sewage-mOTUs-amr.rds")
saveRDS(gs_genera_amr, "R/global-sewage-genera-amr.rds")


