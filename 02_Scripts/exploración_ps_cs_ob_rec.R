#        Exploración phyloseq obesos, cushing y recuperados

library(phyloseq)
library(tidyverse)

readRDS("03_Results/ps_cs_lg_ob_cs_rec.RDS") -> ps_cs_lg_ob_cs_rec


# ==========================================================
ps_cs_lg_ob_cs_rec
sample_data(ps_cs_lg_ob_cs_rec) -> metadatos_df

as_data_frame(metadatos_df) -> metadatos_df
class(metadatos_df)
glimpse(metadatos_df)

metadatos_df <- metadatos_df %>%
  mutate(
    fenotipo = case_when(
      str_detect(host_phenotype, regex("active", ignore_case = TRUE)) ~
        "Cushing",
      str_detect(host_phenotype, regex("remission", ignore_case = TRUE)) ~
        "Remision",
      str_detect(host_phenotype, regex("obesity", ignore_case = TRUE)) ~
        "Obesos",
      TRUE ~ NA_character_
    )
  )


# ==================================================================
as.data.frame(metadatos_df) -> meta
meta$Run -> rownames(meta)

meta <- meta[sample_names(ps_cs_lg_ob_cs_rec), ,drop = FALSE]

sample_data(ps_cs_lg_ob_cs_rec) <- sample_data(meta)

View(sample_data(ps_cs_lg_ob_cs_rec))


 


# ====================================================================================



# riqueza
plot_richness(ps_cs_lg_ob_cs_rec, x = "fenotipo", measures = c("Chao1", "Shannon", "Simpson"), color = "host_sex")

plot_richness(ps_cs_lg_ob_cs_rec, x = "fenotipo", measures = NULL , color = "host_sex")



