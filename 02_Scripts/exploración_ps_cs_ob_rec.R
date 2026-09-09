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
plot_richness(ps_cs_lg_ob_cs_rec, x = "fenotipo", measures = c("Chao1", "Shannon", "Simpson"), color = "Host_age")

plot_richness(ps_cs_lg_ob_cs_rec, x = "fenotipo", measures = NULL , color = "Host_age")




# =======================================================================================

# Valores numéricos de riqueza/diversidad alfa
estimate_richness(ps_cs_lg_ob_cs_rec, measures = c("Observed", "Chao1", "Shannon", "Simpson"))

# Abundancia relativa por muestra
ps_rel_cs <- transform_sample_counts(ps_cs_lg_ob_cs_rec, function(x) x / sum(x))

# Composición taxonómica: barras apiladas
plot_bar(ps_rel_cs, x = "fenotipo", fill = "Phylum")


  
ps_ord_cs <- tax_glom(ps_rel_cs, taxrank = "Order")

plot_bar(ps_ord_cs, x = "fenotipo", fill = "Order")

#Para diversidad beta y separación entre muestras:
  
  # Distancia Bray-Curtis
dist_bc <- distance(ps_rel_cs, method = "bray")

# PCoA
ord_pcoa <- ordinate(ps_rel_cs, method = "PCoA", distance = dist_bc)

# Visualización
plot_ordination(ps_rel_cs, ord_pcoa, color = "fenotipo")

#También puedes usar otros métodos de ordenación:
  
  # NMDS con Bray-Curtis
ord_nmds <- ordinate(ps_rel_cs, method = "NMDS", distance = "bray")
plot_ordination(ps_rel_cs, ord_nmds, color = "fenotipo", shape =
                  "host_sex")

#Otras funciones prácticas:
  
plot_heatmap(ps_rel_cs, sample.label = "Sample_name")

plot_tree(ps_rel_cs, color = "host_phenotype")

plot_network(ps_rel_cs, distance = "bray", max.dist = 0.3)

sample_sums(ps_cs_lg_ob_cs_rec)  # profundidad de secuenciación por muestra
taxa_sums(ps)    # abundancia total por taxón

# ===========================================================================================#





subset_samples(ps_cs_lg_ob_cs_rec, fenotipo == "Cushing") -> ps_cs_activo
ps_cs_activo
#saveRDS(ps_cs_activo, "03_Results/ps_cs_activo.RDS")

subset_samples(ps_cs_lg_ob_cs_rec, fenotipo == "Remision") -> ps_cs_remision
ps_cs_remision
#saveRDS(ps_cs_remision, "03_Results/ps_cs_remision.RDS")

subset_samples(ps_cs_lg_ob_cs_rec, fenotipo == "Obesos") -> ps_obesos
ps_obesos
#saveRDS(ps_obesos, "03_Results/ps_obesos.RDS")


