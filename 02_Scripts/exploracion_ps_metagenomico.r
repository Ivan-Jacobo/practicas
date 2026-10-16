# Explorando phyloseq

library(phyloseq)
readRDS("03_Results/objetos_phyloseq/ps_metagenomico_variante.RDS") -> ps_metagenomico_variante

sample_names(ps_metagenomico_variante) -> Nombres

sample_data(ps_metagenomico_variante) -> df_ps_meta

as.data.frame(df_ps_meta) -> df_ps_meta

intento <- data.frame(
  Nombres = Nombres,
  Grupo   = c(rep("Control", 50), rep("Cushing", 50))
)


rownames(intento) <- Nombres

sample_data(intento) -> sample_data(ps_metagenomico_variante)

View(otu_table(ps_metagenomico_variante))

#saveRDS(ps_metagenomico_variante, "03_Results/objetos_phyloseq/ps_metagenomico_variante.RDS")
        
        
estimate_richness(ps_metagenomico_variante, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
plot_richness(ps_metagenomico_variante, x = "Grupo", measures = c("Observed","Shannon", "Chao1", "Simpson"), color = "Grupo")
head(otu_table(ps_metagenomico_variante))

# Nota:
  # R manda un aviso que dice que no hay singletons y doubletons
    # Esto es debido a que en bracken se utiliza un argumento que filtra las lecturas para cada taxon
    # con un minimo de 10, por lo que todo singleton y doubleton son eliminados con el objetvo de dejar
    # unicamente a los taxones que sea más probable que verdaderamente esten ahí.



# 1) Colapsar conteos (NArm = FALSE conserva los que no tienen género asignado)
ps_genero <- tax_glom(ps_metagenomico_variante, taxrank = "Genus", NArm = FALSE)

# 2) Abundancia relativa por muestra
ps_genero_rel <- transform_sample_counts(ps_genero, function(x) x / sum(x))

# 3) Quedarse con los 15 géneros más abundantes para que el gráfico sea legible
top15 <- names(sort(taxa_sums(ps_genero_rel), decreasing = TRUE))[1:15]
ps_top <- prune_taxa(top15, ps_genero_rel)

plot_bar(ps_top, x = "Sample", fill = "Genus")





# ==============================================000

# Grafico mas bonito

library(ggplot2)


# 1) Colapsar conteos a género, conservando los taxones sin género asignado
ps_genero <- tax_glom(ps_metagenomico_variante, taxrank = "Genus", NArm = FALSE)

# 2) Abundancia relativa (cada muestra suma 1)
ps_genero_rel <- transform_sample_counts(ps_genero, function(x) x / sum(x))

# 3) Top 15 géneros y fusión del resto en un solo taxón "Otros"
top15 <- names(sort(taxa_sums(ps_genero_rel), decreasing = TRUE))[1:15]
otros <- setdiff(taxa_names(ps_genero_rel), top15)

ps_final <- merge_taxa(ps_genero_rel, otros, archetype = 1)  # suma sus abundancias
tax_table(ps_final)[otros[1], ] <- "Otros"                   # renombrar el taxón fusionado

# Verificación: todas las muestras deben sumar 1
summary(sample_sums(ps_final))

# 4) Gráfico
p <- plot_bar(ps_final, x = "Sample", fill = "Genus")
#p$layers <- list()   # quitar la capa original (tiene un borde grueso)

p <- p +
  geom_bar(stat = "identity", position = "stack") +          # borde negro delgado
  guides(fill = guide_legend(ncol = 5)) +
  theme(
    axis.text.x      = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x     = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.35, "cm")
  )

p








# Generación de GBA pero separando controles de cushing

head(sample_data(ps_metagenomico_variante))


# Asi se genera un vector de los nombres de las muestras
sample_names(ps_metagenomico_variante) -> nombres_muestras

# Este es un vector con los nombres de las muestras que cuentan con el patron que me interesa
muestras_a_mantener <- nombres_muestras[startsWith(nombres_muestras, "Control")]

# Este es un objeto phyloseq que contiene unicamente los sujetos que me interesan
ps_metagenomico_controles <- prune_samples(muestras_a_mantener, ps_metagenomico_variante)

ps_metagenomico_controles

#=================
saveRDS(ps_metagenomico_controles, "03_Results/ps_metagenomico_controles.RDS")
#=================



# Este es un vector con los nombres de las muestras que cuentan con el patron que me interesa
muestras_a_mantener <- nombres_muestras[startsWith(nombres_muestras, "CS")]

# Este es un objeto phyloseq que contiene unicamente los sujetos que me interesan
ps_metagenomico_cushing <- prune_samples(muestras_a_mantener, ps_metagenomico_variante)

ps_metagenomico_cushing

#================
saveRDS(ps_metagenomico_cushing, "03_Results/ps_metagenomico_cushing.RDS")
#================


















#===============================================================
# GBA controles
#===============================================================


# 1) Colapsar conteos a género, conservando los taxones sin género asignado
ps_genero_controles <- tax_glom(ps_metagenomico_controles, taxrank = "Genus", NArm = FALSE)

# 2) Abundancia relativa (cada muestra suma 1)
ps_genero_cont_rel <- transform_sample_counts(ps_genero_controles, function(x) x / sum(x))

# 3) Top 15 géneros y fusión del resto en un solo taxón "Otros"
top15 <- names(sort(taxa_sums(ps_genero_cont_rel), decreasing = TRUE))[1:15]
otros <- setdiff(taxa_names(ps_genero_cont_rel), top15)

ps_final_cont <- merge_taxa(ps_genero_cont_rel, otros, archetype = 1)  # suma sus abundancias
tax_table(ps_final_cont)[otros[1], ] <- "Otros"                   # renombrar el taxón fusionado

# Verificación: todas las muestras deben sumar 1
summary(sample_sums(ps_final_cont))

# 4) Gráfico
t <- plot_bar(ps_final_cont, x = "Sample", fill = "Genus")
#p$layers <- list()   # quitar la capa original (tiene un borde grueso)

t <- t +
  geom_bar(stat = "identity", position = "stack") +          # borde negro delgado
  guides(fill = guide_legend(ncol = 5)) +
  theme(
    axis.text.x      = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x     = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.35, "cm")
  )

t

















#===============================================================
# GBA cushing
#===============================================================



# 1) Colapsar conteos a género, conservando los taxones sin género asignado
ps_genero_cushing <- tax_glom(ps_metagenomico_cushing, taxrank = "Genus", NArm = FALSE)

# 2) Abundancia relativa (cada muestra suma 1)
ps_genero_cush_rel <- transform_sample_counts(ps_genero_cushing, function(x) x / sum(x))

# 3) Top 15 géneros y fusión del resto en un solo taxón "Otros"
top15 <- names(sort(taxa_sums(ps_genero_cush_rel), decreasing = TRUE))[1:15]
otros <- setdiff(taxa_names(ps_genero_cush_rel), top15)

ps_final_cush <- merge_taxa(ps_genero_cush_rel, otros, archetype = 1)  # suma sus abundancias
tax_table(ps_final_cush)[otros[1], ] <- "Otros"                   # renombrar el taxón fusionado

# Verificación: todas las muestras deben sumar 1
summary(sample_sums(ps_final_cush))

# 4) Gráfico
g <- plot_bar(ps_final_cush, x = "Sample", fill = "Genus")
#p$layers <- list()   # quitar la capa original (tiene un borde grueso)

g <- g +
  geom_bar(stat = "identity", position = "stack") +          # borde negro delgado
  guides(fill = guide_legend(ncol = 5)) +
  theme(
    axis.text.x      = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x     = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.35, "cm")
  )

g



#================================================================================0

























































# ===================================================================
#Para diversidad beta y separación entre muestras:
  
ps_genero <- tax_glom(ps_metagenomico_variante, taxrank = "Genus", NArm = FALSE)
ps_genero_rel <- transform_sample_counts(ps_genero, function(x) x / sum(x))

head(sample_data(ps_genero_rel))

  # Distancia Bray-Curtis
dist_bc_all <- distance(ps_metagenomico_variante, method = "bray")

# PCoA
ord_pcoa_all <- ordinate(ps_genero_rel, method = "PCoA", distance = dist_bc_all)

all(rownames(ord_pcoa_all$vectors) %in% sample_names(ps_genero_rel))
# sale TRUE
dim(ord_pcoa_all$vectors)                 # debe ser 100 filas
head(rownames(ord_pcoa_all$vectors))      # si sale NULL, ese es el problema
head(attr(dist_bc_all, "Labels"), 60)         # ¿la distancia conserva los nombres?
find("distance")                          # debe salir solo "package:phyloseq"

# Visualización
plot_ordination(ps_genero_rel, ord_pcoa_all, color = "Grupo" )

#También puedes usar otros métodos de ordenación:
  
  # NMDS con Bray-Curtis
ord_nmds <- ordinate(ps_genero_rel, method = "NMDS", distance = "bray")
plot_ordination(ps_genero_rel, ord_nmds, color = "Grupo", shape =
                  "Grupo")

#Otras funciones prácticas:

sample_data(ps_genero_rel)
  
plot_heatmap(ps_genero_rel, sample.label = "Grupo")

plot_tree(ps_genero_rel, color = "Grupo")

plot_network(ps_genero_rel, distance = "bray", max.dist = 0.3)

sample_sums(ps_genero_rel)  # profundidad de secuenciación por muestra
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


