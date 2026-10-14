# Explorando phyloseq

library(phyloseq)
readRDS("03_Results/ps_metagenomico_variante.RDS") -> ps_metagenomico_variante



estimate_richness(ps_metagenomico_variante, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
plot_richness(ps_metagenomico, x = "Grupo", measures = c("Observed","Shannon", "Chao1", "Simpson"), color = "Grupo")


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

merge_samples(ps_metagenomico_variante, group = Grupo$Control) -> ps_m_control









