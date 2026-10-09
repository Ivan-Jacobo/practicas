# Explorando phyloseq

library(phyloseq)
readRDS("03_Results/phyloseq_ps_cushing.RDS") -> ps_cushing

estimate_richness(ps_cushing, measures = c("Observed", "Chao1", "Shannon", "Simpson"))


otu_mat <- as(otu_table(ps_cushing), "matrix")

sum(is.na(otu_mat))
sum(!is.finite(otu_mat))

head(otu_mat)

summary(as.vector(otu_mat))
any(otu_mat < 0, na.rm = TRUE)


which(is.na(otu_mat), arr.ind = TRUE)[1:10, ]





ps_cushing
#
View(sample_data(ps_cushing))
View(tax_table(ps_cushing))





# Abundancia relativa por muestra
met_cush_rel <- transform_sample_counts(ps_cushing, function(x) x / sum(x))

# Composición taxonómica: barras apiladas
plot_bar(met_cush_rel, x = "Grupo", fill = "kingdom")



ps_ord_cs <- tax_glom(ps_rel_cs, taxrank = "Order")

plot_bar(ps_ord_cs, x = "fenotipo", fill = "Order")


