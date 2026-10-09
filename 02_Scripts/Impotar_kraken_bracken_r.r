# Como pasar los datos de kraken y bracken a R

# =============== (1) =======================================

# De kraken se usan los archivos que son .k2report
# De bracken se usan los archivos .bracken

# procesamiento .bracken
  # Solo se carga como base de datos separada por tabuladores


# =============== (2) =======================================


# procesamiento kraken: Añadir los nombres de las columnas
  # colnames(Control01.k2report) <- c(
  # "percentage",
  # "clade_reads",
  # "direct_reads",
  # "minimizer_count",
  # "distinct_minimizers",
  # "rank_code",
  # "taxid",
  # "name"
  # )

# =============== (3) =======================================

# confirmar que taxid y taxonomy id son compatibles
  # Control01.k2report[which(Control01.k2report$taxid == 9), ] 
  # Control01.bracken[which(Control01.bracken$taxonomy_id ==9),] 


# =============== (4) =======================================

# Generar los objetos que contendran la taxonomia y la abundancia y despues hacer un join

# taxonomia <- Control01.k2report %>%
#   select(taxid, rank_code, name)


# abundancia <- Control01.bracken %>%
#  select(taxonomy_id, taxonomy_lvl ,new_est_reads)


# Control01 <- abundancia %>%
#  left_join(
#    taxonomia,
#    by = c("taxonomy_id" = "taxid")
#  )


# head(Control01)

# =============== (5) =======================================

# Revizar las multiples variables que kraken generó
# Suele añadir bastantes categorias que cmplicarían el proceso y que son poco usadas,
# por lo que solo me quedo con las clasicas.

# unique(Control01.k2report$rank_code)


# =============== (6) =======================================

#Filtrar de acuerdo a las variables taxonomicas que me interesan

# C01_f <-Control01.k2report %>%
#  filter(rank_code %in% c("K", "P", "C", "O", "F", "G", "S"))
# View(C01_f)

# =============== (7) =======================================

# Generar la taxonomia
# tax_filt <- C01_f %>%
#  mutate(name = str_trim(name)) %>%
#  mutate(
#    kingdom = if_else(rank_code == "K", name, NA_character_),
#    phylum   = if_else(rank_code == "P", name, NA_character_),
#    class   = if_else(rank_code == "C", name, NA_character_),
#    order   = if_else(rank_code == "O", name, NA_character_),
#    family = if_else(rank_code == "F", name, NA_character_),
#    genus  = if_else(rank_code == "G", name, NA_character_),
#    species = if_else(rank_code == "S", name, NA_character_)
#  ) %>%
#  fill(
#    kingdom, phylum, class, order,
#    family, genus, species,
#    .direction = "down"
#  )

# View(tax_filt)

# unique(tax_filt$family) # Es para ver cuantas familias distintas hay

# print(head(tax_filt))




# ============================

tax_table <- tax_filt %>%
  select(taxid, kingdom, phylum, class, order, family, genus, species) %>%
  distinct(taxid, .keep_all = TRUE)
View(tax_table)


abundancia <- Control01.bracken %>%
  select(taxonomy_id, new_est_reads) %>%
  rename(taxid = taxonomy_id)
View(abundancia)

Control01 <- abundancia %>%
  left_join(tax_table, by = "taxid")
View(Control01)

# ==============================

otu_table <- Control01 %>%
  select(taxid, new_est_reads) %>%
  column_to_rownames("taxid") %>%
  as.matrix()

otu_table <- otu_table(
  otu_table,
  taxa_are_rows = TRUE
)


# ================================

tax_table <- Control01 %>%
  select(taxid, kingdom, phylum, class, order, family, genus, species) %>%
  distinct(taxid, .keep_all = TRUE) %>%
  column_to_rownames("taxid") %>%
  as.matrix()

tax_table <- tax_table(tax_table)

# ================================== 


dim(otu_table)
View(otu_table)

dim(tax_table)
View(tax_table)

all(rownames(otu_table) == rownames(tax_table))


