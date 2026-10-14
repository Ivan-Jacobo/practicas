# ============================================================
# KRAKEN2 + BRACKEN → PHYLOSEQ
# Flujo de trabajo para una sola muestra
# ============================================================

library(dplyr)
library(tidyverse)
library(phyloseq)



# ===================== (1) =====================
# Cargar archivos

Control01.k2report <- read.delim(
  "01_Raw_data/Control01.k2report",
  header = FALSE,
  sep = "\t",
  check.names = FALSE
)

Control01.bracken <- read.delim(
  "01_Raw_data/Control01.bracken",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)


#####



bracken_seleccion <- lapply(names(bracken_lista), function(muestra) {
  df <- bracken_lista[[muestra]]
  
  df %>%
    select(taxonomy_id, name, taxonomy_lvl, new_est_reads) %>%
    rename(!!muestra := new_est_reads)
})

bracken_combinado <- Reduce(
  function(x, y) full_join(
    x, y,
    by = c("taxonomy_id", "name", "taxonomy_lvl")
  ),
  bracken_seleccion
)

 adw

# ===================== (2) =====================
# Nombrar las columnas del reporte de Kraken2

colnames(Control01.k2report) <- c(
  "percentage",
  "clade_reads",
  "direct_reads",
  "minimizer_count",
  "distinct_minimizers",
  "rank_code",
  "taxid",
  "name"
)


# ===================== (3) =====================
# Filtrar los niveles taxonómicos que nos interesan

C01_f <- Control01.k2report %>%
  filter(rank_code %in% c("K", "P", "C", "O", "F", "G", "S"))


# ===================== (4) =====================
# Reconstruir la taxonomía

tax_filt <- C01_f %>%
  mutate(name = str_trim(name)) %>%
  mutate(
    kingdom = if_else(rank_code == "K", name, NA_character_),
    phylum  = if_else(rank_code == "P", name, NA_character_),
    class   = if_else(rank_code == "C", name, NA_character_),
    order   = if_else(rank_code == "O", name, NA_character_),
    family  = if_else(rank_code == "F", name, NA_character_),
    genus   = if_else(rank_code == "G", name, NA_character_),
    species = if_else(rank_code == "S", name, NA_character_)
  ) %>%
  tidyr::fill(
    kingdom, phylum, class, order,
    family, genus, species,
    .direction = "down"
  )


# ===================== (5) =====================
# Extraer abundancias de Bracken

abundancia <- Control01.bracken %>%
  select(taxonomy_id, new_est_reads) %>%
  rename(taxid = taxonomy_id)


# ===================== (6) =====================
# Unir abundancia y taxonomía

Control01 <- abundancia %>%
  left_join(
    tax_filt %>%
      select(taxid, kingdom, phylum, class, order, family, genus, species) %>%
      distinct(taxid, .keep_all = TRUE),
    by = "taxid"
  )


# ===================== (7) =====================
# Crear OTU table

otu <- Control01 %>%
  select(taxid, new_est_reads) %>%
  tibble::column_to_rownames("taxid") %>%
  as.matrix()

otu <- phyloseq::otu_table(
  otu,
  taxa_are_rows = TRUE
)

View(otu)

# ===================== (8) =====================
# Crear tax table

tax <- Control01 %>%
  select(
    taxid, kingdom, phylum, class,
    order, family, genus, species
  ) %>%
  distinct(taxid, .keep_all = TRUE) %>%
  tibble::column_to_rownames("taxid") %>%
  as.matrix()

tax <- phyloseq::tax_table(tax)

View(tax)
# ===================== (9) =====================
# Comprobaciones finales

dim(otu)
dim(tax)

all(rownames(otu) == rownames(tax))




















# ===================================================================== #
# Lo mismo pero ahora para muchos archivos

archivo <- "01_Raw_data/Controles_Cushing.KREPORT"

controles_cushing <- read.delim(
  archivo,
  skip = 102,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

View(controles_cushing)

# Nota: En el archivo .KREPORT me dice que que S1 al S50 son los controles del 1 al 50
# mientras que para los S51 al S100 son cushing del 1 al 50.


directorio <- "01_Raw_data/Bracken"

archivos_bracken <- list.files(
  path = directorio,
  pattern = "^(Control|CS)[0-9]+\\.bracken$",
  full.names = TRUE
)

bracken_lista <- lapply(archivos_bracken, read.delim)

names(bracken_lista) <- sub(
  "\\.bracken$", "",
  basename(archivos_bracken)
)

# Verificar archivos cargados
length(bracken_lista)
head(names(bracken_lista), 101)


#=======================================


library(dplyr)



archivos_bracken <- list.files(
  path = directorio,
  pattern = "^(Control|CS)[0-9]+\\.bracken$",
  full.names = TRUE
)

bracken_lista <- lapply(archivos_bracken, read.delim)

names(bracken_lista) <- sub(
  "\\.bracken$", "",
  basename(archivos_bracken)
)




bracken_seleccion <- lapply(names(bracken_lista), function(muestra) {
  
  df <- bracken_lista[[muestra]]
  
  df %>%
    select(taxonomy_id, name, taxonomy_lvl, new_est_reads) %>%
    rename(!!muestra := new_est_reads)
  
})

bracken_combinado <- Reduce(
  function(x, y) {
    full_join(
      x, y,
      by = c("taxonomy_id", "name", "taxonomy_lvl")
    )
  },
  bracken_seleccion
)

View(bracken_combinado)



# ===================00


# Reemplazar NA por 0 solo en las columnas de abundancia
columnas_muestras <- grep(
  "^(Control|CS)[0-9]+$",
  names(bracken_combinado),
  value = TRUE
)

bracken_combinado <- bracken_combinado %>%
  mutate(
    across(all_of(columnas_muestras), ~ tidyr::replace_na(.x, 0))
  )

# ===============



# ===================== (2) =====================

####
colnames(controles_cushing)
head(controles_cushing[, c("lvl_type", "taxid", "name")])


taxonomia_kreport <- controles_cushing %>%
  select(lvl_type, taxid, name) %>%
  distinct(taxid, .keep_all = TRUE)

# Revisar resultado
dim(taxonomia_kreport)
head(taxonomia_kreport)



# ===================== (3) =====================
# Filtrar los niveles taxonómicos que nos interesan

controles_cushing_f <- controles_cushing %>%
  filter(lvl_type %in% c("K", "P", "C", "O", "F", "G", "S"))

head(data.frame(
  taxid = controles_cushing$taxid,
  lvl_type = controles_cushing$lvl_type,
  name = controles_cushing$name,
  espacios = nchar(controles_cushing$name) -
    nchar(trimws(controles_cushing$name))
), 15)

















# 1)

library(dplyr)
library(stringr)
library(tidyr)

tax_completo <- controles_cushing %>%
  mutate(
    taxid = as.character(taxid),
    name_original = name,
    indent = nchar(name) - nchar(str_trim(name)),
    name = str_trim(name)
  )

# Inicializar columnas taxonómicas
rangos <- c("K", "P", "C", "O", "F", "G", "S")
columnas <- c("kingdom", "phylum", "class", "order",
              "family", "genus", "species")

for (col in columnas) {
  tax_completo[[col]] <- NA_character_
}

# Reconstruir ancestros usando la indentación
pila <- list()

for (i in seq_len(nrow(tax_completo))) {
  
  nivel <- tax_completo$indent[i]
  
  # Retirar ancestros del mismo nivel o inferiores
  pila <- pila[
    vapply(pila, function(x) x$indent < nivel, logical(1))
  ]
  
  # Agregar taxón actual
  pila[[length(pila) + 1]] <- list(
    indent = nivel,
    lvl_type = tax_completo$lvl_type[i],
    name = tax_completo$name[i]
  )
  
  # Asignar los ancestros según su rango
  for (j in seq_along(pila)) {
    rango <- match(pila[[j]]$lvl_type, rangos)
    
    if (!is.na(rango)) {
      tax_completo[[columnas[rango]]][i] <- pila[[j]]$name
    }
  }
}



#2)

tax_filt <- tax_completo %>%
  filter(lvl_type %in% rangos) %>%
  select(
    taxid, name, lvl_type,
    kingdom, phylum, class, order,
    family, genus, species
  )






#3)
bracken_combinado <- bracken_combinado %>%
  mutate(taxonomy_id = as.character(taxonomy_id))

bracken_anotado <- bracken_combinado %>%
  left_join(
    tax_filt,
    by = c("taxonomy_id" = "taxid"),
    suffix = c("_bracken", "_kreport")
  )



# 4)
dim(bracken_anotado)
head(bracken_anotado)
View(bracken_anotado)















# 5)


library(dplyr)
library(tibble)
library(phyloseq)


# Crear matriz de abundancias
otu <- bracken_anotado %>%
  select(
    taxonomy_id,
    matches("^(Control|CS)[0-9]+$")
  ) %>%
  distinct(taxonomy_id, .keep_all = TRUE) %>%
  column_to_rownames("taxonomy_id") %>%
  as.matrix()

# Crear objeto phyloseq
otu <- otu_table(
  otu,
  taxa_are_rows = TRUE
)

# Verificar
dim(otu)
View(otu)

sum(is.na(otu))
sum(!is.finite(otu))



# 6)

# Crear matriz taxonómica
tax <- bracken_anotado %>%
  select(
    taxonomy_id,
    kingdom, phylum, class,
    order, family, genus, species
  ) %>%
  distinct(taxonomy_id, .keep_all = TRUE) %>%
  column_to_rownames("taxonomy_id") %>%
  as.matrix()

# Convertir a tax_table
tax <- tax_table(tax)

# Visualizar
View(as.data.frame(tax))





# 7)

ntaxa(otu)
ntaxa(tax)

all(taxa_names(otu) == taxa_names(tax))




# 8)

# Identificar las muestras
muestras <- sample_names(otu)

# Crear metadatos
metadata <- data.frame(
  Grupo = ifelse(grepl("^Control", muestras), "Control", "Cushing"),
  row.names = muestras
)

# Convertir a sample_data
sampledata <- sample_data(metadata)

# Revisar
View(as.data.frame(sampledata))
table(sampledata$Grupo)

all(sample_names(otu) == rownames(metadata))



# 9)

# Ensamblar objeto phyloseq
ps_metagenomico <- phyloseq(
  otu,
  tax,
  sampledata
)

# Verificar el objeto
ps_metagenomico

saveRDS(ps_metagenomico, "03_Results/ps_metagenomico.RDS")
