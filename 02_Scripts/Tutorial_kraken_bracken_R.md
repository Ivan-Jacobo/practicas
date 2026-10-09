---
title: "Desde kraken y bracken hasta R"
author: "Ivan jacobo"
date: "2026-10-08"
output: html_document
---

## Selección de archivos

### Kraken2

Kraken 2 es un clasificador taxonómico. Examina los k-mers de cada secuencia, los consulta en una base de datos que mapea cada k-mer al ancestro común más bajo (LCA) de todos los genomas que lo contienen, y con eso asigna una etiqueta taxonómica a la lectura.

Al usar Kraken se generan 3 archivos

1) archivo.kraken

2) archivo.k2report

3) archivo.log

De estos 3 archivos el que nos interesa es el .k2report, el cual cuenta con la siguiente estructura

```{r}

print(head(Control01.k2report))

print(head(Control01.k2report$percentage))
# Porcentaje de fragmentos en el clado con raíz en este taxón

print(head(Control01.k2report$clade_reads))
# Fragmentos en el taxón y todos sus descendientes (acumulado)

print(head(Control01.k2report$direct_reads))
# Fragmentos asignados exactamente a este taxón

print(head(Control01.k2report$rank_code))
# U, R, D, K, P, C, O, F, G, S; los intermedios llevan un número (G2, R1, S1)

print(head(Control01.k2report$taxid))
# Identificador taxonomico

print(head(Control01.k2report$name))
# Nombre científico, indentado según la jerarquía

```

De estas columnas, las mas relevantes son las de taxid, rank_code y name.

```{r}
taxonomia <- Control01.k2report %>%
   select(taxid, rank_code, name)
```

### Bracken

**Bracken** (Bayesian Reestimation of Abundance with KrakEN) existe por una limitación de Kraken. Cuando una lectura es idéntica en varias especies, la asignación correcta es el LCA (por ejemplo el género), así que muchas lecturas quedan “varadas” por encima de especie y el conteo por especie queda subestimado. Bracken usa probabilidades derivadas de la propia base de datos de Kraken para redistribuir esas lecturas hacia los taxones inferiores y estimar abundancia a nivel de especie, género o superior

Bracken no corre solo con el reporte. Requiere un archivo de distribución de k-mers específico para la base de datos y la longitud de lectura, generado con `bracken-build`: se clasifica la biblioteca completa contra su propia base de datos, se calcula qué taxón asignaría Kraken a cada lectura perfecta de longitud `-l` extraída de cada genoma, y se resume en `databaseXmers.kmer_distrib`

```{r}
print(head(Control01.bracken))

print(head(Control01.bracken$name))
# Nombre del taxon

print(head(Control01.bracken$taxonomy_id))
# tax id de NCBI

print(head(Control01.bracken$taxonomy_lvl))
# Nivel taxonomico estimado

print(head(Control01.bracken$kraken_assigned_reads))
# Lecturas que kraken asigno al taxon

print(head(Control01.bracken$added_reads))
# Lecturas añadidas a partir de una reestimación

print(head(Control01.bracken$new_est_reads))
# Estimación final

print(head(Control01.bracken$fraction_total_reads))
# proporción de lecturas que fueron añadidas por nec est reads con respecto de las que kraken había estimado
# -> new_est_reads / kraken_assigned_reads


```

Por lo tanto del archivo .bracken unicamente mantenemos las siguientes columnas: taxonomy_id, taxonomy_lvl ,new_est_reads.daw

```{r}
abundancia <- Control01.bracken %>%
  select(taxonomy_id, taxonomy_lvl ,new_est_reads)
```

```{r}
Control01 <- abundancia %>%
  left_join(
    taxonomia,
    by = c("taxonomy_id" = "taxid")
  )
```

Filtrar de acuerdo a las variables taxonomicas que me interesan

```{r}


C01_f <-Control01.k2report %>%
  filter(rank_code %in% c("K", "P", "C", "O", "F", "G", "S"))
 print(head(C01_f))
```

Generar la taxonomia

```{r}

 tax_filt <- C01_f %>%
  mutate(name = str_trim(name)) %>%
  mutate(
    kingdom = if_else(rank_code == "K", name, NA_character_),
    phylum   = if_else(rank_code == "P", name, NA_character_),
    class   = if_else(rank_code == "C", name, NA_character_),
    order   = if_else(rank_code == "O", name, NA_character_),
    family = if_else(rank_code == "F", name, NA_character_),
    genus  = if_else(rank_code == "G", name, NA_character_),
    species = if_else(rank_code == "S", name, NA_character_)
  ) %>%
  fill(
    kingdom, phylum, class, order,
    family, genus, species,
    .direction = "down"
  )
 
 
print(head(tax_filt)) 
```
