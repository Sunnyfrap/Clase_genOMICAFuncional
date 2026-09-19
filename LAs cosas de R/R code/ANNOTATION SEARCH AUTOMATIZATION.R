##GENE ANNOTATION SEARCH AUTOMATIZATION

install.packages("rentrez")

library(rentrez)

# 1. Cargar tus archivos de no anotados reales
catuai_sin_annot <- read.csv("Catuai_Sin_Anotacion_Con_Etiqueta.csv")
cr95_sin_annot <- read.csv("CR95_Sin_Anotacion_Con_Etiqueta.csv")

# 2. Función automática para consultar NCBI usando los gene_id (LOCs) reales
rescatar_anotacion_ncbi <- function(df_no_anotados, nombre_variedad) {
  
  # Extraer los IDs reales de tu columna gene_id
  ids <- unique(df_no_anotados$gene_id)
  
  # Como NCBI a veces se satura si mandas miles de golpe, los procesamos en bloques de 200
  bloques <- split(ids, ceiling(seq_along(ids) / 200))
  resultados_totales <- list()
  
  message(paste("Consultando NCBI en automático para:", nombre_variedad))
  
  for (i in seq_along(bloques)) {
    sub_ids <- bloques[[i]]
    
    # Consulta a la base de datos "gene" de NCBI
    busqueda <- tryCatch({
      entrez_search(db = "gene", term = paste(sub_ids, collapse = " OR "), retmax = length(sub_ids))
    }, error = function(e) { NULL })
    
    if (!is.null(busqueda) && length(busqueda$ids) > 0) {
      resumenes <- entrez_summary(db = "gene", id = busqueda$ids)
      
      # Extraer descripción oficial
      if (is.list(resumenes) && !is.data.frame(resumenes)) {
        temp_df <- do.call(rbind, lapply(resumenes, function(x) {
          data.frame(
            gene_id = ifelse(is.null(x$name), NA, x$name),
            descripcion_ncbi = ifelse(is.null(x$description), NA, x$description),
            stringsAsFactors = FALSE
          )
        }))
        resultados_totales[[i]] <- temp_df
      }
    }
  }
  
  if (length(resultados_totales) == 0) {
    message(paste("No se encontraron registros nuevos en NCBI para", nombre_variedad))
    return(NULL)
  }
  
  tabla_final <- do.call(rbind, resultados_totales)
  return(tabla_final)
}

# 3. Ejecución directa para tus dos archivos de no anotados
anotados_ncbi_Catuai <- rescatar_anotacion_ncbi(catuai_sin_annot, "Catuai")
anotados_ncbi_CR95 <- rescatar_anotacion_ncbi(cr95_sin_annot, "CR95")

# Puedes revisar si NCBI guardó algo útil con:
head(anotados_ncbi_Catuai)
head(anotados_ncbi_CR95)
write.table(anotados_ncbi_Catuai, file = "anotados_ncbi_Catuai.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(anotados_ncbi_CR95, file = "anotados_ncbi_CR95.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
