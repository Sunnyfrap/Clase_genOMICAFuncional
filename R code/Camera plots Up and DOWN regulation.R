
###NOTE TO SELF: BEFORE RUNNING THIS CODE, YOU NEED TO RUN THE PIPELINE UNTIL 
##POINT 28


##############################################################
# 29. CAMERA BARPLOTS
##############################################################

plot.camera.barplot <- function(
    res,
    title,
    top_n = 25
) {
  
  if (nrow(res) == 0) {
    return(
      ggplot() +
        labs(
          title = paste(
            title,
            "- no significant GO terms"
          )
        ) +
        theme_void()
    )
  }
  
  # Top UP y DOWN por separado
  res.up <- res[
    res$Direction == "Up",
    ,
    drop = FALSE
  ]
  
  res.down <- res[
    res$Direction == "Down",
    ,
    drop = FALSE
  ]
  
  res.up <- head(
    res.up[
      order(res.up$FDR),
      ,
      drop = FALSE
    ],
    top_n
  )
  
  res.down <- head(
    res.down[
      order(res.down$FDR),
      ,
      drop = FALSE
    ],
    top_n
  )
  
  res.plot <- rbind(
    res.up,
    res.down
  )
  
  # Up = positivo
  # Down = negativo
  res.plot$value <- ifelse(
    res.plot$Direction == "Up",
    -log10(res.plot$FDR),
    log10(res.plot$FDR)
  )
  
  # Ordenar términos
  res.plot$TERM <- factor(
    res.plot$TERM,
    levels = res.plot$TERM[
      order(res.plot$value)
    ]
  )
  
  ggplot(
    res.plot,
    aes(
      x = value,
      y = TERM,
      fill = Direction
    )
  ) +
    geom_col(
      width = 0.7
    ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.5
    ) +
    scale_fill_manual(
      values = c(
        Up = "#66A61E",
        Down = "#FFD92F"
      )
    ) +
    labs(
      title = title,
      x = "-log10(FDR)",
      y = "GO term",
      fill = "Direction"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 9
      ),
      legend.position = "top"
    )
}


##############################################################
# INDIVIDUAL PLOTS
##############################################################

camera.catuai.plot <- plot.camera.barplot(
  sig_results$Catuai,
  "Catuai",
  top_n = 50
)

camera.CR95.plot <- plot.camera.barplot(
  sig_results$CR95,
  "CR95",
  top_n = 50
)

print(camera.catuai.plot)
print(camera.CR95.plot)


##############################################################
# GO TERMS COMPARTIDOS (CORREGIDO DE RAÍZ)
##############################################################

# Sincronizamos las variables 'top' directamente de tu objeto sig_results
catuai.top <- sig_results$Catuai
cr95.top <- sig_results$CR95

# Detectar automáticamente si debemos usar 'GOID' o 'TERM' para cruzar los datos
goid_col <- if ("GOID" %in% colnames(catuai.top)) "GOID" else "TERM"

shared.go <- intersect(
  catuai.top[[goid_col]],
  cr95.top[[goid_col]]
)

# Mensaje de control en tu consola para auditoría visual rápida
cat("\n--- Verificación de Datos ---")
cat("\nTérminos en Catuai:", nrow(catuai.top))
cat("\nTérminos en CR95:", nrow(cr95.top))
cat("\nTérminos compartidos detectados:", length(shared.go), "\n\n")


##############################################################
# QUEDARNOS SOLO CON LOS COMPARTIDOS Y COMBINAR
##############################################################

if (length(shared.go) > 0) {
  
  catuai.shared <- catuai.top[
    catuai.top[[goid_col]] %in% shared.go,
    ,
    drop = FALSE
  ]
  
  cr95.shared <- cr95.top[
    cr95.top[[goid_col]] %in% shared.go,
    ,
    drop = FALSE
  ]
  
  # Asignar cultivares antes de unir las tablas
  catuai.shared$cultivar <- "Catuai"
  cr95.shared$cultivar <- "CR95"
  
  shared.camera <- rbind(
    catuai.shared,
    cr95.shared
  )
  
  ##############################################################
  # VALOR Y ORDENAMIENTO DE FACTORES
  ##############################################################
  
  shared.camera$value <- ifelse(
    shared.camera$Direction == "Up",
    -log10(shared.camera$FDR),
    log10(shared.camera$FDR)
  )
  
  # Ordenar los GO Terms según la magnitud de su FDR
  term.order <- unique(
    shared.camera$TERM[
      order(abs(shared.camera$value))
    ]
  )
  
  shared.camera$TERM <- factor(
    shared.camera$TERM,
    levels = rev(term.order)
  )
  
  ##############################################################
  # GRÁFICA COMPARTIDA
  ##############################################################
  
  shared.camera.plot <- ggplot(
    shared.camera,
    aes(
      x = value,
      y = TERM,
      fill = cultivar
    )
  ) +
    geom_col(
      position = position_dodge(
        width = 0.8
      ),
      width = 0.7
    ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.5
    ) +
    scale_fill_manual(
      values = c(
        Catuai = "#F8766D",
        CR95 = "#00BFC4"
      )
    ) +
    labs(
      title = paste("CAMERA:", length(shared.go), "Shared GO terms"),
      x = "-log10(FDR)",
      y = "GO term",
      fill = "Cultivar"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 8
      ),
      legend.position = "top"
    )
  
  print(shared.camera.plot)
  
} else {
  
  # Si el intersect falla, creamos una gráfica vacía de advertencia informativa
  shared.camera.plot <- ggplot() +
    labs(title = "CAMERA - No shared GO terms detected (Check column names)") +
    theme_void()
  
  print(shared.camera.plot)
}
