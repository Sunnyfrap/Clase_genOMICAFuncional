

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

camera.catuai.plot <-
  plot.camera.barplot(
    sig_results$Catuai,
    "Catuai",
    top_n = 50
  )

camera.CR95.plot <-
  plot.camera.barplot(
    sig_results$CR95,
    "CR95",
    top_n = 50
  )


print(camera.catuai.plot)
print(camera.CR95.plot)



camera.catuai.plot <-
  plot.camera.barplot(
    sig_results$Catuai,
    "Catuai",
    top_n = 1000
  )

camera.cr95.plot <-
  plot.camera.barplot(
    sig_results$CR95,
    "CR95",
    top_n = 1000
  )



##############################################################
# GO TERMS COMPARTIDOS
##############################################################

shared.go <- intersect(
  catuai.top$GOID,
  cr95.top$GOID
)


##############################################################
# QUEDARNOS SOLO CON LOS COMPARTIDOS
##############################################################

catuai.shared <- catuai.top[
  catuai.top$GOID %in% shared.go,
  ,
  drop = FALSE
]

cr95.shared <- cr95.top[
  cr95.top$GOID %in% shared.go,
  ,
  drop = FALSE
]


##############################################################
# COMBINAR
##############################################################

shared.camera <- rbind(
  catuai.shared,
  cr95.shared
)

shared.camera$cultivar <- c(
  rep("Catuai", nrow(catuai.shared)),
  rep("CR95", nrow(cr95.shared))
)


##############################################################
# VALOR
##############################################################

shared.camera$value <- ifelse(
  shared.camera$Direction == "Up",
  -log10(shared.camera$FDR),
  log10(shared.camera$FDR)
)


##############################################################
# ORDENAR GO TERMS
##############################################################

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
# GRÁFICA
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
    title = "CAMERA: Shared GO terms",
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

shared.camera.plot
##############################################################
# SHOW
##############################################################


print(shared.camera.plot)

