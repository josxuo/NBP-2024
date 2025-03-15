#################
### BCS theme ###
#################

theme_bcs <- function() {
  theme(
    # backgrounds
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid.major = element_line(color = "#0f5733", linewidth = 0.2, linetype = 5),
    panel.grid.minor = ggplot2::element_blank(),
    
    # text
    text = element_text(color = "#0A3C23"),
    axis.text = element_text(color = "#0A3C23", family = "Archivo", size = 15),
    axis.title = element_text(color = "#0A3C23", family = "connectdisplay", size = 18),
    plot.title = element_text(color = "#0A3C23", , family = "connectdisplay", size = 28),
    plot.subtitle = element_text(color = "#0A3C23", family = "Archivo-Italic", size = 22),
    plot.caption = element_text(color = "#0A3C23", family = "Archivo", size = 12),
    
    # lines and borders
    axis.line = element_line(color = "#0A3C23"),
    axis.ticks = element_line(color = "#0A3C23"), 
    panel.border = element_rect(color = "#0A3C23", fill = NA),
    
    # legends
    legend.background = element_rect(color = "white"),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(color = "#0A3C23"),
    legend.title = element_text(color = "#0A3C23"),
    
    # facets
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color = "#0A3C23", face = "bold")
  )
}
