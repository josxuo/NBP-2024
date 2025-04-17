#################
### BCS theme ###
#################

bcs_colors <- c(
  "dark green" = "#0A3C23",
  "cream" = "#FAF5F0",
  "yellow green" = "#E6FF55",
  "peach" = "#FFB98C",
  "bright green" = "#36BA3A"
)


theme_bcs <- function() {
  theme(
    # backgrounds
    panel.background = element_rect(fill = NA, color = NA),
    plot.background = element_rect(fill = NA, color = NA), 
    panel.grid.major = element_line(color = "grey", linewidth = 0.2, linetype = 1),
    panel.grid.minor = ggplot2::element_blank(),
    
    # text
    text = element_text(color = "#0A3C23"),
    axis.text = element_text(color = "#0A3C23", family = "Archivo", size = 15),
    axis.title = element_text(color = "#0A3C23", family = "connectdisplay", size = 18),
    plot.title = element_text(color = "#0A3C23", , family = "connectdisplay", 
                              hjust = 0.5, size = 28, lineheight = 0.4),
    plot.subtitle = element_text(color = "#0A3C23", family = "Archivo-Italic", size = 22, hjust = 0.5),
    plot.caption = element_text(color = "#0A3C23", family = "Archivo", size = 12),
    
    # lines and borders
    axis.line = element_line(color = "#0A3C23"),
    axis.ticks = element_line(color = "#0A3C23"), 
    panel.border = element_rect(color = "#0A3C23", fill = NA),
    
    # legends
    legend.position = "right",
    legend.justification = "top",
    legend.box = "vertical",
    legend.box.spacing = unit(0.5, "cm"),
    legend.background = element_rect(fill = NA, color = NA),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(color = "#0A3C23", family = "Archivo", size = 12),
    legend.title = element_text(color = "#0A3C23", family = "connectdisplay", size = 15),
    
    # facets
    strip.background = element_rect(fill = NA),
    strip.text = element_text(color = "#0A3C23", family = "connectdisplay", size = 15)
  )
}

chart_colors <- c(
"dark green" = "#0A3C23", #(Deep Forest Green - Brand Color)
"bright green" = "#36BA3A", #(Bright Green - Brand Color)
"peach" = "#FFB98C", # (Warm Peach - Brand Color)
#"yellow green" = "#E6FF55", # (Lime Yellow - Brand Color)
"evergreen" = "#195D37", # (Rich Evergreen - Darker complement of #0A3C23)
"mint green" = "#72D672", # (Soft Mint Green - Lighter complement of #36BA3A)
"coral" = "#FF8E54", # (Vibrant Coral - Bolder complement of #FFB98C)
#"electric yellow-green" = "#CCFF33", # (Electric Yellow-Green - Variation of #E6FF55)
"teal green" = "#478F5A", # (Muted Teal-Green - Mid-tone between #0A3C23 and #36BA3A)
"burnt orange" = "#D97B3E", #(Burnt Orange - Earthy counterbalance to #FFB98C)
"spring green" = "#B2E632",# (Fresh Spring Green - Brighter, fresher than #E6FF55)
"deep teal" = "#005F48", # (Deep Teal - Contrast against light greens/yellows)
"pastel green" = "#88D498",# (Pastel Green - Softer accent for balance)
"golden peach" = "#FFC477",# (Golden Peach - Harmonizing extension of #FFB98C)
"gold" = "#FFD700"# (Gold - Rich highlight color for emphasis)
)
