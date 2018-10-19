library(reshape2)
library(ggplot2)
library(ggthemes)
library(showtext)
library(dplyr)
library(cowplot)

font_add(
  "Helvetica",
  regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
  italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
  bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf",
  bolditalic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-BdIt.otf"
)
showtext_auto()


irx.plot <- function(x) {
  
  if (!is.null(x)) {
    irx.plot.data <- subset(irx.roi, genotype == x)
  } else {
    irx.plot.data <- irx.roi
  }
  
  if (!is.null(x)) {
    irx.polygon <-
      ggplot(data = irx.plot.data,
             aes(
               x = ROIx,
               y = ROIy,
               group = interaction(object, number)
             )) +
      geom_polygon(
        colour = "black",
        aes(fill = Perim., linetype = object),
        alpha = 0.5,
        size = 0.5
      ) +
      scale_fill_viridis_c(name = "Perimeter [µm]", limits = c(0, 125)) +
      guides(linetype = guide_legend(title = "")) +
      xlim(0, 1.05) +
      ylim(0, 1.05) +
      labs(x = "Circularity",
           y = "Relative Distance") +
      annotate("text", x = 0.025, y = 1, label = x, family = "Helvetica", size = 5) +
      theme(
        text = element_text(family = "Helvetica"),
        legend.position = c(0.025, 0.75),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
      )
  } else {
    irx.polygon <-
      ggplot(data = irx.plot.data,
             aes(
               x = ROIx,
               y = ROIy,
               group = interaction(object, number)
             )) +
      geom_polygon(
        colour = "black",
        aes(fill = genotype, linetype = object),
        alpha = 0.5,
        size = 0.5
      ) +
      scale_fill_few() +
      guides(linetype = guide_legend(order = 1)) +
      xlim(0, 1.05) +
      ylim(0, 1.05) +
      labs(x = "Circularity",
           y = "Relative Distance") +
      theme(
        text = element_text(family = "Helvetica"),
        legend.position = c(0.025, 0.8),
        legend.title = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
      )
  }
  
  irx.polygon.densx <-
    ggplot(data = irx.plot.data,
           aes(Circ., linetype = object)) +
    geom_density(adjust = 3) +
    xlim(0, 1.05) +
    # ylim(0, 1.05) +
    # theme_void() +
    theme(
      legend.position = "none",
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank()
    )
  
  irx.polygon.densy <-
    ggplot(data = irx.plot.data,
           aes(rel.distance, linetype = object)) +
    geom_density(adjust = 3) +
    xlim(0, 1.05) +
    # ylim(0, 1.05) +
    # theme_void() +
    theme(
      legend.position = "none",
      plot.margi = unit(c(0, 0, 1.25, 0), "cm"),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank()
    ) +
    coord_flip()
  
  first.col <-
    plot_grid(
      irx.polygon,
      irx.polygon.densx,
      ncol = 1,
      rel_heights = c(10, 1),
      align = "v"
    )
  second.col <-
    plot_grid(
      irx.polygon.densy,
      NULL,
      ncol = 1,
      rel_heights = c(10, 1),
      align = "h"
    )
  scatter.grid <-
    plot_grid(
      first.col,
      second.col,
      ncol = 2,
      rel_widths = c(10, 1),
      align = "hv"
    )
  scatter.grid
}

pdf("irx.pdf", width = 10, height = 10)
irx.plot("col-0")
irx.plot("ccr1-3")
irx.plot("ccoaomt1")
irx.plot("cad4xcad5")
irx.plot("fah1")
irx.plot("omt1")
irx.plot("ccr1xfah1")
irx.plot("4cl1")
irx.plot(NULL)
dev.off()