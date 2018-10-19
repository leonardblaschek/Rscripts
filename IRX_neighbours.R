library(reshape2)
library(ggplot2)
library(ggthemes)
library(showtext)
library(dplyr)
library(cowplot)
library(agricolae)

font_add(
  "Helvetica",
  regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
  italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
  bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf",
  bolditalic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-BdIt.otf"
)
showtext_auto()

irx.n <- read.csv("file:///home/leonard/Documents/Uni/PhD/IRX/IRX_neighbours.csv")

###############################
# functions for scaling and statistics
###############################
scale_this <- function(x){
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

tukey <- function(x) {
  aov1 <- aov(data = x, value ~ genotype)
  groups <- HSD.test(aov1, "genotype", alpha = 0.05)
  groups$groups[["genotype"]] <- rownames(groups$groups)
  groups$groups[["value"]] <- 0
  return(groups[["groups"]])
}

###############################
# create reference at the heigth of the cambium for distance calculation
###############################
irx.n.ref <-
  subset(irx.n, object == "ref", select = c(1, 2, 3, 8, 9, 29))
irx.n.ref$ref.y1 <- irx.n.ref$Y
irx.n.ref$ref.y2 <- irx.n.ref$Y
irx.n.ref$ref.x1 <- irx.n.ref$X - (irx.n.ref$Length / 2)
irx.n.ref$ref.x2 <- irx.n.ref$X + (irx.n.ref$Length / 2)

###############################
# merge reference into data frame
###############################
irx.n <-
  full_join(
    subset(irx.n, object != "ref"),
    irx.n.ref[, c(1, 2, 3, 7, 8, 9, 10)],
    by = c("genotype", "replicate", "technical")
  )

###############################
# calculate difference of the centre of each vessel to the cambium
###############################
irx.n$Distance <-
  apply(irx.n[, c("X", "Y", "ref.x1", "ref.x2", "ref.y1", "ref.y2")],
        1 ,
        function(x) {
          a <- c(x[1], x[2])
          b <- c(x[3], x[5])
          c <- c(x[4], x[6])
          v1 <- b - c
          v2 <- a - b
          m <- cbind(v1, v2)
          d <- abs(det(m)) / sqrt(sum(v1 * v1))
          d
        })

irx.n <- irx.n %>%
  group_by(genotype, replicate, number, object) %>%
  mutate(perim.P = sum(Length[adj.object=="P"]) / Perim.,
         perim.F = sum(Length[adj.object=="F"]) / Perim.,
         perim.PX = sum(Length[adj.object=="PX"]) / Perim.,
         perim.PMX = sum(Length[adj.object=="PMX"]) / Perim.,
         perim.SMX = sum(Length[adj.object=="SMX"]) / Perim.)

irx.n <- irx.n %>%  
  filter(., adj.object == "") %>%
  select(., c(1:5,12,15,16,20,34:39))

irx.lm <- irx.n %>% lm(Circ. ~ perim.F, data = .)

irx.n.melt <- melt(irx.n, id = c("genotype", "replicate", "technical", "object"))

###############################
# calculate z-scores
###############################
irx.n.melt <- irx.n.melt %>%
  group_by(variable) %>%
  mutate(value.scaled = scale_this(value))

###############################
# order from innermost to outermost cell type
###############################
irx.roi$object <-
  ordered(irx.roi$object, levels = c("PX", "PMX", "SMX"))

###############################
# remove number variable for overview plot
###############################
irx.n.melt <- subset(irx.n.melt, variable != "number")

###############################
# Tukey-HSD test 
###############################
# irx.n.letters <- irx.n.melt %>%
#   group_by(object, variable) %>%
#   do(data.frame(tukey(.)))

###############################
# shape and position overview plot
###############################
irx.overview <-
  ggplot(data = irx.n.melt, aes(x = genotype, y = value)) +
  geom_violin(draw_quantiles = 0.5, adjust = 1.5) +
  geom_jitter(
    aes(fill = value.scaled),
    shape = 21,
    width = 0.1,
    alpha = 0.9,
    size = 2,
    stroke = 0.25
  ) +
  # geom_label(data = irx.n.letters,
  #            aes(label = groups),
  #            fill = rgb(1, 1, 1, 0.75),
  #            hjust = 1,
  #            label.size = 0,
  #            family = "Helvetica") +
  scale_fill_distiller(palette = "RdBu", name = "Z-score by\ncolumn") +
  scale_y_continuous(expand = expand_scale(mult = c(0.21,0.05))) +
  # scale_x_discrete(
  #   labels = rev(c(
  #     "Col-0",
  #     # expression(italic("4cl1")),
  #     # expression(italic("4cl2")),
  #     # expression(paste(italic("4cl1"), "x", italic("4cl2"))),
  #     expression(italic("ccoaomt1")),
  #     expression(italic("fah1")),
  #     expression(italic("omt1")),
  #     expression(italic("ccr1")),
  #     # expression(paste(italic("ccr1"), "x", italic("fah1"))),
  #     # expression(italic("cad4")),
  #     # expression(italic("cad5")),
  #     expression(paste(italic("cad4"), "x", italic("cad5")))
  #   ))
  # ) +
  facet_grid(object ~ variable, scales = "free_x") +
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Helvetica"),
    strip.text = element_text(hjust = 0, face = "italic"),
    # axis.line.y = element_line(size = 0.75, lineend = "square"),
    axis.ticks = element_line(
      size = 0.25,
      lineend = "square",
      color = "black"
    ),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(
      colour = "black",
      size = 10,
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.25),
    panel.spacing = unit(1.5, "mm"),
    # plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.position = "bottom",
    # legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.height = unit(4, "mm"),
    legend.key.width = unit(30, "mm"),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  coord_flip()

pdf("irx_n_overview.pdf", width = 10, height = 10)
irx.overview
dev.off()

