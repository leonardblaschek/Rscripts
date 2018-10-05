library(reshape2)
library(ggplot2)
library(ggthemes)
library(showtext)
library(dplyr)
library(cowplot)
library(agricolae)

font_add("Helvetica", regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf",
         bolditalic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-BdIt.otf")
showtext_auto()

irx.data <-
  read.csv("file:///home/leonard/Documents/Uni/PhD/IRX/irx_measurements.csv", sep = "\t")
# irx.data <- subset(irx.data, replicate == 1 & technical == 3)
irx.ref <- subset(irx.data, object == "ref", select = c(1, 2, 3, 6, 7, 27))
irx.ref$ref.y1 <- irx.ref$Y
irx.ref$ref.y2 <- irx.ref$Y
irx.ref$ref.x1 <- irx.ref$X - (irx.ref$Length / 2)
irx.ref$ref.x2 <- irx.ref$X + (irx.ref$Length / 2)

irx.data <-
  merge(subset(irx.data, object != "ref"),
        irx.ref[, c(1, 2, 3, 7, 8, 9, 10)],
        by = c("genotype", "replicate", "technical"))

irx.data$Distance <-
  apply(irx.data[, c("X", "Y", "ref.x1", "ref.x2", "ref.y1", "ref.y2")],
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
irx.data$number <- row(irx.data)
irx.roi <- irx.data


irx.roix <- data.frame(irx.roi %>% select(genotype, replicate, technical, object, number, contains("ROIx")))
irx.roix <- melt(irx.roix, id = c("genotype", "replicate", "technical", "object", "number"), value.name = "ROIx")
irx.roix <- cbind(irx.roix, colsplit(irx.roix$variable, "x", c("variable", "point")))

irx.roiy <- data.frame(irx.roi %>% select(genotype, replicate, technical, object, number, contains("ROIy")))
irx.roiy <- melt(irx.roiy, id = c("genotype", "replicate", "technical", "object", "number"), value.name = "ROIy")
irx.roiy <- cbind(irx.roiy, colsplit(irx.roiy$variable, "x", c("variable", "point")))

irx.roi <- cbind(irx.roix[, c(1:5, 7, 9)], irx.roiy[, 7])
colnames(irx.roi)[8] <- "ROIy"
irx.roi <- irx.roi[apply(irx.roi[c(6,8)],1,function(z) !any(z==0)),]

irx.data <- irx.data[, c("genotype", "replicate", "technical", "object", "number", "Area", "Perim.", "Width", "Height", "Circ.", "Distance")]
irx.data$technical <- as.factor(as.character(irx.data$technical))
irx.data$object <- ordered(irx.data$object, levels = c("PX", "PMX", "SMX"))
irx.melt <- melt(irx.data, id = c("genotype", "replicate", "technical", "object"))

irx.data <- irx.data %>%
  group_by(genotype, replicate, technical) %>%
  mutate(Distance = Distance / max(Distance, na.rm = TRUE))

irx.roi <- merge(irx.roi, irx.data, by = c("genotype", "replicate", "technical", "object", "number"))
irx.roi <- irx.roi[order(irx.roi$genotype, irx.roi$replicate, irx.roi$technical, irx.roi$number, irx.roi$point), ]

irx.roi$geno.level <- as.numeric(irx.roi$genotype)
irx.data$geno.level <- as.numeric(irx.data$genotype)
irx.roi <- irx.roi %>%
  group_by(genotype, replicate, technical, object, number) %>%
  mutate(ROIy = (Circ. + ((ROIy - mean(ROIy, na.rm = TRUE))) / 2500),
         ROIx = (geno.level + ((ROIx - mean(ROIx, na.rm = TRUE))) / 1250))

irx.roi$object <- ordered(irx.roi$object, levels = c("PX", "PMX", "SMX"))

irx.avg <- irx.data %>%
  group_by(genotype, object) %>%
  mutate(ROIy = mean(Circ.))
colnames(irx.avg)[12] <- "ROIx"
irx.avg <- irx.avg[!duplicated(irx.avg$ROIy),]

file.remove("stats_circ.csv")

print.HSD.circ <- function(x) {
  aov1 <- aov(Circ. ~ genotype, data = x)
  groups <- HSD.test(aov1, "genotype", alpha = 0.05)
  groups[["groups"]][["object"]] <- unique(x$object)
  groups[["groups"]][["Circ."]] <- 0
  write.table(
    groups[["groups"]],
    file = "stats_circ.csv",
    append = TRUE,
    sep = ",",
    col.names = FALSE
  )
}

# set statistical letters for circularity
irx.data %>%
  group_by(object) %>%
  do(data.frame(print.HSD.circ(.)))
letters.circ <-
  read.csv("file:///home/leonard/R/Output/autopheno/stats_circ.csv",
           header = FALSE)
colnames(letters.circ) <-
  c("genotype", "ROIy", "group", "object")
letters.circ$ROIx <- as.numeric(letters.circ$genotype)
# letters.circ$cell <-
#   as.factor(paste(letters.circ$genotype, letters.circ$object))

irx.plot <- function() {
irx.polygon <- ggplot(data = irx.roi, aes(x = ROIx, y = ROIy, group = interaction(object, number))) +
  geom_polygon(colour = "black", aes(fill = Distance), alpha = 0.5, size = 0.1) +
  geom_point(data = irx.avg, aes(x = ROIx, y = ROIy), 
             fill = NA, group = NA, size = 15, shape = 45, colour = "red", alpha = 0.75) +
  # stat_summary(data = irx.avg, aes(x = ROIx, y = ROIy), fun.y = mean, fun.ymin = mean, fun.ymax = mean,
  #              geom = "crossbar", width = 1, size = 0.5, fatten = 0) +
  scale_fill_viridis_c() +
  guides(linetype = guide_legend(title="")) +
  scale_x_continuous(limits = c(0.9, 3.1), breaks = c(1, 2, 3), labels = rev(c("col-0", "ccr1-3", "ccoaomt1"))) +
  ylim(0, 1.05) +
  theme_base() +
  labs(x = "Genotype", y = "Circularity") +
  theme(text = element_text(family = "Helvetica"),
        legend.position = "right",
        plot.margin = unit(c(0,0,0,0), "cm"),
        plot.background = element_blank()) +
  geom_text(
    data = letters.circ,
    aes(label = group),
    group = NA,
    family = "Helvetica",
    angle = 0,
    colour = "black",
    hjust = 0.5,
    size = 4,
    position = position_dodge(width = 0.8)
  ) +
  facet_wrap(~ object)

# irx.polygon.densx <- ggplot(data = irx.data, aes(Circ., linetype = object)) +
#   geom_density() +
#   xlim(0, 1.05) +
#   # ylim(0, 1.05) +
#   # theme_void() +
#   theme(legend.position = "none",
#         plot.margin = unit(c(0,0,0,0), "cm"),
#         axis.title=element_blank(),
#         axis.text=element_blank(),
#         axis.line=element_blank(),
#         axis.ticks=element_blank())
# 
# irx.polygon.densy <- ggplot(data = irx.data, aes(Distance, linetype = object)) +
#   geom_density() +
#   xlim(0, 1.05) +
#   # ylim(0, 1.05) +
#   # theme_void() +
#   theme(legend.position = "none",
#         plot.margi = unit(c(0,0,1.25,0), "cm"),
#         axis.title=element_blank(),
#         axis.text=element_blank(),
#         axis.line=element_blank(),
#         axis.ticks=element_blank()) +
#   coord_flip()
# 
# first.col <- plot_grid(irx.polygon, irx.polygon.densx, ncol = 1, rel_heights = c(10, 1), align = "v")
# second.col <- plot_grid(irx.polygon.densy, NULL, ncol = 1, rel_heights = c(10, 1), align = "h")
# scatter.grid <- plot_grid(first.col, second.col, ncol = 2, rel_widths = c(10, 1), align = "hv")
# scatter.grid
irx.polygon
}

pdf("irx_facet.pdf", width = 12, height = 4)
irx.plot()
dev.off()
