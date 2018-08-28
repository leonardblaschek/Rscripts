library(sysfonts)
library(showtext)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(reshape2)
library(colorspace)
library(agricolae)
library(plyr)
library(dplyr)
library(rowr)
library(cowplot)
library(Hmisc)
library(corrplot)
library(PerformanceAnalytics)
library(colorspace)

font_add("Helvetica", regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf",
         bolditalic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-BdIt.otf")
showtext_auto()


fw.data <-
  read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol paper/foodweb.csv")
fw.data$od <- fw.data$ODx255 / 255

fw.avg <-
  ddply(fw.data,
        c("genotype", "cell.type", "adj.cell.type"),
        summarise,
        od = mean(od, na.rm = TRUE))

fw.cast <-
  dcast(subset(fw.avg), genotype ~ cell.type + adj.cell.type, mean, value.var = "od")
cast.avg <- ddply(
  fw.cast,
  c("genotype"),
  summarise,
  "CB -> XF" = XF_CB / XF_XF,
  "V -> XF" = XF_V / XF_XF,
  "IF -> XF" = XF_IF / XF_XF,
  "PA -> PX" = PX_PA / PX_PX,
  "V -> PX" = PX_V / PX_PX,
  "LP -> IF" = IF_LP / IF_IF,
  "CB -> IF" = IF_CB / IF_IF,
  "XF -> IF" = IF_XF / IF_IF,
  "IF -> LP" = LP_IF/ LP_LP,
  "CB -> V" = V_CB / V_V,
  "PA -> V" = V_PA / V_V,
  "XF -> V" = V_XF / V_V,
  "PX -> V" = V_PX / V_V
)
cast.avg <- melt(cast.avg, id = c("genotype"))

ratio.avg <- ddply(cast.avg, c("genotype", "variable"), summarise,
                   mean = mean(value, na.rm = TRUE),
                   sd = sd(value, na.rm = TRUE))
fw.cast <- data.matrix(fw.cast[, -1])
fw.corr <- rcorr(fw.cast)

cast.avg$genotype <-
  ordered(
    cast.avg$genotype,
    levels = c(
      "col-0",
      "4cl1",
      "4cl2",
      "4cl1x4cl2",
      "ccoaomt1",
      "fah1",
      "omt1",
      "ccr1",
      "cad4",
      "cad5",
      "cad4xcad5"
    )
  )

ratio.avg$genotype <-
  ordered(
    ratio.avg$genotype,
    levels = rev(c(
      "col-0",
      "4cl1",
      "4cl2",
      "4cl1x4cl2",
      "ccoaomt1",
      "fah1",
      "omt1",
      "ccr1",
      "cad4",
      "cad5",
      "cad4xcad5"
    ))
  )

label_ccr <- data.frame(mean = 4, sd =  0, genotype = factor("ccr1", levels = c(
  "col-0",
  "4cl1",
  "4cl2",
  "4cl1x4cl2",
  "ccoaomt1",
  "fah1",
  "omt1",
  "ccr1",
  "cad4",
  "cad5",
  "cad4xcad5"
)), variable = "XF -> IF")

fw.cast <- melt(data.frame(fw.cast))
regressions <- lm(as.matrix(fw.cast[, 2]) ~ fw.cast$variable)

interaction.avg <-
  ggplot(data = ratio.avg, aes(x = genotype, y = mean, fill = mean, ymin = mean - sd, ymax = mean + sd)) + 
  geom_hline(yintercept = 1, linetype = 1, size = 0.5) +
  # geom_line(group = 1) +
  geom_errorbar(width = 0.2) +
  geom_point(size = 3, shape = 21) + 
  # scale_fill_distiller(palette = "RdBu", limits = c(0.4, 1.6), direction = -1, na.value = "#b2182b") +
  scale_fill_gradientn(colours = c("#2166ac", "#67a9cf", "#d1e5f0", "white", "#fddbc7", "#ef8a62", "#b2182b"), limits = c(0.4,1.6),na.value = "#b2182b") +
  scale_x_discrete(
    labels = rev(c(
      "col-0",
      expression(italic("4cl1")),
      expression(italic("4cl2")),
      expression(italic("4cl1x4cl2")),
      expression(italic("ccoaomt1")),
      expression(italic("fah1")),
      expression(italic("omt1")),
      expression(italic("ccr1")),
      expression(italic("cad4")),
      expression(italic("cad5")),
      expression(italic("cad4xcad5"))
    ))
  ) +
  scale_y_continuous(limits = c(-0.5, 5)) +
  theme_minimal() + 
  theme(text = element_text(size = 14, family = "Helvetica"),
        axis.ticks = element_line(
          size = 0.25,
          lineend = "square",
          color = "black"
        ),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", size = 0.25),
        panel.spacing.x = unit(1.5, "mm"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        # legend.background = element_rect(fill = "grey95", color = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        strip.text = element_text(
          vjust = 0.1,
          hjust = 0,
          face = "italic"
        )
  ) +
  geom_text(data = label_ccr, label = paste("17", sprintf('\u2192'), sep = ""), size = 3) +
  facet_wrap( ~ variable, ncol = 13) +
  coord_flip()

pdf("corrplots.pdf", width = 10, height = 10)
corrplot.mixed(fw.corr$r,
               tl.col = "black")
# chart.Correlation(fw.cast, method = "pearson")
dev.off()

pdf("foodweb.pdf", width = 11, height = 2)
interaction.avg
dev.off()