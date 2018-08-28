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
  read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol/foodweb.csv")
fw.data$od <- fw.data$ODx255 / 255

fw.avg <-
  ddply(fw.data,
        c("genotype", "cell.type", "adj.cell.type", "replicate"),
        summarise,
        od = mean(od, na.rm = TRUE))

fw.cast <-
  dcast(fw.avg, genotype + replicate ~ cell.type + adj.cell.type, mean, value.var = "od")
cast.avg <- ddply(
  fw.cast,
  c("genotype", "replicate"),
  summarise,
  "CB -> XF" = XF_CB / XF_XF,
  "MX -> XF" = XF_V / XF_XF,
  "IF -> XF" = XF_IF / XF_XF,
  "PA -> PX" = PX_PA / PX_PX,
  "MX -> PX" = PX_V / PX_PX,
  "LP -> IF" = IF_LP / IF_IF,
  "CB -> IF" = IF_CB / IF_IF,
  "XF -> IF" = IF_XF / IF_IF,
  "IF -> LP" = LP_IF/ LP_LP,
  "CB -> MX" = V_CB / V_V,
  "PA -> MX" = V_PA / V_V,
  "XF -> MX" = V_XF / V_V,
  "PX -> MX" = V_PX / V_V
)
cast.avg <- melt(cast.avg[, -2], id = c("genotype"))

ratio.avg <- ddply(cast.avg, c("genotype", "variable"), summarise,
                   mean = mean(value, na.rm = TRUE),
                   sd = sd(value, na.rm = TRUE))

fw.avg <-
  ddply(fw.data,
        c("genotype", "cell.type", "adj.cell.type"),
        summarise,
        od = mean(od, na.rm = TRUE))

fw.cast <-
  dcast(fw.avg, genotype ~ cell.type + adj.cell.type, mean, value.var = "od")
cast.avg <- ddply(
  fw.cast,
  c("genotype"),
  summarise,
  "CB -> XF" = XF_CB / XF_XF,
  "MX -> XF" = XF_V / XF_XF,
  "IF -> XF" = XF_IF / XF_XF,
  "PA -> PX" = PX_PA / PX_PX,
  "MX -> PX" = PX_V / PX_PX,
  "LP -> IF" = IF_LP / IF_IF,
  "CB -> IF" = IF_CB / IF_IF,
  "XF -> IF" = IF_XF / IF_IF,
  "IF -> LP" = LP_IF/ LP_LP,
  "CB -> MX" = V_CB / V_V,
  "PA -> MX" = V_PA / V_V,
  "XF -> MX" = V_XF / V_V,
  "PX -> MX" = V_PX / V_V
)
cast.avg <- melt(cast.avg, id = c("genotype"))
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

label_ccr <- data.frame(mean = 4.5, sd =  0, genotype = factor("ccr1", levels = c(
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

slopes <- data.frame(coef(summary(lm(fw.cast[, 12] ~ fw.cast[, 8])))["fw.cast[, 8]", "Estimate"], row.names = "slope")
colnames(slopes)[1] <- "PX -> MX"
slopes[, "MX -> PX"] <- coef(summary(lm(fw.cast[, 9] ~ fw.cast[, 13])))["fw.cast[, 13]", "Estimate"]
slopes[, "MX -> XF"] <- coef(summary(lm(fw.cast[, 17] ~ fw.cast[, 13])))["fw.cast[, 13]", "Estimate"]
slopes[, "XF -> MX"] <- coef(summary(lm(fw.cast[, 14] ~ fw.cast[, 18])))["fw.cast[, 18]", "Estimate"]
slopes[, "XF -> IF"] <- coef(summary(lm(fw.cast[, 4] ~ fw.cast[, 18])))["fw.cast[, 18]", "Estimate"]
slopes[, "IF -> XF"] <- coef(summary(lm(fw.cast[, 16] ~ fw.cast[, 2])))["fw.cast[, 2]", "Estimate"]
slopes[, "IF -> LP"] <- coef(summary(lm(fw.cast[, 5] ~ fw.cast[, 2])))["fw.cast[, 2]", "Estimate"]
slopes[, "LP -> IF"] <- coef(summary(lm(fw.cast[, 3] ~ fw.cast[, 6])))["fw.cast[, 6]", "Estimate"]
slopes[, "CB -> XF"] <- NA
slopes[, "PA -> PX"] <- NA
slopes[, "CB -> IF"] <- NA
slopes[, "CB -> MX"] <- NA
slopes[, "PA -> MX"] <- NA

intercepts <- data.frame(coef(summary(lm(fw.cast[, 12] ~ fw.cast[, 8])))["(Intercept)", "Estimate"], row.names = "intercept")
colnames(intercepts)[1] <- "PX -> MX"
intercepts[, "MX -> PX"] <- coef(summary(lm(fw.cast[, 9] ~ fw.cast[, 13])))["(Intercept)", "Estimate"]
intercepts[, "MX -> XF"] <- coef(summary(lm(fw.cast[, 17] ~ fw.cast[, 13])))["(Intercept)", "Estimate"]
intercepts[, "XF -> MX"] <- coef(summary(lm(fw.cast[, 14] ~ fw.cast[, 18])))["(Intercept)", "Estimate"]
intercepts[, "XF -> IF"] <- coef(summary(lm(fw.cast[, 4] ~ fw.cast[, 18])))["(Intercept)", "Estimate"]
intercepts[, "IF -> XF"] <- coef(summary(lm(fw.cast[, 16] ~ fw.cast[, 2])))["(Intercept)", "Estimate"]
intercepts[, "IF -> LP"] <- coef(summary(lm(fw.cast[, 5] ~ fw.cast[, 2])))["(Intercept)", "Estimate"]
intercepts[, "LP -> IF"] <- coef(summary(lm(fw.cast[, 3] ~ fw.cast[, 6])))["(Intercept)", "Estimate"]
intercepts[, "CB -> XF"] <- NA
intercepts[, "PA -> PX"] <- NA
intercepts[, "CB -> IF"] <- NA
intercepts[, "CB -> MX"] <- NA
intercepts[, "PA -> MX"] <- NA

r.squared <- data.frame(summary(lm(fw.cast[, 12] ~ fw.cast[, 8]))["adj.r.squared"], row.names = "r.squared")
colnames(r.squared)[1] <- "PX -> MX"
r.squared[, "MX -> PX"] <- summary(lm(fw.cast[, 9] ~ fw.cast[, 13]))["adj.r.squared"]
r.squared[, "MX -> XF"] <- summary(lm(fw.cast[, 17] ~ fw.cast[, 13]))["adj.r.squared"]
r.squared[, "XF -> MX"] <- summary(lm(fw.cast[, 14] ~ fw.cast[, 18]))["adj.r.squared"]
r.squared[, "XF -> IF"] <- summary(lm(fw.cast[, 4] ~ fw.cast[, 18]))["adj.r.squared"]
r.squared[, "IF -> XF"] <- summary(lm(fw.cast[, 16] ~ fw.cast[, 2]))["adj.r.squared"]
r.squared[, "IF -> LP"] <- summary(lm(fw.cast[, 5] ~ fw.cast[, 2]))["adj.r.squared"]
r.squared[, "LP -> IF"] <- summary(lm(fw.cast[, 3] ~ fw.cast[, 6]))["adj.r.squared"]
r.squared[, "CB -> XF"] <- NA
r.squared[, "PA -> PX"] <- NA
r.squared[, "CB -> IF"] <- NA
r.squared[, "CB -> MX"] <- NA
r.squared[, "PA -> MX"] <- NA

linregs <- rbind(slopes, intercepts, r.squared)
linregs$coef <- rownames(linregs)
linregs <- melt(linregs, id = "coef")
linregs <- dcast(linregs, variable ~ coef)
linregs$genotype <- "col-0"
linregs$variable <- as.character(linregs$variable)
linregs$variable <- as.factor(as.character(linregs$variable, levels = levels(ratio.avg$varaible)))

ratio.avg <- merge(ratio.avg, linregs, by = c("variable", "genotype"), all = TRUE)

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
  scale_y_continuous(limits = c(-0.5, 4.5)) +
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
  geom_text(data = label_ccr, label = paste("17", sprintf('\u2192'), sep = ""), size = 3, hjust = 1, colour = "#b2182b") +
  # geom_text(aes(x = genotype, y = - 7), label = 
  #   ifelse(is.na(ratio.avg$slope), 
  #          "", 
  #          ifelse(ratio.avg$intercept > 0, paste("y=", round(ratio.avg$slope, 1), "\n+", round(ratio.avg$intercept, 2), "\nR²=", round(ratio.avg$r.squared, 2), sep = ""),
  #                 paste("y=", round(ratio.avg$slope, 1), "\n", round(ratio.avg$intercept, 2), "\nR²=", round(ratio.avg$r.squared, 2), sep = ""))
  #          ), 
  # family = "Helvetica", size = 3, hjust = 0, vjust = 1) +
  facet_wrap( ~ variable, nrow = 1) +
  coord_flip()

linregs.col <- ggplot(data = linregs, aes(x = variable, y = slope, fill = slope)) +
  geom_point(size = 5, shape = 21) +
  scale_fill_viridis_c()

pdf("corrplots.pdf", width = 10, height = 10)
corrplot.mixed(fw.corr$r,
               tl.col = "black")
plot(fw.cast[, 17] ~ fw.cast[, 13], xlim = c(0,0.5), ylim = c(0, 0.5), ylab = "XF_MX", xlab = "MX_MX", main = "MX -> XF")
abline(lm(fw.cast[, 17] ~ fw.cast[, 13]))
plot(fw.cast[, 16] ~ fw.cast[, 2], xlim = c(0,0.5), ylim = c(0, 0.5), ylab = "XF_IF", xlab = "IF_IF", main = "IF -> XF")
abline(lm(fw.cast[, 16] ~ fw.cast[, 2]))
dev.off()

pdf("foodweb.pdf", width = 13, height = 2)
interaction.avg
linregs.col

dev.off()

write.csv(linregs[-c(9:13), -5], file = "regressions.csv", row.names = FALSE)
