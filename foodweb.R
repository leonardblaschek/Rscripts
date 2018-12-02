###############################
# If anyone is ever trying to understend this hot mess: I am sorry.
# I hope I have time to clean this up soon.
###############################

library(showtext)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(Hmisc)
library(corrplot)
library(PerformanceAnalytics)
library(forcats)

font_add("Helvetica", regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf",
         bolditalic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-BdIt.otf")
showtext_auto()

fw.data <-
  read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol/foodweb.csv")
fw.data <- mutate(fw.data, cell.type = recode(cell.type, "V" = "MX"),
                  adj.cell.type = recode(adj.cell.type, "V" = "MX"))


###############################
# Pixel values of the measured image were multiplied by 255 for better visualisation
###############################
fw.data$od <- fw.data$ODx255 / 255

fw.avg <-
  ddply(fw.data,
        c("genotype", "cell.type", "adj.cell.type", "replicate"),
        summarise,
        od = mean(od, na.rm = TRUE))

fw.cast <-
  dcast(fw.avg, genotype + replicate ~ cell.type + adj.cell.type, mean, value.var = "od")

###############################
# These are the values for fig. SX, showing the impact of individal muations on cell cooperativity.
# 
# Expressing how the neighbouring cell wall of cell type B influences the cell wall of cell type A,
# compared to cell walls of cell type A next to cell type A.
# 
# Averaged for each REPLICATE AND GENOTYPE.
###############################
cast.avg <- ddply(
  fw.cast,
  c("genotype", "replicate"),
  summarise,
  "CB -> XF" = XF_CB / XF_XF,
  "MX -> XF" = XF_MX / XF_XF,
  "IF -> XF" = XF_IF / XF_XF,
  "PA -> PX" = PX_PA / PX_PX,
  "MX -> PX" = PX_MX / PX_PX,
  "LP -> IF" = IF_LP / IF_IF,
  "CB -> IF" = IF_CB / IF_IF,
  "XF -> IF" = IF_XF / IF_IF,
  "IF -> LP" = LP_IF/ LP_LP,
  "CB -> MX" = MX_CB / MX_MX,
  "PA -> MX" = MX_PA / MX_MX,
  "XF -> MX" = MX_XF / MX_MX,
  "PX -> MX" = MX_PX / MX_MX
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

###############################
# These are the values for fig. SX, showing the impact of individal muations on cell cooperativity.
# 
# Expressing how the neighbouring cell wall of cell type B influences the cell wall of cell type A,
# compared to cell walls of cell type A next to cell type A.
# 
# Averaged for each GENOTYPE
###############################
cast.avg <- ddply(
  fw.cast,
  c("genotype"),
  summarise,
  "CB -> XF" = XF_CB / XF_XF,
  "MX -> XF" = XF_MX / XF_XF,
  "IF -> XF" = XF_IF / XF_XF,
  "PA -> PX" = PX_PA / PX_PX,
  "MX -> PX" = PX_MX / PX_PX,
  "LP -> IF" = IF_LP / IF_IF,
  "CB -> IF" = IF_CB / IF_IF,
  "XF -> IF" = IF_XF / IF_IF,
  "IF -> LP" = LP_IF/ LP_LP,
  "CB -> MX" = MX_CB / MX_MX,
  "PA -> MX" = MX_PA / MX_MX,
  "XF -> MX" = MX_XF / MX_MX,
  "PX -> MX" = MX_PX / MX_MX
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

###############################
# creating the geom_text dataframe for the annotation of the ccr1 in XF -> IF
###############################
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

###############################
# Summarising the correlation parameters in a useable table.
# 
# The correlations express how the absorbance of a cell wall of cell type A next to cell wall of cell type B
# correlates with the absorbance in cell walls of cell type B next to cell type B
###############################
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

r <- data.frame(fw.corr$r["PX_PX", "MX_PX"], row.names = "r")
colnames(r)[1] <- "PX -> MX"
r[, "MX -> PX"] <- fw.corr$r["MX_MX", "PX_MX"]
r[, "MX -> XF"] <- fw.corr$r["MX_MX", "XF_MX"]
r[, "XF -> MX"] <- fw.corr$r["XF_XF", "MX_XF"]
r[, "XF -> IF"] <- fw.corr$r["XF_XF", "IF_XF"]
r[, "IF -> XF"] <- fw.corr$r["IF_IF", "XF_IF"]
r[, "IF -> LP"] <- fw.corr$r["IF_IF", "LP_IF"]
r[, "LP -> IF"] <- fw.corr$r["LP_LP", "IF_LP"]
r[, "CB -> XF"] <- NA
r[, "PA -> PX"] <- NA
r[, "CB -> IF"] <- NA
r[, "CB -> MX"] <- NA
r[, "PA -> MX"] <- NA

P <- data.frame(fw.corr$P["PX_PX", "MX_PX"], row.names = "P")
colnames(P)[1] <- "PX -> MX"
P[, "MX -> PX"] <- fw.corr$P["MX_MX", "PX_MX"]
P[, "MX -> XF"] <- fw.corr$P["MX_MX", "XF_MX"]
P[, "XF -> MX"] <- fw.corr$P["XF_XF", "MX_XF"]
P[, "XF -> IF"] <- fw.corr$P["XF_XF", "IF_XF"]
P[, "IF -> XF"] <- fw.corr$P["IF_IF", "XF_IF"]
P[, "IF -> LP"] <- fw.corr$P["IF_IF", "LP_IF"]
P[, "LP -> IF"] <- fw.corr$P["LP_LP", "IF_LP"]
P[, "CB -> XF"] <- NA
P[, "PA -> PX"] <- NA
P[, "CB -> IF"] <- NA
P[, "CB -> MX"] <- NA
P[, "PA -> MX"] <- NA

linregs <- rbind(slopes, intercepts, r.squared, r, P)
linregs$coef <- rownames(linregs)
linregs <- melt(linregs, id = "coef")
linregs <- dcast(linregs, variable ~ coef)
linregs$genotype <- "col-0"
linregs$variable <- as.character(linregs$variable)
linregs$variable <- as.factor(as.character(linregs$variable, levels = levels(ratio.avg$variable)))
ratio.avg <- merge(ratio.avg, linregs, by = c("variable", "genotype"), all = TRUE)

###############################
# table of regression (r.squared, slope, intercept) and correlation (r, P) parameters
###############################)
write.csv(linregs[-c(9:13), -7], file = "regressions.csv", row.names = FALSE)

###############################
# plot fig. SX
###############################
interaction.avg <-
  ggplot(data = ratio.avg, aes(
    x = genotype, 
    y = mean, 
    fill = mean, 
    ymin = mean - sd, 
    ymax = mean + sd
  )) + 
  geom_hline(
    yintercept = 1, 
    linetype = 1, 
    size = 0.5
  ) +
  geom_errorbar(width = 0.2) +
  geom_point(size = 3, shape = 21) + 
  scale_fill_gradientn(
    colours = c("#2166ac", "#67a9cf", "#d1e5f0", "white", "#fddbc7", "#ef8a62", "#b2182b"), 
    limits = c(0.4,1.6),
    na.value = "#b2182b"
  ) +
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
  labs(y = "Absorbance ratio") +
  theme(
    text = element_text(size = 14, family = "Helvetica"),
    axis.ticks = element_line(
      size = 0.25,
      lineend = "square",
      color = "black"
    ),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      colour = "black"
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      fill = NA, 
      color = "black", 
      size = 0.25
    ),
    panel.spacing.x = unit(1.5, "mm"),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    strip.text = element_text(
      vjust = 0.1,
      hjust = 0,
      face = "italic"
    )
  ) +
  geom_text(data = label_ccr, 
            label = paste("17", sprintf('\u2192'), sep = ""), 
            size = 3, 
            hjust = 1, 
            colour = "#b2182b") +
  facet_wrap( ~ variable, nrow = 1) +
  coord_flip()

###############################
# plot overall correlation matrix
###############################
pdf("corrplots.pdf", width = 10, height = 10)
# corrplot.mixed(fw.corr$r,
#                tl.col = "black")
chart.Correlation(fw.cast, histogram=TRUE, pch=21)
dev.off()

pdf("foodweb.pdf", width = 13, height = 2)
interaction.avg

dev.off()



