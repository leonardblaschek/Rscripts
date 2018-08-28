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


poplar <- read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol/poplar_foodweb.csv")
poplar$replicate <- factor(poplar$replicate)
poplar$technical <- factor(poplar$technical)
poplar$number <- row(poplar)
poplar$cell.type <- plyr::revalue(poplar$cell.type, c(
  "F" = "Fibre",
  "V" = "Vessel",
  "R" = "Ray",
  "CB" = "Cambium"
))
poplar$adj.cell.type <- plyr::revalue(poplar$adj.cell.type , c(
  "F" = "Fibre",
  "V" = "Vessel",
  "R" = "Ray",
  "CB" = "Cambium",
  "PA" = "Parenchyma"
))

poplar$Distance <-
  apply(poplar[, c("X", "Y", "ref.x1", "ref.x2", "ref.y1", "ref.y2")],
        1 ,
        function(x) {
          a <- c(x[1], x[2])
          b <- c((x[3] / 5.9), (x[5] / 5.9))
          c <- c((x[4] / 5.9), (x[6] / 5.9))
          v1 <- b - c
          v2 <- a - b
          m <- cbind(v1, v2)
          d <- abs(det(m)) / sqrt(sum(v1 * v1))
          d
        })

poplar$diff <- poplar$OD.stained - poplar$OD.unstained
poplar.bg <- subset(poplar, cell.type == "Cambium" & adj.cell.type == "Cambium",
                    select = c("genotype", "replicate", "technical", "diff"))
colnames(poplar.bg)[4] <- "OD.bg"

poplar.bg <- poplar.bg %>%
  group_by(genotype, replicate, technical) %>%
  mutate(OD.bg = mean(OD.bg, na.rm = TRUE))

poplar <- merge(poplar, unique(poplar.bg), by =c("genotype", "replicate", "technical"))

poplar$diff.adj <- poplar$diff - poplar$OD.bg
poplar$genotype <- ordered(poplar$genotype, levels = c("WT", "c4h", "ccr"))

fw.pop.avg <-
  ddply(poplar,
        c("genotype", "cell.type", "adj.cell.type", "replicate"),
        summarise,
        od = mean(diff.adj, na.rm = TRUE))

fw.pop.cast <-
  dcast(fw.pop.avg, genotype + replicate ~ cell.type + adj.cell.type, mean, value.var = "od")
cast.pop.avg <- ddply(
  fw.pop.cast,
  c("genotype", "replicate"),
  summarise,
  "CB -> V" = Vessel_Cambium / Vessel_Vessel,
  "R -> V" = Vessel_Ray / Vessel_Vessel,
  "F -> V" = Vessel_Fibre / Vessel_Vessel,
  "PA -> V" = Vessel_Parenchyma / Vessel_Vessel,
  "CB -> F" = Fibre_Cambium / Fibre_Fibre,
  "V -> F" = Fibre_Vessel / Fibre_Fibre,
  "R -> F" = Fibre_Ray / Fibre_Fibre,
  "PA -> F" = Fibre_Parenchyma / Fibre_Fibre,
  "CB -> R" = Ray_Cambium/ Ray_Ray,
  "V -> R" = Ray_Vessel / Ray_Ray,
  "F -> R" = Ray_Fibre / Ray_Ray
)
cast.pop.avg <- melt(cast.pop.avg[, -2], id = c("genotype"))

ratio.avg <- ddply(cast.pop.avg, c("genotype", "variable"), summarise,
                   mean = mean(value, na.rm = TRUE),
                   sd = sd(value, na.rm = TRUE))

fw.pop.avg <-
  ddply(poplar,
        c("genotype", "cell.type", "adj.cell.type"),
        summarise,
        od = mean(diff.adj, na.rm = TRUE))

fw.pop.cast <-
  dcast(fw.pop.avg, genotype ~ cell.type + adj.cell.type, mean, value.var = "od")
cast.pop.avg <- ddply(
  fw.pop.cast,
  c("genotype"),
  summarise,
  "CB -> V" = Vessel_Cambium / Vessel_Vessel,
  "R -> V" = Vessel_Ray / Vessel_Vessel,
  "F -> V" = Vessel_Fibre / Vessel_Vessel,
  "PA -> V" = Vessel_Parenchyma / Vessel_Vessel,
  "CB -> F" = Fibre_Cambium / Fibre_Fibre,
  "V -> F" = Fibre_Vessel / Fibre_Fibre,
  "R -> F" = Fibre_Ray / Fibre_Fibre,
  "PA -> F" = Fibre_Parenchyma / Fibre_Fibre,
  "CB -> R" = Ray_Cambium/ Ray_Ray,
  "V -> R" = Ray_Vessel / Ray_Ray,
  "F -> R" = Ray_Fibre / Ray_Ray
)
cast.pop.avg <- melt(cast.pop.avg, id = c("genotype"))
fw.pop.cast <- data.matrix(fw.pop.cast[, -1])
# fw.corr <- rcorr(fw.pop.cast)

cast.pop.avg$genotype <-
  ordered(
    cast.pop.avg$genotype,
    levels = c(
      "WT",
      "c4h",
      "ccr"
    )
  )

ratio.avg$genotype <-
  ordered(
    ratio.avg$genotype,
    levels = rev(c(
      "WT",
      "c4h",
      "ccr"
    ))
  )


slopes <- data.frame(coef(summary(lm(fw.pop.cast[, 12] ~ fw.pop.cast[, 3])))["fw.pop.cast[, 3]", "Estimate"], row.names = "slope")
colnames(slopes)[1] <- "F -> V"
slopes[, "V -> R"] <- coef(summary(lm(fw.pop.cast[, 10] ~ fw.pop.cast[, 15])))["fw.pop.cast[, 15]", "Estimate"]
slopes[, "V -> F"] <- coef(summary(lm(fw.pop.cast[, 6] ~ fw.pop.cast[, 15])))["fw.pop.cast[, 15]", "Estimate"]
slopes[, "R -> V"] <- coef(summary(lm(fw.pop.cast[, 14] ~ fw.pop.cast[, 9])))["fw.pop.cast[, 9]", "Estimate"]
slopes[, "R -> F"] <- coef(summary(lm(fw.pop.cast[, 5] ~ fw.pop.cast[, 9])))["fw.pop.cast[, 9]", "Estimate"]
slopes[, "F -> R"] <- coef(summary(lm(fw.pop.cast[, 8] ~ fw.pop.cast[, 3])))["fw.pop.cast[, 3]", "Estimate"]
slopes[, "CB -> R"] <- NA
slopes[, "CB -> V"] <- NA
slopes[, "CB -> F"] <- NA
slopes[, "PA -> V"] <- NA
slopes[, "PA -> F"] <- NA

intercepts <- data.frame(coef(summary(lm(fw.pop.cast[, 12] ~ fw.pop.cast[, 3])))["(Intercept)", "Estimate"], row.names = "intercept")
colnames(intercepts)[1] <- "F -> V"
intercepts[, "V -> R"] <- coef(summary(lm(fw.pop.cast[, 10] ~ fw.pop.cast[, 15])))["(Intercept)", "Estimate"]
intercepts[, "V -> F"] <- coef(summary(lm(fw.pop.cast[, 6] ~ fw.pop.cast[, 15])))["(Intercept)", "Estimate"]
intercepts[, "R -> V"] <- coef(summary(lm(fw.pop.cast[, 14] ~ fw.pop.cast[, 9])))["(Intercept)", "Estimate"]
intercepts[, "R -> F"] <- coef(summary(lm(fw.pop.cast[, 5] ~ fw.pop.cast[, 9])))["(Intercept)", "Estimate"]
intercepts[, "F -> R"] <- coef(summary(lm(fw.pop.cast[, 8] ~ fw.pop.cast[, 3])))["(Intercept)", "Estimate"]
intercepts[, "CB -> R"] <- NA
intercepts[, "CB -> V"] <- NA
intercepts[, "CB -> F"] <- NA
intercepts[, "PA -> V"] <- NA
intercepts[, "PA -> F"] <- NA

r.squared <- data.frame(summary(lm(fw.pop.cast[, 12] ~ fw.pop.cast[, 3]))["adj.r.squared"], row.names = "r.squared")
colnames(r.squared)[1] <- "F -> V"
r.squared[, "V -> R"] <- summary(lm(fw.pop.cast[, 10] ~ fw.pop.cast[, 15]))["adj.r.squared"]
r.squared[, "V -> F"] <- summary(lm(fw.pop.cast[, 6] ~ fw.pop.cast[, 15]))["adj.r.squared"]
r.squared[, "R -> V"] <- summary(lm(fw.pop.cast[, 14] ~ fw.pop.cast[, 9]))["adj.r.squared"]
r.squared[, "R -> F"] <- summary(lm(fw.pop.cast[, 5] ~ fw.pop.cast[, 9]))["adj.r.squared"]
r.squared[, "F -> R"] <- summary(lm(fw.pop.cast[, 8] ~ fw.pop.cast[, 3]))["adj.r.squared"]
r.squared[, "CB -> R"] <- NA
r.squared[, "CB -> V"] <- NA
r.squared[, "CB -> F"] <- NA
r.squared[, "PA -> V"] <- NA
r.squared[, "PA -> F"] <- NA

linregs <- rbind(slopes, intercepts, r.squared)
linregs$coef <- rownames(linregs)
linregs <- melt(linregs, id = "coef")
linregs <- dcast(linregs, variable ~ coef)
linregs$genotype <- "WT"
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
  scale_fill_gradientn(colours = c("#2166ac", "#67a9cf", "#d1e5f0", "white", "#fddbc7", "#ef8a62", "#b2182b"), limits = c(-1,3),na.value = "#b2182b") +
  # scale_x_discrete(
  #   labels = rev(c(
  #     "WT",
  #     expression(italic("c4h")),
  #     expression(italic("ccr")),
  #     expression(italic("4cl1x4cl2")),
  #     expression(italic("ccoaomt1")),
  #     expression(italic("fah1")),
  #     expression(italic("omt1")),
  #     expression(italic("ccr1")),
  #     expression(italic("cad4")),
  #     expression(italic("cad5")),
  #     expression(italic("cad4xcad5"))
  #   ))
  # ) +
  scale_y_continuous(limits = c(-1, 8)) +
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
  # geom_text(aes(x = genotype, y = - 7), label = 
            #   ifelse(is.na(ratio.avg$slope), 
            #          "", 
            #          ifelse(ratio.avg$intercept > 0, paste("y=", round(ratio.avg$slope, 1), "\n+", round(ratio.avg$intercept, 2), "\nR²=", round(ratio.avg$r.squared, 2), sep = ""),
            #                 paste("y=", round(ratio.avg$slope, 1), "\n", round(ratio.avg$intercept, 2), "\nR²=", round(ratio.avg$r.squared, 2), sep = ""))
            #          ), 
            # family = "Helvetica", size = 3, hjust = 0, vjust = 1) +
  facet_wrap( ~ variable, nrow = 3) +
  coord_flip()

linregs.col <- ggplot(data = linregs, aes(x = variable, y = slope, fill = slope)) +
  geom_point(size = 5, shape = 21) +
  scale_fill_viridis_c()

# pdf("corrplots.pdf", width = 10, height = 10)
# corrplot.mixed(fw.corr$r,
#                tl.col = "black")
# plot(fw.pop.cast[, 17] ~ fw.pop.cast[, 13], xlim = c(0,0.5), ylim = c(0, 0.5), ylab = "XF_MX", xlab = "MX_MX", main = "MX -> XF")
# abline(lm(fw.pop.cast[, 17] ~ fw.pop.cast[, 13]))
# plot(fw.pop.cast[, 16] ~ fw.pop.cast[, 2], xlim = c(0,0.5), ylim = c(0, 0.5), ylab = "XF_IF", xlab = "IF_IF", main = "IF -> XF")
# abline(lm(fw.pop.cast[, 16] ~ fw.pop.cast[, 2]))
# dev.off()

pdf("foodweb_pop.pdf", width = 5, height = 5)
interaction.avg
plot(fw.pop.cast[, 6] ~ fw.pop.cast[,15], xlim = c(0,1), ylim = c(0,1))
linregs.col
dev.off()

write.csv(linregs[-c(9:13), -5], file = "regressions.csv", row.names = FALSE)
