library(ggplot2)
library(ggthemes)
library(plyr)
library(reshape2)
library(showtext)
library(cowplot)

font_add("Helvetica",
         regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf")
showtext_auto()

lm_eqn = function(m) {
  l <- list(
    a = format(coef(m)[1], digits = 2),
    b = format(abs(coef(m)[2]), digits = 2),
    r2 = format(summary(m)$r.squared, digits = 3)
  )
  
  
  if (coef(m)[2] >= 0)  {
    eq <-
      substitute(italic(y) == a + b ~ italic(x) * "," ~ italic(R) ^ 2 ~ "= " ~
                   r2, l)
  } else {
    eq <-
      substitute(italic(y) == a - b ~ italic(x) * "," ~ italic(R) ^ 2 ~ "= " ~
                   r2, l)
  }
  
  as.character(as.expression(eq))
  
}

staining <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol paper/17-10_timelines/video_reworked.csv"
  )
staining$time <- (staining$slice * 15) - 15
staining$OD <- ifelse(staining$OD > 1 , NA, staining$OD)
# staining$OD <- 255/(10^staining$OD)
staining <- subset(staining, time < 256)
staining$time <- as.factor(as.character(staining$time))
staining$thickness <-
  ordered(as.factor(as.character(staining$thickness)),
          levels = c("12", "25", "50", "100", "150"))

staining.time <- paste(sort(as.integer(levels(staining$time))))
staining$time <-
  ordered(as.factor(as.character(staining$time)), levels = staining.time)
staining[1:50 + rep(seq(0, ((nrow(
  staining
) - 50) - 50), by = 300), each = 50), 2] <- "IF"
staining[51:100 + rep(seq(0, (nrow(staining) - 50), by = 300), each = 50), 2] <-
  "MX"
staining[101:150 + rep(seq(0, (nrow(staining) - 50), by = 300), each = 50), 2] <-
  "XF"
staining[151:200 + rep(seq(0, (nrow(staining) - 50), by = 300), each = 50), 2] <-
  "LP"
staining[201:250 + rep(seq(0, (nrow(staining) - 50), by = 300), each = 50), 2] <-
  "PX"
staining[251:300 + rep(seq(0, (nrow(staining) - 50), by = 300), each = 50), 2] <-
  "PH"
staining$cell.type <-
  ordered(staining$cell.type, levels = c("PH", "IF", "LP", "MX", "PX", "XF"))
staining$hue <- ((staining$H + 128) / 255) * 360
staining$point <- NA
pointcount <- data.frame(1:50)
colnames(pointcount) <- 'point'
staining[1:50 + rep(seq(0, ((nrow(
  staining
) - 50)), by = 50), each = 50), 10] <- pointcount$point
staining.bg <- subset(staining, time == "0", select = 1:5)
staining.bg$OD.bg <- staining.bg$OD
staining.bg$point <- NA
staining.bg[1:50 + rep(seq(0, ((
  nrow(staining.bg) - 50
)), by = 50), each = 50), 7] <- pointcount$point
staining.bg$OD <- NULL
staining <-
  merge(
    staining,
    staining.bg,
    by = c("genotype", "replicate", "cell.type", "thickness", "point")
  )
staining$OD.adj <- round((staining$OD - staining$OD.bg), digits = 4)

staining.cml <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol paper/17-10_timelines/video_CML.csv"
  )
staining.cml$time <- (staining.cml$slice * 15) - 15
staining.cml$OD <- ifelse(staining.cml$OD > 1 , NA, staining.cml$OD)
# staining.cml$OD <- 255/(10^staining.cml$OD)
staining.cml <- subset(staining.cml, time < 256)
staining.cml$time <- as.factor(as.character(staining.cml$time))
staining.cml$thickness <-
  ordered(as.factor(as.character(staining.cml$thickness)),
          levels = c("12", "25", "50", "100", "150"))
staining.cml.time <-
  paste(sort(as.integer(levels(staining.cml$time))))
staining.cml$time <-
  ordered(as.factor(as.character(staining.cml$time)), levels = staining.cml.time)
staining.cml$hue <- ((staining.cml$H + 128) / 255) * 360

staining.cml$point <- NA
pointcount <- data.frame(1:50)
colnames(pointcount) <- 'point'
staining.cml[1:50 + rep(seq(0, ((
  nrow(staining.cml) - 50
)), by = 50), each = 50), 10] <- pointcount$point
staining.cml.bg <- subset(staining.cml, time == "0", select = 1:5)
staining.cml.bg$OD.bg <- staining.cml.bg$OD
staining.cml.bg$point <- NA
staining.cml.bg[1:50 + rep(seq(0, ((
  nrow(staining.cml.bg) - 50
)), by = 50), each = 50), 7] <- pointcount$point
staining.cml.bg$OD <- NULL
staining.cml <-
  merge(
    staining.cml,
    staining.cml.bg,
    by = c("genotype", "replicate", "cell.type", "thickness", "point")
  )
staining.cml$OD.adj <-
  round((staining.cml$OD - staining.cml$OD.bg), digits = 4)

staining <- rbind(staining, subset(staining.cml, replicate != 3))
fading.cml <- subset(fading, cell.type == 'CML', select = c(1:12))
fading.cml[, 9] <- revalue(fading.cml[, 9], c('15' = '240'))
staining <- rbind(staining, fading.cml)
staining <-
  staining[order(
    staining$genotype,
    staining$cell.type,
    staining$thickness,
    staining$replicate,
    staining$slice,
    staining$point
  ),]
staining <-
  subset(
    staining,
    replicate != '1_old' &
      replicate != '2_old' &
      replicate != '3_old' & replicate != 'A' &
      genotype == 'col-0' & cell.type != "PH" & cell.type != 'CML'
  )


staining.pre <-
  ddply(
    staining,
    c("genotype", "time", "cell.type", "thickness", "replicate"),
    summarise,
    mean.hue1 = mean(hue, na.rm = TRUE),
    mean.OD.adj1 = mean(OD.adj, na.rm = TRUE),
    mean.OD.bg1 = mean(OD.bg, na.rm = TRUE),
    mean.OD1 = mean(OD, na.rm = TRUE),
    SD.hue1 = sd(hue, na.rm = TRUE),
    SD.OD.adj1 = sd(OD.adj, na.rm = TRUE)
  )

staining.pre.50 <-
  subset(staining.pre, thickness == "50", select = c(1, 2, 3, 5, 7))
colnames(staining.pre.50)[5] <- "ref"
staining.pre <-
  merge(
    staining.pre,
    staining.pre.50,
    by = c("genotype", "time", "cell.type", "replicate"),
    all = TRUE
  )
staining.pre$relative.OD <-
  staining.pre$mean.OD.adj1 / staining.pre$ref
staining.pre$relative.OD[is.na(staining.pre$relative.OD)] <- 1

# staining.avg <- ddply(staining, c("genotype", "time", "cell.type", "thickness"), summarise,
#     mean.hue = mean(hue, na.rm = TRUE),
#     mean.OD.adj = mean(OD.adj, na.rm = TRUE),
#     mean.OD.bg = mean(OD.bg, na.rm = TRUE),
#     mean.OD = mean(OD, na.rm = TRUE),
#     SD.hue = sd(hue, na.rm = TRUE),
#     SD.OD.adj = sd(OD.adj, na.rm = TRUE))

staining.avg <-
  ddply(
    staining.pre,
    c("genotype", "time", "cell.type", "thickness"),
    summarise,
    mean.hue = mean(mean.hue1, na.rm = TRUE),
    mean.OD.adj = mean(mean.OD.adj1, na.rm = TRUE),
    mean.OD.rel = mean(relative.OD, na.rm = TRUE),
    mean.OD = mean(mean.OD1, na.rm = TRUE),
    SD.OD.rel = sd(relative.OD, na.rm = TRUE),
    SD.hue = sd(mean.hue1, na.rm = TRUE),
    SD.OD.adj = sd(mean.OD.adj1, na.rm = TRUE)
  )

stain.hue <-
  ggplot(
    aes(
      y = mean.hue,
      x = time,
      ymin = mean.hue - SD.hue,
      ymax = mean.hue + SD.hue,
      group = thickness
    ),
    data = staining.avg
  ) +
  # geom_errorbar(width = 0.3, size = 0.2) +
  geom_ribbon(aes(fill = thickness), alpha = 0.75) +
  geom_point(
    size = 1.8,
    stroke = 0.2,
    aes(fill = thickness),
    shape = 21
  ) +
  labs(x = "Time [s]", y = "Hue") +
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Helvetica"),
    # axis.line.y = element_line(size = 0.75, lineend = "square"),
    axis.ticks.y = element_line(
      size = 0.25,
      lineend = "square",
      color = "black"
    ),
    axis.ticks.x = element_line(
      size = 0.25,
      lineend = "square",
      color = "black"
    ),
    axis.title.y = element_text(size = 14, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(
      size = 12,
      angle = 45,
      vjust = 1,
      hjust = 0.9,
      colour = "black"
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.25),
    panel.spacing.x = unit(1.5, "mm"),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.position = c(0.12, 0.8),
    # legend.background = element_rect(fill = "grey95", color = NA),
    legend.margin = margin(0, 10, 0, 0),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    strip.text = element_blank()
  ) +
  scale_x_discrete(breaks = c("0", "60", "120", "180", "240")) +
  scale_y_continuous(breaks = c(300, 350, 400, 450),
                     labels = c("300", "350", "40", "90")) +
  #     scale_fill_grey(start = 0, end = 1) +
  scale_fill_brewer(palette = "YlGnBu", name = "Thickness [mm]") +
  facet_wrap(~ cell.type, ncol = 7) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
pdf("time_staining_pre.pdf", width = 12, height = 2)
stain.hue
dev.off()

stain.od <-
  ggplot(
    aes(
      y = mean.OD.adj,
      x = time,
      ymin = mean.OD.adj - SD.OD.adj,
      ymax = mean.OD.adj + SD.OD.adj,
      group = thickness
    ),
    data = staining.avg
  ) +
  # geom_errorbar(width = 0.3, size = 0.2) +
  geom_ribbon(aes(fill = thickness), alpha = 0.75) +
  geom_point(
    size = 1.8,
    stroke = 0.2,
    aes(fill = thickness),
    shape = 21
  ) +
  labs(x = "", y = "Absorbance") +
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Helvetica"),
    legend.position = "none",
    # axis.line.y = element_line(size = 0.75, lineend = "square"),
    axis.ticks.y = element_line(
      size = 0.25,
      lineend = "square",
      color = "black"
    ),
    axis.title.y = element_text(size = 14, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(
      size = 12,
      angle = 45,
      vjust = 1,
      hjust = 0.9,
      colour = "black"
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.25),
    panel.spacing.x = unit(1.5, "mm"),
    strip.text = element_text(
      vjust = 0.1,
      hjust = 0,
      face = "italic"
    ),
    plot.margin = unit(c(0, 0, -0.4, 0), "cm")
  ) +
  scale_x_discrete(breaks = c()) +
  scale_y_continuous(limits = c(-0.05, 0.6),
                     labels = c(" 0.0", " 0.2", " 0.4", " 0.6")) +
  scale_fill_brewer(palette = "YlGnBu", name = "Thickness [mm]") +
  facet_wrap(~ cell.type, ncol = 7)
pdf("time_staining_OD_pre.pdf",
    width = 12,
    height = 2)
stain.od
dev.off()

# plot grid
pdf("staining_grid.pdf", height = 4, width = 10)
stain.plot <-
  plot_grid(
    stain.od,
    stain.hue,
    labels = c('K', 'L'),
    ncol = 1,
    nrow = 2,
    label_fontfamily = "Helvetica",
    rel_heights = c(1, 1.15),
    # label_size = 20,
    hjust = 0,
    vjust = 1
  )
stain.plot
dev.off()

# p <- ggplot(aes(y = mean.OD.rel, x = time, ymin= mean.OD.rel, ymax= mean.OD.rel+SD.OD.rel, group = thickness),
#         data = subset(staining.avg, time == 180 & cell.type != "PH")) +
#     geom_bar(aes(fill = thickness), stat = "identity", position=position_dodge(width = 0.85), width = 0.75, color = "black") +
#     geom_errorbar(colour = "black",position=position_dodge(width = 0.85), width = 0.2) +
#     labs(x = "", y = "Relative Absorbance") +
#     theme_minimal() +
#     theme(text = element_text(size = 10, family = "Helvetica"), panel.grid.minor = element_blank(),
#         legend.position = "none", legend.title = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
#     ylim(0, 2) +
#     scale_x_discrete(breaks =c()) +
# #     scale_fill_grey(start = 0, end = 1) +
#     scale_fill_brewer(palette = "YlGnBu", name = "Thickness [mm]") +
#     facet_wrap(~cell.type, ncol = 7)
# pdf("time_staining_OD_rel.pdf", width = 12, height = 2)
# p
# dev.off()

# hue<-subset(staining,  time == 240 & cell.type == 'IF' & genotype == 'col-0')$hue
# OD<-subset(staining,  time == 240 & cell.type == 'IF' & genotype == 'col-0')$OD.adj
# OD.bg<-subset(staining,  time == 240 & cell.type == 'IF' & genotype == 'col-0')$OD.bg
# thickness<-subset(staining,  time == 240 & cell.type == 'IF' & genotype == 'col-0')$thickness
# corr <- round(cor(OD,hue), digits = 2)
#
# p <- ggplot(aes(y = OD.adj, x = hue, fill = thickness), data = subset(staining,  time == 240 & cell.type == 'IF' & genotype == 'col-0')) +
#     stat_ellipse(geom = "polygon",aes(fill = thickness, color = NULL), alpha = 0.75) +
#     geom_point(size = 1.8, stroke = 0.05, shape = 21, alpha = 0.75)+
#     geom_ribbon(stat = 'smooth', method = "lm", se = TRUE, alpha = 0.1, aes(color = NULL), fill = "grey25") +
#     labs(x = "Hue", y = "Absorbance") +
#     theme_minimal() +
#     theme(text = element_text(size = 10, family = "Helvetica"), panel.grid.minor = element_blank(),
#         legend.position = c(0.15,0.87), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
#     ylim(-0.1, 0.8) +
#     xlim(230, 370) +
#     scale_fill_brewer(palette = "YlGnBu", name = "Thickness [mm]") +
#     annotate("text", x = 300, y = 0.8,
#         label = lm_eqn(lm(OD.adj ~ hue, subset(staining,  time == 240 & cell.type == 'IF' & genotype == 'col-0'))),
#         colour = "black", size = 5, parse = TRUE, family = 'Helvetica', hjust = 0) +
#     annotate("text", x = 300, y = 0.75,
#         label = paste(expression(italic("r")), " == ", corr),
#         colour = "black", size = 5, parse = TRUE, family = 'Helvetica', hjust = 0)
# pdf("hue_OD_thickness_color_IF.pdf")
# p
# dev.off()
#
# hue<-subset(staining,  time == 240 & cell.type == 'MX' & genotype == 'col-0')$hue
# OD<-subset(staining,  time == 240 & cell.type == 'MX' & genotype == 'col-0')$OD.adj
# OD.bg<-subset(staining,  time == 240 & cell.type == 'MX' & genotype == 'col-0')$OD.bg
# thickness<-subset(staining,  time == 240 & cell.type == 'MX' & genotype == 'col-0')$thickness
# corr <- round(cor(OD,hue), digits = 2)
#
# p <- ggplot(aes(y = OD.adj, x = hue, fill = thickness), data = subset(staining,  time == 240 & cell.type == 'MX' & genotype == 'col-0')) +
#     stat_ellipse(geom = "polygon",aes(fill = thickness, color = NULL), alpha = 0.75) +
#     geom_point(size = 1.8, stroke = 0.05, shape = 21, alpha = 0.75)+
#     geom_ribbon(stat = 'smooth', method = "lm", se = TRUE, alpha = 0.1, aes(color = NULL), fill = "grey25") +
#     labs(x = "Hue", y = "Absorbance") +
#     theme_minimal() +
#     theme(text = element_text(size = 10, family = "Helvetica"), panel.grid.minor = element_blank(),
#         legend.position = c(0.15,0.87), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
#     ylim(-0.1, 0.8) +
#     xlim(230, 370) +
#     scale_fill_brewer(palette = "YlGnBu", name = "Thickness [mm]") +
#     annotate("text", x = 300, y = 0.8,
#         label = lm_eqn(lm(OD.adj ~ hue, subset(staining,  time == 240 & cell.type == 'MX' & genotype == 'col-0'))),
#         colour = "black", size = 5, parse = TRUE, family = 'Helvetica', hjust = 0) +
#     annotate("text", x = 300, y = 0.75,
#         label = paste(expression(italic("r")), " == ", corr),
#         colour = "black", size = 5, parse = TRUE, family = 'Helvetica', hjust = 0)
# pdf("hue_OD_thickness_color_MX.pdf")
# p
# dev.off()
#
#
# p <- ggplot(aes(y = OD.adj, x = hue, fill = thickness), data = subset(staining,  time == 180 & genotype == 'col-0' &
#         cell.type != 'PH' & cell.type != 'CML')) +
#     stat_ellipse(geom = "polygon",aes(fill = thickness), alpha = 0.75, color = "black", size = 0.1) +
#     geom_point(size = 1, stroke = 0.05, shape = 21, alpha = 0.5)+
#     labs(x = "Hue", y = "Absorbance") +
#     theme_minimal() +
#     theme(text = element_text(size = 12, family = "Helvetica"), panel.grid.minor = element_blank(),
#         legend.position = c(0.06,0.785), legend.key.size = unit(0.75, 'lines'),
#         strip.text = element_text(vjust = 0.1, hjust = 0, face = "italic"),
#         legend.title = element_text(size = 10)) +
#     ylim(-0.15, 0.8) +
#     xlim(230, 370) +
#     scale_fill_brewer(palette = "YlGnBu", name = "Thickness [mm]") +
#     facet_wrap(~cell.type, ncol = 5)
# pdf("hue_OD_thickness_color_facet.pdf", height = 4, width = 10)
# p
# dev.off()
#
staining$thickness <- paste(staining$thickness, "µm")
staining$thickness <-
  ordered(staining$thickness,
          levels = c("12 µm", "25 µm", "50 µm", "100 µm", "150 µm"))
p <-
  ggplot(
    aes(y = OD.adj, x = hue, fill = cell.type),
    data = subset(
      staining,
      time == 240 & genotype == 'col-0' &
        (cell.type == 'IF' | cell.type == 'MX')
    )
  ) +
  stat_ellipse(
    geom = "polygon",
    aes(fill = cell.type),
    alpha = 0.75,
    color = "black",
    size = 0.1
  ) +
  geom_point(
    size = 1,
    stroke = 0.05,
    shape = 21,
    alpha = 0.5
  ) +
  labs(x = "Hue", y = "Absorbance") +
  theme_minimal() +
  theme(
    text = element_text(size = 15, family = "Helvetica"),
    strip.text = element_text(hjust = 0, face = "italic"),
    # axis.line.y = element_line(size = 0.75, lineend = "square"),
    axis.ticks = element_line(
      size = 0.25,
      lineend = "square",
      color = "black"
    ),
    axis.title = element_text(size = 12),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(
      colour = "black",
      size = 10,
      angle = 0,
      vjust = 1,
      hjust = 0.5
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.25),
    panel.spacing = unit(1.5, "mm"),
    # plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.position = c(0.03, 0.95),
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.height = unit(4, "mm")
  ) +
  ylim(-0.15, 0.8) +
  xlim(230, 370) +
  #     scale_fill_brewer(palette = "RdYlBu", name = "Cell type") +
  scale_fill_grey(start = 0.1, end = 0.9, name = "Cell type") +
  facet_wrap( ~ thickness, ncol = 5)
pdf("hue_OD_thickness_color_facet_inverse.pdf",
    height = 4,
    width = 10)
p
dev.off()
