library(ggplot2)
library(ggthemes)
library(plyr)
library(reshape2)
library(showtext)

font_add("Helvetica", regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf")
showtext_auto()

fading <-
  read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol paper/fading_WT.csv")
fading$time <- (fading$slice * 15) - 15
fading$OD <- ifelse(fading$OD > 1 , NA, fading$OD)
# fading$OD <- 255/(10^fading$OD)
fading$time <- as.factor(as.character(fading$time))
fading$thickness <-
  ordered(as.factor(as.character(fading$thickness)),
          levels = c("12", "25", "50", "100", "150"))
fading.time <- paste(sort(as.integer(levels(fading$time))))
fading$time <-
  ordered(as.factor(as.character(fading$time)), levels = fading.time)
fading[1:50 + rep(seq(0, ((nrow(
  fading
) - 50) - 50), by = 350), each = 50), 2] <- "IF"
fading[51:100 + rep(seq(0, (nrow(fading) - 50), by = 350), each = 50), 2] <-
  "CML"
fading[101:150 + rep(seq(0, (nrow(fading) - 50), by = 350), each = 50), 2] <-
  "MX"
fading[151:200 + rep(seq(0, (nrow(fading) - 50), by = 350), each = 50), 2] <-
  "XF"
fading[201:250 + rep(seq(0, (nrow(fading) - 50), by = 350), each = 50), 2] <-
  "PX"
fading[251:300 + rep(seq(0, (nrow(fading) - 50), by = 350), each = 50), 2] <-
  "LP"
fading[301:350 + rep(seq(0, (nrow(fading) - 50), by = 350), each = 50), 2] <-
  "PH"

# subtract unstained background
fading$cell.type <- ordered(fading$cell.type)
fading$hue <- ((fading$H + 128) / 255) * 360
fading$point <- NA
pointcount <- data.frame(1:50)
colnames(pointcount) <- 'point'
fading[1:50 + rep(seq(0, ((nrow(
  fading
) - 50)), by = 50), each = 50), 10] <- pointcount$point
fading.bg <- subset(fading, time == "0", select = 1:5)
fading.bg$OD.bg <- fading.bg$OD
fading.bg$point <- NA
fading.bg[1:50 + rep(seq(0, ((nrow(
  fading.bg
) - 50)), by = 50), each = 50), 7] <- pointcount$point
fading.bg$OD <- NULL
fading <-
  merge(fading,
        fading.bg,
        by = c("genotype", "replicate", "cell.type", "thickness", "point"))
fading$OD.adj <- round((fading$OD - fading$OD.bg), digits = 4)

# subtract diff in PH to adjust for clearing/bleaching
fading.bg.ph <-
  ddply(
    subset(fading, cell.type == "PH" &
             slice == 2, select = c(1, 2, 3, 4, 12)),
    c("genotype", "replicate", "thickness"),
    summarise,
    OD.bg.ph = mean(OD.adj, na.rm = TRUE)
  )
fading.bg.ph$cell.type <- NULL
fading <-
  merge(
    fading,
    fading.bg.ph,
    all = TRUE,
    by = c("genotype", "replicate", "thickness")
  )
fading$OD.adj <- fading$OD.adj - fading$OD.bg.ph
fading <- subset(fading, replicate != 4)

fading.pre <-
  ddply(
    fading,
    c("genotype", "time", "cell.type", "thickness", "replicate"),
    summarise,
    mean.hue1 = mean(hue, na.rm = TRUE),
    mean.OD.adj1 = mean(OD.adj, na.rm = TRUE),
    mean.OD1 = mean(OD, na.rm = TRUE),
    SD.hue1 = sd(hue, na.rm = TRUE),
    SD.OD.adj1 = sd(OD.adj, na.rm = TRUE)
  )

fading.avg <-
  ddply(
    fading.pre,
    c("genotype", "time", "cell.type", "thickness"),
    summarise,
    mean.hue = mean(mean.hue1, na.rm = TRUE),
    mean.OD.adj = mean(mean.OD.adj1, na.rm = TRUE),
    mean.OD = mean(mean.OD1, na.rm = TRUE),
    SD.hue = sd(mean.hue1, na.rm = TRUE),
    SD.OD.adj = sd(mean.OD.adj1, na.rm = TRUE)
  )


fade.hue <- ggplot(
  aes(y = mean.hue1, x = time, fill = thickness),
  data = subset(
    fading.pre,
    genotype == 'col-0' & cell.type != 'CML' & cell.type != 'PH'
  )
) +
  geom_vline(xintercept = 1,
             size = 15,
             color = "grey95") +
  geom_vline(xintercept = 3,
             size = 15,
             color = "grey95") +
  scale_y_continuous(breaks = c(300, 350, 400, 450),
                     labels = c("300", "350", "40", "90")) +
  scale_x_discrete(
    breaks = c(0, 15, 30, 45),
    labels = c("unstained", "stained", "faded", "restained")
  ) +
  #     geom_errorbar(position="dodge", width=0.2, size=0.2) +
  #     geom_line(size = 0.3, stat = "identity", aes(color=thickness, group=thickness)) +
  #     geom_point(size=2, stroke=0.25, shape=21, position = position_dodge(width = 0.6)) +
  geom_boxplot(size = 0.1, width = 0.5, position = position_dodge(width = 0.75)) +
  labs(x = "", y = "Hue") +
  scale_fill_brewer(palette = "YlGnBu", name = "Thickness [mm]") +
  scale_color_brewer(palette = "YlGnBu") +
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Helvetica"),
    # axis.line.y = element_line(size = 0.75, lineend = "square"),
    axis.ticks.y = element_line(
      size = 0.25,
      lineend = "square",
      color = "black"
    ),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(
      size = 12,
      colour = "black",
      angle = 45,
      vjust = 1,
      hjust = 0.9
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.25),
    panel.spacing.x = unit(1.5, "mm"),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.position = "none",
    strip.text = element_blank()
  ) +
  facet_wrap(~ cell.type, ncol = 7)
pdf("thickness_fade_hue.pdf",
    width = 8,
    height = 2)
fade.hue
dev.off()


fade.od <-
  ggplot(
    aes(
      y = mean.OD.adj,
      x = time,
      ymin = mean.OD.adj,
      ymax = mean.OD.adj + SD.OD.adj,
      fill = thickness
    ),
    data = subset(
      fading.avg,
      genotype == 'col-0' & cell.type != 'CML' & cell.type != 'PH'
    )
  ) +
  geom_vline(xintercept = 1,
             size = 15,
             color = "grey95") +
  geom_vline(xintercept = 3,
             size = 15,
             color = "grey95") +
  geom_errorbar(position = position_dodge(width = 0.75),
                width = 0.5,
                size = 0.1) +
  #     geom_line(size = 0.3, stat = "identity", alpha = 0.4, aes(group=1)) +
  #     geom_point(size=2.5, stroke=0.25, shape=21, fill='#41b6c4') +
  geom_bar(
    width = 0.5,
    stat = "identity",
    position = position_dodge(width = 0.75),
    color = 'black',
    size = 0.1
  ) +
  labs(x = "", y = "Absorbance") +
  scale_fill_brewer(palette = "YlGnBu", name = "Thickness [mm]") +
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Helvetica"),
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
    plot.margin = unit(c(0, 0, -0.4, 0), "cm"),
    legend.position = c(0.06, 0.9),
    # legend.background = element_rect(fill = "grey95", color = NA),
    legend.margin = margin(0, 10, 0, 0),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.size = unit(3, "mm")
  ) +
  #     ylim(0, 0.4) +
  scale_x_discrete(breaks = c()) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6), labels = c(" 0.0", " 0.2", " 0.4", " 0.6"), limits = c(0,0.6)) +
  facet_wrap(~ cell.type, ncol = 7) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
pdf("thickness_fade_OD.pdf", width = 8, height = 1.75)
fade.od
dev.off()

# plot grid
# pdf("fading_thickness_grid.pdf", height = 4, width = 10)
# fade.plot <- plot_grid(
#   fade.od,
#   fade.hue,
#   labels = c('R', 'S'),
#   ncol = 1,
#   nrow = 2,
#   label_fontfamily = "Helvetica",
#   rel_heights = c(1, 1.25),
#   # label_size = 20,
#   hjust = 0,
#   vjust = 1
# )
# fade.plot
# dev.off()
# 
