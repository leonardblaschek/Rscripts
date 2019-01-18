library(png)
library(grid)
library(sysfonts)
library(showtext)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
# library(reshape2)
library(agricolae)
library(plyr)
library(dplyr)
# library(rowr)
library(cowplot)


####### import Helvetica Neue #### ####
font_add(
  "Helvetica",
  regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
  italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
  bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf",
  bolditalic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-BdIt.otf"
)
showtext_auto()

#### import staining measurements from video stills ####
staining <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/17-10_timelines/video_reworked.csv"
  )


#### slices are in 15 s intervals, starting at 0 ####
staining$time <- (staining$slice * 15) - 15
staining$OD <- ifelse(staining$OD > 1 , NA, staining$OD)
staining$time <- as.factor(as.character(staining$time))
staining$thickness <-
  ordered(as.factor(as.character(staining$thickness)),
          levels = c("12", "25", "50", "100", "150"))
staining.time <- paste(sort(as.integer(levels(staining$time))))
staining$time <-
  ordered(as.factor(as.character(staining$time)), levels = staining.time)


#### set cell types according to measurement order ####
staining[1:50 + rep(seq(0, (nrow(staining) - 50), by = 300), each = 50), 2] <- 
  "IF"
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


#### calculate correct hue value ####
staining$hue <- ((staining$H + 128) / 255) * 360


#### subtract background ####
staining$point <- NA
pointcount <- data.frame(1:50)
colnames(pointcount) <- 'point'
staining[1:50 + rep(seq(0, ((nrow(staining) - 50)), by = 50), each = 50), 10] <- pointcount$point

staining.bg <- subset(staining, time == "0", select = 1:5)
staining.bg$OD.bg <- staining.bg$OD
staining.bg$point <- NA
staining.bg[1:50 + rep(seq(0, ((nrow(staining.bg) - 50)), by = 50), each = 50), 7] <- pointcount$point
staining.bg$OD <- NULL

staining <-
  merge(
    staining,
    staining.bg,
    by = c("genotype", "replicate", "cell.type", "thickness", "point")
  )
staining$OD.adj <- round((staining$OD - staining$OD.bg), digits = 4)

staining.bg.ph <-
  ddply(
    subset(staining, cell.type == "PH", select = c(1, 2, 3, 4, 12)),
    c("genotype", "replicate", "thickness"),
    summarise,
    OD.bg.ph = mean(OD.adj, na.rm = TRUE)
  )
staining.bg.ph$cell.type <- NULL

staining <-
  merge(
    staining,
    staining.bg.ph,
    all = TRUE,
    by = c("genotype", "replicate", "thickness")
  )
staining$OD.adj <- staining$OD.adj - staining$OD.bg.ph

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
      genotype == 'col-0' & cell.type != "PH"
  )

staining$thickness <- paste(staining$thickness, "µm")
staining$thickness <-
  ordered(staining$thickness,
          levels = c("12 µm", "25 µm", "50 µm", "100 µm", "150 µm"))

#### import fading measurements ####
fading <-
  read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol/fading_WT.csv")


#### slices are in 15 s intervals, starting at 0 ####
fading$time <- (fading$slice * 15) - 15
fading$OD <- ifelse(fading$OD > 1 , NA, fading$OD)
fading$time <- as.factor(as.character(fading$time))
fading$thickness <-
  ordered(as.factor(as.character(fading$thickness)),
          levels = c("12", "25", "50", "100", "150"))
fading.time <- paste(sort(as.integer(levels(fading$time))))
fading$time <-
  ordered(as.factor(as.character(fading$time)), levels = fading.time)


#### subtract background ####
fading[1:50 + rep(seq(0, (nrow(fading) - 50), by = 350), each = 50), 2] <- 
  "IF"
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

fading$hue <- ((fading$H + 128) / 255) * 360

fading$cell.type <- ordered(fading$cell.type)

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

timeline_stain <- subset(staining, select = c(2, 3, 4, 9, 10, 12))
timeline_fade <- subset(fading, cell.type != "CML" & time != 0 & time != 15, select = c(2, 3, 4, 9, 10, 12))
timeline_fade$thickness <- revalue(timeline_fade$thickness, c(
  "12" = "12 µm",
  "25" = "25 µm",
  "50" = "50 µm",
  "100" = "100 µm",
  "150" = "150 µm")
)

#### manipulate values for correct placement on the continuous scale ####
timeline_fade$time <- revalue(timeline_fade$time, c(
  "30" = "241",
  "45" = "242")
)
timeline_fade$replicate <- as.factor(as.character(timeline_fade$replicate))

timeline <- merge(timeline_stain, timeline_fade, all = TRUE)
timeline$time <- as.numeric(as.character(timeline$time))
timeline <- subset(timeline, time < 245 & cell.type != "PH")
timeline$time[timeline$time == 241] <- 270
timeline$time[timeline$time == 242] <- 300
names(timeline)[6] <- "od"

#### average by replicate and then by genotype ####
timeline.pre <-
  ddply(timeline,
        c("replicate", "thickness", "time", "cell.type"),
        summarise,
        mean.hue.pre = mean(hue, na.rm = TRUE),
        mean.od.pre = mean(od, na.rm = TRUE))

timeline.avg <-
  ddply(timeline.pre,
        c("thickness", "time", "cell.type"),
        summarise,
        mean.hue = mean(mean.hue.pre),
        sd.hue = sd(mean.hue.pre),
        mean.od = mean(mean.od.pre),
        sd.od = sd(mean.od.pre))

#### plot the absorbance timeline ####
tl.od <-
  ggplot(timeline.avg,
         aes(x = time, y = mean.od, ymin = mean.od - sd.od, ymax = mean.od + sd.od, group = thickness)
  ) +
  geom_errorbar(data = subset(timeline.avg, time == 270 | time == 300), width = 10, size = 0.2) +
  geom_ribbon(data = subset(timeline.avg, time != 270 & time != 300), aes(fill = thickness), alpha = 0.5) +
  geom_vline(xintercept = 255, linetype = 2) +
  geom_point(aes(fill = thickness), shape = 21, size = 2, stroke = 0.2) +
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Helvetica"),
    axis.ticks.y = element_line(
      size = 0.25,
      lineend = "square",
      color = "black"
    ),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.25),
    panel.spacing.x = unit(1.5, "mm"),
    plot.margin = unit(c(0.4, 0, 0, 0.1), "cm"),
    legend.position = c(0.194, 1.2),
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(2, "mm"),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    strip.text = element_text(
      vjust = 0.1,
      hjust = 0,
      face = "italic"
    )
  ) +
  labs(y = "Absorbance") +
  scale_fill_brewer(palette = "YlGnBu", name = "Thickness [mm]") +
  scale_x_continuous(breaks = c(0, 60, 120, 180, 240, 270, 300), labels = c("0 s", "60 s", "120 s", "180 s", "240 s", "24 h", "re-stained")) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4), labels = c(" 0.0", " 0.2", " 0.4")) +
  facet_wrap(~ cell.type, ncol = 5) 

pdf("timeline_od.pdf", height = 4, width = 10)
tl.od
dev.off()

#### plot the hue timeline ####
tl.hue <-
  ggplot(timeline.avg,
         aes(x = time, y = mean.hue, ymin = mean.hue - sd.hue, ymax = mean.hue + sd.hue, group = thickness)
  ) +
  geom_errorbar(data = subset(timeline.avg, time == 270 | time == 300), width = 10, size = 0.2) +
  geom_ribbon(data = subset(timeline.avg, time != 270 & time != 300), aes(fill = thickness), alpha = 0.5) +
  geom_vline(xintercept = 255, linetype = 2) +
  geom_point(aes(fill = thickness), shape = 21, size = 2, stroke = 0.2) +
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Helvetica"),
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
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(
      size = 12,
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      colour = "black"
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.25),
    panel.spacing.x = unit(1.5, "mm"),
    plot.margin = unit(c(0, 0, 0, 0.1), "cm"),
    legend.position = "none",
    strip.text = element_blank()
  ) +
  labs(y = "Hue") +
  scale_fill_brewer(palette = "YlGnBu", name = "Thickness [mm]") +
  scale_x_continuous(breaks = c(0, 60, 120, 180, 240, 270, 300), labels = c("0 s", "60 s", "120 s", "180 s", "240 s", "24 h", "re-stained")) +
  scale_y_continuous(breaks = c(300, 350, 400, 450), labels = c("300", "350", "40", "90")) +
  facet_wrap(~ cell.type, ncol = 5) 
pdf("timeline_hue.pdf", height = 4, width = 10)
tl.hue
dev.off()

#### plot the timeline grid ####
pdf("timeline_grid.pdf", width = 10, height = 4)
tl.grid <- plot_grid(tl.od, 
                     tl.hue,
                     labels = c('(c)', '(d)'),
                     ncol = 1,
                     nrow = 2,
                     label_fontfamily = "Helvetica",
                     rel_heights = c(1, 1.225),
                     label_size = 18,
                     hjust = 0,
                     vjust = 1
)
tl.grid
dev.off()

#### plot the ellipses for absorbance and hue in MX and IF against the section thickness ####
thick.ellipse <-
  ggplot(
    aes(x = OD.adj, y = hue, fill = cell.type),
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
  labs(y = "Hue", x = "Absorbance") +
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Helvetica"),
    strip.text = element_text(hjust = 0, face = "italic"),
    axis.ticks = element_line(
      size = 0.25,
      lineend = "square",
      color = "black"
    ),
    axis.title = element_text(size = 12),
    axis.text.y = element_text(size = 12, colour = "black"),
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
    legend.position = c(0.17, 0.15),
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.spacing = unit(0.25, "mm"),
    legend.key.height = unit(4, "mm"),
    legend.key.width = unit(4, "mm"),
    plot.margin = unit(c(0, 0, 0, 0.1), "cm")
  ) +
  xlim(-0.15, 0.8) +
  ylim(230, 370) +
  scale_fill_grey(start = 0.1, end = 0.9, name = "Cell type", labels = c(" IF", " MX")) +
  facet_wrap( ~ thickness, ncol = 5)

pdf("hue_OD_thickness_color_facet_inverse.pdf",
    height = 2,
    width = 10)
thick.ellipse
dev.off()


#### import images ####
timeline_imgs <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/timeline/cropped/montage.png"))


#### plot the whole figure ####
pdf("fig4.pdf", height = 12, width = 10)
plot_grid(
  timeline_imgs,
  thick.ellipse,
  tl.grid,
  labels = c("(a)", "(b)", ""),
  label_fontfamily = "Helvetica",
  nrow = 3,
  ncol = 1,
  label_size = 18,
  rel_heights = c(1.2, 0.5, 1),
  hjust = 0,
  vjust = 1
)
dev.off()

#### reduce figure size via ghostscript ####
system("gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/printer -dNOPAUSE -dQUIET -dBATCH -dDetectDuplicateImages -dCompressFonts=true -r150 -sOutputFile=stain_fade_grid.pdf fade_stain.pdf")