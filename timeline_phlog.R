library(png)
library(grid)
library(sysfonts)
library(showtext)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(reshape2)
library(agricolae)
library(plyr)
library(dplyr)
library(rowr)
library(cowplot)

# import Helvetica Neue
font_add(
  "Helvetica",
  regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
  italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
  bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf",
  bolditalic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-BdIt.otf"
)
showtext_auto()

timeline <-
  read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol paper/timeline_reworked.csv")

timeline[1:50 + rep(seq(0, (nrow(timeline) - 50), by = 300), each = 50), 3] <-
  "IF"
timeline[51:100 + rep(seq(0, (nrow(timeline) - 50), by = 300), each = 50), 3] <-
  "MX"
timeline[101:150 + rep(seq(0, (nrow(timeline) - 50), by = 300), each = 50), 3] <-
  "XF"
timeline[151:200 + rep(seq(0, (nrow(timeline) - 50), by = 300), each = 50), 3] <-
  "PX"
timeline[201:250 + rep(seq(0, (nrow(timeline) - 50), by = 300), each = 50), 3] <-
  "LP"
timeline[251:300 + rep(seq(0, (nrow(timeline) - 50), by = 300), each = 50), 3] <-
  "PH"

timeline$slice <- as.factor(as.character(timeline$slice))
timeline$thickness <- ordered(as.factor(as.character(timeline$thickness)), levels = c("12", "25", "50", "100", "150"))

timeline.slices <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol paper/timeline_reworked_slices.csv"
  )
timeline.slices$slice <-
  as.factor(as.character(timeline.slices$slice))
timeline.slices$thickness <-
  as.factor(as.character(timeline.slices$thickness))
timeline.slices$time <- as.numeric(as.character(revalue(timeline.slices$time, c("24h" = "241", "24hs" = "242"))))

timeline <- merge(timeline, timeline.slices, all = TRUE)
timeline.hue <- timeline[!timeline$value < 1,]
names(timeline.hue)[6] <- "hue"
timeline <- timeline[!timeline$value > 1,]
names(timeline)[6] <- "od"
timeline <- merge(timeline, timeline.hue,
                  by = c("replicate", "thickness", "slice", "cell.type", "point", "time"))

timeline$cell.type <- as.factor(as.character(timeline$cell.type))
timeline$point <- as.factor(as.character(timeline$point))
timeline$replicate <- as.factor(as.character(timeline$replicate))

timeline.bg <- subset(timeline, time == "0")
timeline.bg$od.bg <- timeline.bg$od

timeline <-
  merge(timeline,
        timeline.bg[, c(1, 2, 4, 5, 9)],
        by = c("replicate", "thickness", "cell.type", "point"))

timeline$od <- timeline$od - timeline$od.bg
timeline$hue <- ((timeline$hue + 128) / 256) * 360

timeline.ph <- subset(timeline, cell.type == "PH")
timeline.ph$od.ph <- timeline.ph$od
timeline.ph$od.ph[timeline.ph$od.ph > 0.2] <- NA
timeline.ph <- ddply(timeline.ph, c("replicate", "thickness", "slice"), summarise, od.ph = mean(od.ph, na.rm = TRUE))
timeline.ph$od.ph[is.nan(timeline.ph$od.ph)] <- 0

timeline <- merge(timeline,
                  timeline.ph,
                  by = c("replicate", "thickness", "slice"), all = TRUE)
timeline$od <- timeline$od - timeline$od.ph

timeline <- subset(timeline, time < 245 & cell.type != "PH" & cell.type != "LP")
timeline$time[timeline$time == 241] <- 270
timeline$time[timeline$time == 242] <- 285

timeline.pre <-
  ddply(timeline,
        c("replicate", "thickness", "time", "cell.type"),
        summarise,
        mean.hue.pre = mean(hue),
        mean.od.pre = mean(od))

timeline.avg <-
  ddply(timeline.pre,
        c("thickness", "time", "cell.type"),
        summarise,
        mean.hue = mean(mean.hue.pre),
        sd.hue = sd(mean.hue.pre),
        mean.od = mean(mean.od.pre),
        sd.od = sd(mean.od.pre))

tl.od <-
  ggplot(timeline.avg,
    aes(x = time, y = mean.od, ymin = mean.od - sd.od, ymax = mean.od + sd.od)
  ) +
  geom_errorbar(data = subset(timeline.avg, time == 270 | time == 285), width = 5) +
  geom_ribbon(data = subset(timeline.avg, time != 270 & time != 285), aes(fill = thickness), alpha = 0.5) +
  # geom_ribbon(data = subset(timeline.avg, time == 270 | time == 285), aes(fill = thickness), alpha = 0.5) +
  geom_vline(xintercept = 255, linetype = 2) +
  # geom_line(data = subset(timeline.avg, time != 270 & time != 285), aes(colour = thickness), size = 1) +
  # geom_point(data = subset(timeline.avg, time == 270 | time == 285), aes(fill = thickness), size = 3, shape = 21) +
  geom_point(aes(fill = thickness), shape = 21, size = 2) +
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
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, colour = "black"),
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
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.position = "top",
    # legend.background = element_rect(fill = "grey95", color = NA),
    legend.margin = margin(0, 10, 0, 0),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    strip.text = element_text(
      vjust = 0.1,
      hjust = 0,
      face = "italic"
    )
  ) +
  scale_fill_brewer(palette = "YlGnBu", name = "Thickness [mm]") +
  # scale_fill_manual(values = c("#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84")) +
  scale_colour_brewer(palette = "YlGnBu", name = "Thickness [mm]") +
  # scale_colour_manual(values = c("#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84")) +
  scale_x_continuous(breaks = c(0, 60, 120, 180, 240, 270, 285), labels = c("0 s", "60 s", "120 s", "180 s", "240 s", "24 h", "re-stained")) +
  facet_wrap(~ cell.type, ncol = 4)

pdf("timeline_od.pdf", height = 4, width = 10)
tl.od
dev.off()

tl.hue <-
  ggplot(timeline.pre,
         aes(x = time, y = mean.hue.pre, fill = thickness)
  ) +
  geom_vline(xintercept = 415, linetype = 2) +
  geom_line(data = subset(timeline.pre, time != 430 & time != 460)) +
  geom_point(shape = 21, size = 3) +
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
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, colour = "black"),
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
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.position = "top",
    # legend.background = element_rect(fill = "grey95", color = NA),
    legend.margin = margin(0, 10, 0, 0),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    strip.text = element_text(
      vjust = 0.1,
      hjust = 0,
      face = "italic"
    )
  ) +
  scale_fill_brewer(palette = "YlGnBu", name = "Thickness [mm]") +
  scale_colour_brewer(palette = "YlGnBu", name = "Thickness [mm]") +
  scale_x_continuous(breaks = c(0, 100, 200, 300, 400, 430, 460), labels = c("0 s", "100 s", "200 s", "300 s", "400 s", "24 h", "re-stained")) +
  facet_wrap(~ cell.type, ncol = 3)

pdf("timeline_hue.pdf", height = 4, width = 10)
tl.hue
dev.off()