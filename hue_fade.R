library(ggplot2)
library(ggthemes)
library(plyr)
library(showtext)
library(cowplot)

font_add("Helvetica", regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf")
showtext_auto()

fade <-
  read.csv("/home/leonard/Documents/Uni/Phloroglucinol paper/fading_reworked.csv")
fade$time <-
  ifelse(grepl("1", fade$time), "0", ifelse(grepl("2", fade$time), "24",
                                            ifelse(
                                              grepl("3", fade$time), "48", ifelse(grepl("4", fade$time), "96", ifelse(grepl("5",
                                                                                                                            fade$time), "168", fade$time))
                                            )))
fade$cell.type2[((fade$cell.type > 0 &
                    fade$cell.type < 51) | (fade$cell.type > 250 &
                                              fade$cell.type < 301) |
                   (fade$cell.type > 500 & fade$cell.type < 551) | (fade$cell.type >
                                                                      750 &
                                                                      fade$cell.type < 801) |
                   (fade$cell.type > 1000 & fade$cell.type < 1051)
)] <- "IF"
fade$cell.type2[((fade$cell.type > 50 &
                    fade$cell.type < 101) | (fade$cell.type >
                                               300 &
                                               fade$cell.type < 351) |
                   (fade$cell.type > 550 & fade$cell.type < 601) |
                   (fade$cell.type > 800 &
                      fade$cell.type < 851) | (fade$cell.type > 1050 & fade$cell.type <
                                                 1101)
)] <- "MX"
fade$cell.type2[((fade$cell.type > 100 &
                    fade$cell.type < 151) | (fade$cell.type >
                                               350 &
                                               fade$cell.type < 401) |
                   (fade$cell.type > 600 & fade$cell.type < 651) |
                   (fade$cell.type > 850 &
                      fade$cell.type < 901) | (fade$cell.type > 1100 & fade$cell.type <
                                                 1151)
)] <- "XF"
fade$cell.type2[((fade$cell.type > 150 &
                    fade$cell.type < 201) | (fade$cell.type >
                                               400 &
                                               fade$cell.type < 451) |
                   (fade$cell.type > 650 & fade$cell.type < 701) |
                   (fade$cell.type > 900 &
                      fade$cell.type < 951) | (fade$cell.type > 1150 & fade$cell.type <
                                                 1201)
)] <- "LP"
fade$cell.type2[((fade$cell.type > 200 &
                    fade$cell.type < 251) | (fade$cell.type >
                                               450 &
                                               fade$cell.type < 501) |
                   (fade$cell.type > 700 & fade$cell.type < 751) |
                   (fade$cell.type > 950 &
                      fade$cell.type < 1001) | (fade$cell.type > 1200 & fade$cell.type <
                                                  1251)
)] <- "CB"
fade$cell.type <- fade$cell.type2
fade$cell.type2 <- NULL
fade$hue <- ((fade$H + 128) / 255) * 360
fade$time <-
  ordered(as.factor(as.character(fade$time)), levels = c("0", "24", "48",
                                                         "96", "168"))
fade$OD.bg <- NULL
fade.bg <-
  ddply(
    subset(fade, cell.type == "CB"),
    c("genotype", "time", "cell.type",
      "replicate"),
    summarise,
    OD.bg = mean(OD)
  )
fade.bg$cell.type <- NULL
fade <-
  merge(fade,
        fade.bg,
        all = TRUE,
        by = c("genotype", "replicate", "time"))
fade$OD.adj <- round((fade$OD - fade$OD.bg), digits = 4)
fade$OD.adj[fade$OD.adj > 1] <- NA
fade.avg <-
  ddply(
    fade,
    c("genotype", "time", "cell.type"),
    summarise,
    mean.hue = mean(hue, na.rm = TRUE),
    mean.OD = mean(OD.adj, na.rm = TRUE),
    SD.hue = sd(hue, na.rm = TRUE),
    SD.OD = sd(OD.adj, na.rm = TRUE)
  )

b <-
  ggplot(
    subset(fade.avg, genotype == 'col-0'),
    aes(
      y = mean.hue,
      x = time,
      ymin = mean.hue - SD.hue,
      ymax = mean.hue + SD.hue
    )
  ) +
  scale_y_continuous(breaks = c(330, 360, 390, 420),
                     labels = c("330", "360", "30", "60")) +
  #     geom_errorbar(position="dodge", width=0.2, size=0.2) +
  geom_ribbon(fill = 'grey50',
              alpha = 0.5,
              group = 1) +
  geom_line(size = 0.3,
            stat = "identity",
            alpha = 0.4,
            aes(group = 1)) +
  geom_point(
    size = 3,
    stroke = 0.5,
    shape = 21,
    fill = 'white'
  ) +
  #geom_tufteboxplot(outlier.alpha = 0, outlier.size = 0.1, lwd = 0.2, width = 0.5, position = position_dodge(width = 0.6)) +
  labs(x = "Time [h]", y = "Hue") +
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Helvetica"),
    axis.ticks.y = element_line(size = 0.25, lineend = "square", color = "grey35"),
    axis.ticks.x = element_line(size = 0.25, lineend = "square", color = "grey35"),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(
      size = 10,
      angle = 0,
      vjust = 0.5,
      hjust = 0.5
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "grey35", size = 0.25),
    panel.spacing.x = unit(1.5, "mm"),
    legend.position = "none",
    strip.text = element_blank(),
    plot.margin = unit(c(0,0,0,0), "cm")
  ) +
  facet_wrap( ~ cell.type, ncol = 6)
pdf("time_fade_hue.pdf", width = 8, height = 2)
b
dev.off()


a <-
  ggplot(
    aes(
      y = mean.OD,
      x = time,
      ymin = mean.OD - SD.OD,
      ymax = mean.OD + SD.OD
    ),
    data = subset(fade.avg, genotype == 'col-0')
  ) +
  #     geom_errorbar(position="dodge", width=0.2, size=0.2) +
  geom_ribbon(fill = 'grey50',
              alpha = 0.5,
              group = 1) +
  geom_line(size = 0.3,
            stat = "identity",
            alpha = 0.4,
            aes(group = genotype)) +
  geom_point(
    size = 3,
    stroke = 0.5,
    shape = 21,
    fill = 'white'
  ) +
  #     geom_tufteboxplot(outlier.alpha = 0, outlier.size = 0.1, lwd = 0.2, width = 0.5,
  #     position = position_dodge(widåth = 0.6)) +
  labs(x = "", y = "Absorbance") +
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Helvetica"),
    axis.ticks.y = element_line(size = 0.25, lineend = "square", color = "grey35"),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(
      size = 10,
      angle = 45,
      vjust = 1.1,
      hjust = 1
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "grey35", size = 0.25),
    panel.spacing.x = unit(1.5, "mm"),
    legend.position = "none",
    strip.text = element_text(
      vjust = 0.1,
      hjust = 0,
      face = "italic"
    ),
    plot.margin = unit(c(0,0,-0.4,0), "cm")
  ) +
  #     ylim(0, 0.4) +
  scale_x_discrete(breaks = c()) +
  scale_y_continuous(breaks = c(0,0.1, 0.2, 0.3), labels = c(" 0.0", " 0.1", " 0.2", "0.3")) +
  facet_wrap( ~ cell.type, ncol = 6)
pdf("time_fade_OD.pdf", width = 8, height = 2)
a
dev.off()

# plot grid
pdf("fade_grid.pdf", height = 3.5, width = 10)
plot_grid(a,b, labels = c('A','B'), ncol=1, nrow = 2, label_fontfamily = "Helvetica", rel_heights = c(1,1.1))
dev.off()
