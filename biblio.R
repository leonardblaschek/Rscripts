library(sysfonts)
library(showtext)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(dplyr)
library(rowr)
library(plyr)
# library(sunburstR)
library(rPython)
library(ggsunburst)
library(cowplot)
library(png)
library(grid)

# import Helvetica Neue
font_add("Helvetica",
         regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf")
showtext_auto()

bib <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol paper/18-06_draft/bibliography/scopus_export_mod.csv"
  )

bib$snbrst <-
  paste(bib$used, bib$relative.quantity, bib$hand.sections, sep = "-")
bib$snbrst <- revalue(bib$snbrst,
                      c(
                        "no--" = "no",
                        "yes-no-" = "yes-no",
                        "yes-yes-" = "yes-yes"
                      ))

snbrst <- as.data.frame(table(bib$snbrst))
burst <- bib[, 24]


p <-
  ggplot(data = bib, aes(x = used, fill = relative.quantity)) + geom_bar() + theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Helvetica"),
    # axis.line.x = element_line(size = 0.75, lineend = "square"),
    # axis.line.y = element_line(size = 0.75, lineend = "square"),
    axis.ticks.y = element_line(
      size = 0.25,
      lineend = "square",
      color = "grey35"
    ),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    strip.text = element_text(
      vjust = 0.1,
      hjust = 0,
      face = "italic",
      angle = 90
    )
  ) +
  scale_fill_brewer(palette = "Set1", na.value = "grey") +
  facet_wrap(~ Source.title, strip.position = "left")
pdf("bib.pdf")
p
dev.off()

q <-
  sunburstR::sunburst(
    snbrst,
    percent = TRUE,
    count = TRUE,
    colors = c("#2c7bb6", "#d7191c"),
    withD3 = TRUE
  )

ac <-
  ggplot(data = subset(bib, used == "yes"), aes(x = Source, y = acidity)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.2,
              shape = 16,
              alpha = 0.25,
              size = 3) +
  # scale_fill_brewer(palette = "YlGnBu") +
  scale_fill_manual(values = c("#253494", "#ffffcc")) +
  labs(y = "Acidity [M HCl]") + 
  scale_y_continuous(limits = c(0,6), breaks = c(0, 2, 4, 6), labels = c("  0", "  2", "  4", "  6")) +
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Helvetica"),
    legend.position = "none",
    # axis.line.y = element_line(size = 0.75, lineend = "square"),
    axis.ticks.y = element_line(size = 0.25, lineend = "square", color = "black"),
    axis.title.y = element_text(size =10),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.5),
    plot.margin = unit(c(1,1,1,1), "mm"))

th <-
  ggplot(data = subset(bib, hand.sections == "no"), aes(x = Source, y = thickness)) +
  geom_boxplot(width = 0.5, outlier.alpha = 0) +
  geom_jitter(width = 0.2,
              shape = 16,
              alpha = 0.25,
              size = 3) +
  # scale_fill_brewer(palette = "YlGnBu") +
  # scale_fill_manual(values = c("#253494", "#ffffcc")) +
  scale_y_continuous(position = "right") +
  labs(y = "Thickness [µm]") +
  theme_minimal() +
  theme(
    text = element_text(size = 14, family = "Helvetica"),
    legend.position = "none",
    # axis.line.y = element_line(size = 0.75, lineend = "square"),
    axis.ticks.y = element_line(size = 0.25, lineend = "square", color = "black"),
    axis.title.y = element_text(size =10),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.5),
    plot.margin = unit(c(1,1,1,1), "mm"))

bubble <- readPNG("~/R/Output/wiesner/bubble_biblio.png")
p <- rasterGrob(bubble)
pdf("acidity.pdf", width = 6, height = 4)
plot_grid(
  p,
  ac,
  th,
  labels = c('A', 'B', 'C'),
  ncol = 3,
  nrow = 1,
  label_fontfamily = "Helvetica",
  rel_widths = c(4, 1, 1)
)
dev.off()