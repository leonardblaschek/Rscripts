library(tidyverse)
library(sysfonts)
library(showtext)
library(ggthemes)
library(grid)
library(patchwork)
library(ggbeeswarm)

# import Helvetica
font_add(
  "Helvetica",
  regular = "/prop_fonts/01. Helvetica   [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
  italic = "/prop_fonts/01. Helvetica   [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
  bold = "/prop_fonts/01. Helvetica   [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf",
  bolditalic = "/prop_fonts/01. Helvetica   [1957 - Max Miedinger]/HelveticaNeueLTStd-BdIt.otf"
)
showtext_auto()

# generating plot theme
theme_leo <- function(base_size = 8,
                      base_family = "Helvetica") {
  theme_minimal(
    base_size = base_size,
    base_family = base_family
  ) %+replace%
    theme(
      strip.text = element_text(hjust = 0, face = "italic"),
      axis.ticks = element_line(
        size = 0.2,
        lineend = "square",
        color = "black"
      ),
      axis.text.x = element_text(
        size = 8,
        colour = "black", # flipped coords
        margin = margin(1, 1, 1, 1)
      ),
      axis.text.y = element_text(
        colour = "black",
        size = 8,
        angle = 0,
        vjust = 0.5,
        hjust = 1,
        margin = margin(1, 1, 1, 1)
      ),
      axis.title = element_text(
        colour = "black",
        size = 8
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", size = 0.2),
      panel.spacing = unit(1.5, "mm"),
      legend.position = "none",
      legend.text = element_text(size = 8),
      legend.key.height = unit(4, "mm"),
      legend.margin = unit(c(0, 0, 0, 0), "mm"),
      plot.margin = unit(c(2, 2, 2, 2), "mm"),
      complete = TRUE
    )
}

update_geom_defaults("line", list(size = 0.2))
update_geom_defaults("segment", list(size = 0.2))

#### load data ####
bib <-
  read.csv(
    "file:///home/leonard/Documents/Uni/Phloroglucinol/18-06_draft/bibliography/scopus_export_mod_2020.csv"
  ) %>%
  unite("mode", used, relative.quantity, hand.sections, remove = FALSE)

bib_summary <- bib %>%
  summarise(
    total = nrow(bib),
    used = nrow(filter(bib, used == "yes")) / total,
    rel_quant = nrow(filter(bib, used == "yes" & 
                              relative.quantity == "yes")) / total,
    hand_sec = nrow(filter(bib, used == "yes" & 
                             relative.quantity == "yes" &
                             hand.sections == "yes")) / total,
    )
  
  

#### bar graph showing use between the queried journals ####
p <-
  ggplot(data = bib, aes(x = used, fill = relative.quantity)) + geom_bar() + theme_base() +
  theme(
    text = element_text(size = 14, family = "Helvetica"),
    plot.background = element_blank(),
    axis.ticks.y = element_line(
      size = 0.25,
      lineend = "square",
      color = "grey35"
    ),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(
      vjust = 0.1,
      hjust = 0,
      face = "italic",
      angle = 0
    )
  ) +
  scale_fill_brewer(
    palette = "Set1",
    na.value = "grey",
    labels = c("not used", "presence", "quantity")
  ) +
  facet_wrap( ~ Source.title, strip.position = "top")
pdf("bib.pdf")
p
dev.off()

#### boxplot for acidity conditions ####
ac <-
  ggplot(data = subset(bib, used == "yes"), aes(x = Source, y = acidity)) +
  # geom_boxplot(width = 0.5, outlier.alpha = 0) +
  geom_quasirandom(
    aes(fill = mode),
    # width = 0.2,
    shape = 21,
    alpha = 0.75,
    size = 1.5,
    stroke = 0.2
  ) +
  geom_violin(draw_quantiles = 0.5,
              fill = "white",
              colour = "black",
              alpha = 0,
              width = 0.5,
              size = 0.2,
              scale = "width"
  ) +
  scale_fill_manual(values = c("#253494", "#ffffcc")) +
  labs(y = "Surhet [M HCl]") +
  scale_y_continuous(
    limits = c(0, 11),
    breaks = c(0, 2, 4, 6, 8, 10)
  ) +
  scale_fill_manual(values = c("yes_no_no" = "#41b6c4ff", 
                                 "yes_no_" = "#41b6c4ff",
                                 "yes_no_yes" = "#41b6c4ff",
                                 "yes_yes_no" = "#a1dab4ff",
                                 "yes_yes_" = "#a1dab4ff",
                                 "yes_yes_yes" = "#ffffccff"
                                 )) +
  theme_leo() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
  ) +
  coord_flip()

#### boxplot for setion thickness ####
th <-
  ggplot(data = subset(bib, hand.sections == "no"), aes(x = Source, y = thickness)) +
  # geom_boxplot(width = 0.5, outlier.alpha = 0) +
  geom_quasirandom(
    aes(fill = mode),
    # width = 0.2,
    shape = 21,
    alpha = 0.75,
    size = 1.5,
    stroke = 0.2
  ) +
  geom_violin(draw_quantiles = 0.5,
              fill = "white",
              colour = "black",
              alpha = 0,
              width = 0.5,
              size = 0.2,
              scale = "width"
  ) +
  scale_y_continuous(limits = c(0, 350)) +
  scale_fill_manual(values = c("yes_no_no" = "#41b6c4ff", 
                                 "yes_yes_no" = "#a1dab4ff",
                                 "yes_yes_yes" = "#ffffccff")) +
  labs(y = "Snittens tjocklek [µm]") +
  theme_leo() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
  ) +
  coord_flip()

pdf("bib_figs.pdf", width = 8.25 / 2.54, height = 2)
ac / th
dev.off()