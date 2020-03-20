library(tidyverse)
library(ggthemes)
library(showtext)
library(cowplot)
library(ggbeeswarm)
library(patchwork)


#### import Helvetica Neue ####
font_add(
  "Helvetica",
  regular = "/prop_fonts/01. Helvetica   [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
  italic = "/prop_fonts/01. Helvetica   [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
  bold = "/prop_fonts/01. Helvetica   [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf",
  bolditalic = "/prop_fonts/01. Helvetica   [1957 - Max Miedinger]/HelveticaNeueLTStd-BdIt.otf"
)
showtext_auto()

#### generating plot theme ####
theme_leo <- function(base_size = 6,
                      base_family = "Helvetica") {
  theme_minimal(
    base_size = base_size,
    base_family = base_family
  ) %+replace%
    theme(
      strip.text = element_text(hjust = 0, face = "italic"),
      axis.ticks = element_line(
        size = 0.125,
        lineend = "square",
        color = "black"
      ),
      axis.text.x = element_text(
        size = 6,
        colour = "black", # flipped coords
        margin = margin(1, 1, 1, 1)
      ),
      axis.text.y = element_text(
        colour = "black",
        size = 6,
        angle = 0,
        vjust = 0.5,
        hjust = 1,
        margin = margin(1, 1, 1, 1)
      ),
      axis.title = element_text(
        colour = "black",
        size = 6
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", size = 0.25),
      panel.spacing = unit(1.5, "mm"),
      legend.position = "none",
      legend.text = element_text(size = 6),
      legend.key.height = unit(4, "mm"),
      plot.title = element_text(
        size = 6,
        hjust = 0
      ),
      complete = TRUE
    )
}

equal_breaks <- function(n = 3, s = 0.05, ...) {
  function(value) {
    # rescaling
    d <- s * diff(range(value)) / (1 + 2 * s)
    round(seq(0, max(value) - d, length = n), digits = 2)
  }
}

my_breaks <- function(value) {
  if (
    max(value) < 0.005) {
    seq(0, 0.004, length.out = 3)
  } else if (
    max(value) < 0.05) {
    seq(0, 0.02, length.out = 3)
  } else if (
    max(value) < 2) {
    seq(0, 1, length.out = 3)
  } else {
    seq(0, 100, length.out = 3)
  }
}

my_limits <- function(value) {
  if (
    max(value) < 0.005) {
    c(0, 0.0045)
  } else if (
    max(value) < 0.05) {
    c(0, 0.025)
  } else if (
    max(value) < 2) {
    c(-0.25, 1.1)
  } else {
    c(-10, 120)
  }
}

ggtext_size <- 6 / (14 / 5)
twocol <- 18 / 2.54
onehalfcol <- 14 / 2.54
onecol <- 9 / 2.54

#### import data ####

lac_data <- read_csv("/home/leonard/Dropbox/Review/2020_rewrite/Data/lac_data.csv") %>%
  select(1, 2, 8, 10:24) %>%
  pivot_longer(-c("Kingdom", "Species", "Temperature.optimum", "PI"),
    names_to = "variable", values_to = "value"
  ) %>%
  separate(variable, c("var", "substrate"), extra = "merge", remove = TRUE) %>%
  mutate(
    value = as.numeric(value),
    Temperature.optimum = as.numeric(Temperature.optimum),
    Kingdom = ordered(Kingdom, levels = c("Fungi", "Prokaryota", "Plantae"))
  ) %>%
  filter(substrate %in% c("ABTS", "DMP", "SGZ"))

#### Km figure ####

lac.fig.km <- ggplot(
  data = filter(lac_data, var == "Km"),
  aes(x = Kingdom, y = value, fill = Kingdom)
) +
  geom_quasirandom(
    aes(colour = Kingdom),
    dodge.width = 0.6,
    shape = 16,
    alpha = 0.75,
    size = 1
  ) +
  geom_violin(
    draw_quantiles = 0.5,
    fill = "white",
    colour = "black",
    alpha = 0.5,
    width = 0.4,
    position = position_dodge(width = 0.6),
    size = 0.2,
    scale = "width"
  ) +
  scale_y_log10(limits = c(0.5, 10000)) +
  labs(y = "Km [µM]") +
  scale_colour_brewer(palette = "Set1") +
  theme_leo() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  facet_wrap(~substrate)

#### pH figure ####

lac.fig.ph <- ggplot(data = filter(lac_data, var == "pH"), aes(x = Kingdom, y = value, fill = Kingdom)) +
  geom_quasirandom(
    aes(colour = Kingdom),
    dodge.width = 0.6,
    shape = 16,
    alpha = 0.75,
    size = 1
  ) +
  geom_violin(
    draw_quantiles = 0.5,
    fill = "white",
    colour = "black",
    alpha = 0.5,
    width = 0.4,
    position = position_dodge(width = 0.6),
    size = 0.2,
    scale = "width"
  ) +
  labs(y = "pH optimum") +
  scale_colour_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 10), breaks = c(0, 3, 6, 9)) +
  theme_leo() +
  theme(strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  facet_wrap(~substrate)

#### Temp figure ####

lac.fig.temp <- ggplot(data = lac_data, aes(x = Kingdom, y = Temperature.optimum, fill = Kingdom)) +
  geom_quasirandom(
    aes(colour = Kingdom),
    dodge.width = 0.6,
    shape = 16,
    alpha = 0.75,
    size = 1
  ) +
  geom_violin(
    draw_quantiles = 0.5,
    fill = "white",
    colour = "black",
    alpha = 0.5,
    width = 0.4,
    position = position_dodge(width = 0.6),
    size = 0.2,
    scale = "width"
  ) +
  annotate("text",
           x = c(1, 2, 3),
           y = c(5, 15, 25),
           label = c("Fungi", "Prokaryota", "Plantae"),
           family = "Helvetica",
           size = ggtext_size
  ) +
  annotate("segment",
           x = c(1, 2, 3),
           xend = c(1, 2, 3),
           y = c(22, 27, 37),
           yend = c(8, 18, 28),
           colour = "black",
           size = 0.2
  ) +
  labs(y = "Temperature optimum [°C]") +
  scale_colour_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 95)) +
  theme_leo() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

#### pI figure ####

lac.fig.pi <- ggplot(data = lac_data, aes(x = Kingdom, y = PI, fill = Kingdom)) +
  geom_quasirandom(
    aes(colour = Kingdom),
    dodge.width = 0.6,
    shape = 16,
    alpha = 0.75,
    size = 1
  ) +
  geom_violin(
    draw_quantiles = 0.5,
    fill = "white",
    colour = "black",
    alpha = 0.5,
    width = 0.4,
    position = position_dodge(width = 0.6),
    size = 0.2,
    scale = "width"
  ) +
  labs(y = "Isoelectric point") +
  scale_colour_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 10.5), breaks = c(0, 3, 6, 9)) +
  theme_leo() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

#### Arrange panels ####

pdf("lac_fig.pdf", width = onecol, height = onecol * 0.75)
lac.fig.temp + lac.fig.km + lac.fig.pi + lac.fig.ph +
  plot_layout(widths = c(1 , 2.5))
dev.off()
