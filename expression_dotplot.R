library(tidyverse)
library(ggthemes)
library(showtext)
library(ggrepel)


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

expression <- read.csv("file:///home/leonard/Dropbox/Review/figures_PO-review/expression_fig/expression.csv")
expression$PO <- ordered(expression$PO, levels = c("PPO", "LAC", "PRX"))

expression.avg <- expression %>%
  group_by(PO) %>%
  mutate(total = sum(count, na.rm = TRUE)) %>%
  group_by(PO, compartment, method) %>%
  mutate(relative = count/total)

expression.count <- expression %>%
  uncount(weights = count) %>%
  group_by(PO, compartment) %>%
  arrange(method) %>%
  mutate(number = row_number())

ex_plot <- ggplot(expression.count, aes(x = compartment, y = number, fill = PO)) +
  geom_jitter(data = filter(expression.count, method == "pred"),
              aes(x = compartment, y = number, fill = PO),
              alpha = 0) +
  annotate("rect",
           ymin = -Inf,
           ymax = Inf,
           xmin = c(0.5, 2.5, 4.5, 6.5),
           xmax = c(1.5, 3.5, 5.5, 7.5),
           fill = "grey95",
           colour = NA) +
  # geom_vline(xintercept = 1,
  #            size = 7.5,
  #            color = "grey95") +
  # geom_vline(xintercept = 3,
  #            size = 7.5,
  #            color = "grey95") +
  # geom_vline(xintercept = 5,
  #            size = 7.5,
  #            color = "grey95") +
  # geom_vline(xintercept = 7,
  #            size = 7.5,
  #            color = "grey95") +
  geom_jitter(data = filter(expression.count, method == "pred"), aes(x = compartment, y = number, fill = PO),
              width = 0.2,
              shape = 21,
              alpha = 0.5,
              size = 0.5,
              stroke = 0.1) +
  geom_jitter(data = filter(expression.count, method == "exp"), aes(x = compartment, y = number, fill = PO),
              width = 0.2,
              shape = 21,
              alpha = 0.75,
              size = 1,
              stroke = 0.2) +
  # geom_text_repel(data = filter(expression.count, 
  #                               PO == "LAC" & 
  #                                 compartment %in% c("Vacuole", "Thylakoid")),
  #                 label = c("Experimental",
  #                           "Prediction"),
  #                 family = "Helvetica",
  #                 size = ggtext_size,
  #                 # direction = "x",
  #                 hjust = 0,
  #                 vjust = 0.5,
  #                 segment.size = 0.2,
  #                 min.segment.length = 0,
  #                 nudge_y = 20) +
  labs(y = "Reported protein localisations") +
  theme_leo() +
  theme(axis.title.y = element_blank()) +
  scale_fill_few() +
  facet_wrap(~ PO, nrow = 1) +
  coord_flip()

pdf("expression_plot.pdf", height = onecol / 3, width = onecol)
ex_plot
dev.off()
