#### setup ####

library(tidyverse)
library(showtext)
library(colorspace)
library(ggrepel)

pal_ostwald_cont <- c(
  "#155DA7",
  "#0A75B9",
  # "#0C89C9",
  "#FED32F",
  # "#F8A63A",
  "#EF663A",
  "#ED4137"
)

pal_ostwald_disc <- c(
  "#275d95",
  "#e8c245",
  "#d25952"
)

pal_ostwald_disc_long <- c(
  "#8fab1d",
  "#2d7d73",
  "#1d566f",
  "#275d95",
  "#e8c245",
  "#d25952"
)

pal_flame <- c(
  "#1F265D",
  "#203F8F",
  "#404678",
  "#897B88",
  "#E7BBA2",
  "#FFEDC3",
  "#FEFEFE"
)

pal_flame_disc <- c(
  "#ffedc3",
  "#897b88",
  "#1f265d"
)

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
      # panel.border = element_rect(fill = NA, color = "black", size = 0.25),
      axis.line = element_line(size = 0.2),
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

ggtext_size <- 6 / (14 / 5)
twocol <- 18 / 2.54
onehalfcol <- 14 / 2.54
onecol <- 9 / 2.54

#### import data ####

data <- read_csv("/home/leonard/Documents/Uni/PhD/Laccase_activity/2020-06_activity/2020-06_activity.csv",
  skip = 1
) %>%
  mutate(
    variable = case_when(
      mean == "NaN" ~ "hue",
      TRUE ~ "absorbance"
    ),
    value = case_when(
      mean == "NaN" ~ median,
      TRUE ~ mean
    ),
    time = slice * interval - interval
  ) %>%
  select(-mean, -median, -slice, -interval)

#### set cell types according to measurement order ####
data[1:20 + rep(seq(0, (nrow(data) - 20), by = 140), each = 20), "cell.type"] <-
  "IF"
data[21:40 + rep(seq(0, (nrow(data) - 20), by = 140), each = 20), "cell.type"] <-
  "LP"
data[41:60 + rep(seq(0, (nrow(data) - 20), by = 140), each = 20), "cell.type"] <-
  "MX"
data[61:80 + rep(seq(0, (nrow(data) - 20), by = 140), each = 20), "cell.type"] <-
  "PX"
data[81:100 + rep(seq(0, (nrow(data) - 20), by = 140), each = 20), "cell.type"] <-
  "XF"
data[101:120 + rep(seq(0, (nrow(data) - 20), by = 140), each = 20), "cell.type"] <-
  "PH"
data[121:140 + rep(seq(0, (nrow(data) - 20), by = 140), each = 20), "cell.type"] <-
  "BG"

#### number measurement points ####
data <- data %>%
  group_by(genotype, replicate, substrate, treatment, date, cell.type, time, variable) %>%
  mutate(point = row_number())

#### subtract water column ####
data <- data %>%
  group_by(time, genotype, replicate, substrate, treatment, date) %>%
  mutate(value = case_when(
    variable == "absorbance" ~
    value - mean(value[variable == "absorbance" & cell.type == "BG"]),
    TRUE ~ value
  ))

#### subtract 0h background ####
data <- data %>%
  group_by(genotype, replicate, substrate, treatment, date, cell.type) %>%
  mutate(value = case_when(
    variable == "absorbance" ~
    value - value[variable == "absorbance" & time == 20],
    TRUE ~ value
  ))

avg_data <- data %>%
  group_by(genotype, replicate, substrate, treatment, date, cell.type, variable, time) %>%
  summarise(value = mean(value))
  
ABTS_plot <- ggplot(
  filter(
    data,
    variable == "absorbance" &
      # time < 150 &
      cell.type != "BG" &
      substrate == "ABTS"
  ),
  aes(
    x = time,
    y = value,
    colour = genotype
  )
) +
  geom_line(aes(group = interaction(point, genotype, cell.type)),
    size = 0.2,
    alpha = 0.25
  ) +
  geom_smooth() +
  scale_colour_manual(values = pal_ostwald_disc[c(1, 3)]) +
  facet_wrap(~ cell.type) +
  theme_leo() +
  theme(legend.position = "bottom")

pdf("ABTS_plot.pdf", width = onecol, height = onecol)
ABTS_plot
dev.off()

DAF_plot <- ggplot(
  filter(
    data,
    variable == "absorbance" &
      time < 200 &
      cell.type != "BG" &
      substrate == "DAF"
  ),
  aes(
    x = time,
    y = value,
    colour = genotype,
    linetype = treatment
  )
) +
  geom_line(aes(group = interaction(point, genotype, treatment, cell.type)),
            size = 0.2,
            alpha = 0.1
  ) +
  geom_text_repel(data = filter(
    avg_data,
    variable == "absorbance" &
      time == 120 &
      cell.type == "MX" &
      treatment == "none" &
      substrate == "DAF"
  ),
  aes(label = genotype),
  family = "Helvetica",
  size = ggtext_size,
  xlim = c(140, 160),
  direction = "y",
  min.segment.length = 0,
  segment.size = 0.2) +
  geom_smooth(size = 0.5) +
  scale_colour_manual(values = pal_ostwald_disc_long[c(1,4,5,6)]) +
  coord_cartesian(clip = "off") +
  scale_linetype_manual(values = c(2, 1)) +
  facet_wrap(~ cell.type) +
  theme_leo() +
  theme(plot.margin = unit(c(0.1, 1, 0.1, 0.1), "cm"))

pdf("DAF_plot.pdf", width = onecol, height = onecol)
DAF_plot
dev.off()

DAB_plot <- ggplot(
  filter(
    data,
    variable == "absorbance" &
      time < 200 &
      cell.type != "BG" &
      substrate == "DAB"
  ),
  aes(
    x = time,
    y = value,
    colour = genotype,
    linetype = treatment
  )
) +
  geom_line(aes(group = interaction(point, genotype, treatment, cell.type)),
            size = 0.2,
            alpha = 0.1
  ) +
  geom_text_repel(data = filter(
    avg_data,
    variable == "absorbance" &
      time == 180 &
      cell.type == "MX" &
      treatment == "none" &
      substrate == "DAB"
  ),
  aes(label = genotype),
  family = "Helvetica",
  size = ggtext_size,
  xlim = c(180, 260),
  direction = "y",
  min.segment.length = 0,
  segment.size = 0.2) +
  geom_smooth(size = 0.5) +
  scale_colour_manual(values = pal_ostwald_disc_long[c(1,4,5,6)]) +
  coord_cartesian(clip = "off") +
  scale_linetype_manual(values = c(2, 1)) +
  facet_wrap(~ cell.type) +
  theme_leo() +
  theme(plot.margin = unit(c(0.1, 1, 0.1, 0.1), "cm"))

pdf("DAB_plot.pdf", width = onecol, height = onecol)
DAB_plot
dev.off()

