library(tidyverse)
library(showtext)

#### import Helvetica ####
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

pal_ostwald_disc <- c(
  "#275d95",
  "#e8c245",
  "#d25952"
)
#### AtLAC4 ####
DOPE_AtLAC4 <- read_csv("/home/leonard/Dropbox/Review/2020_rewrite/Modelling/LAC4/No_ligand/2020-04-01_dopehr-loopmodel/DOPE_AtLAC4.csv",
  col_names = FALSE
) %>%
  add_column(
    "structure" =
      c(
        "template",
        "model_refined1",
        "model_dopehr",
        "model_dopehr_refined"
      ),
    .before = 1
  ) %>%
  pivot_longer(-structure, names_to = "position", values_to = "DOPE") %>%
  mutate(position = as.numeric(str_remove(position, "X"))) 

DOPE_scores <- ggplot(
  DOPE_AtLAC4,
  aes(
    x = position,
    y = DOPE,
    colour = structure
  )
) +
  annotate("rect",
    xmin = c(292, 320, 351, 399),
    xmax = c(309, 339, 365, 421),
    ymin = -Inf,
    ymax = Inf,
    colour = NA,
    fill = "grey90"
  ) +
  geom_line(
    size = 0.4,
    alpha = 0.75
  ) +
  scale_colour_brewer(palette = "Set2") +
  theme_leo() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

pdf("DOPE_AtLAC4.pdf", width = 4, height = 2)
DOPE_scores
dev.off()

#### AtLAC11 ####
DOPE_AtLAC11 <- read_csv("/home/leonard/Dropbox/Review/2020_rewrite/Modelling/LAC11/No_ligand/2020-04-06_dopehr_loopmodel/DOPE_AtLAC11.csv",
                        col_names = FALSE
) %>%
  add_column(
    "structure" =
      c(
        "template",
        "model_3_4",
        "model_5_4",
        "model_5_3"
      ),
    .before = 1
  ) %>%
  pivot_longer(-structure, names_to = "position", values_to = "DOPE") %>%
  mutate(position = as.numeric(str_remove(position, "X"))) 

DOPE_scores <- ggplot(
  DOPE_AtLAC11,
  aes(
    x = position,
    y = DOPE,
    colour = structure
  )
) +
  # annotate("rect",
  #          xmin = c(292, 320, 351, 399),
  #          xmax = c(309, 339, 365, 421),
  #          ymin = -Inf,
  #          ymax = Inf,
  #          colour = NA,
  #          fill = "grey90"
  # ) +
  geom_line(
    size = 0.4,
    alpha = 0.75
  ) +
  scale_colour_brewer(palette = "Set2") +
  theme_leo() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

pdf("DOPE_AtLAC11.pdf", width = 4, height = 2)
DOPE_scores
dev.off()
