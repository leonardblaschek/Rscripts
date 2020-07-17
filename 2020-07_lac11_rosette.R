#### setup ####
library(tidyverse)
library(jsonlite)

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

# test <- fromJSON("/home/leonard/Documents/Uni/PhD/Phenotyping/2020-05_lac11_segregating/34_3.jpg_13.json")
# test$observations$hue_circular_mean
# test$observations$area$value

#### read-in function ####
read_plus <- function(flnm) {
  json <- fromJSON(flnm)
  json$observations$area$value %>%
    as_tibble() %>%
    rename("raw_area" = "value") %>%
    bind_cols(as_tibble(json$observations$hue_circular_mean$value)) %>%
    rename("mean_hue" = "value") %>%
    mutate(filename = basename(flnm))
}

#### load data ####
rosette_files <-
  list.files(
    path = "/home/leonard/Documents/Uni/PhD/Phenotyping/2020-05_lac11_segregating/",
    pattern = "*.json",
    recursive = FALSE,
    full.names = TRUE
  )

size_marker <- read_csv("/home/leonard/Documents/Uni/PhD/Phenotyping/2020-05_lac11_segregating/size_marker_trays.csv",
  col_names = c("filename", "space")
) %>%
  separate(filename, into = c("dpg", "tray"), sep = "_") %>%
  mutate(
    tray = str_remove(tray, fixed(".jpg"))
  )

rosette_area <- map(rosette_files, read_plus) %>%
  bind_rows() %>%
  separate(filename, into = c("dpg", "tray", "plant"), sep = "_") %>%
  mutate(
    tray = str_remove(tray, fixed(".jpg")),
    plant = str_remove(plant, fixed(".json"))
  ) %>%
  left_join(size_marker) %>%
  left_join(read_csv("/home/leonard/Documents/Uni/PhD/Phenotyping/2020-05_lac11_segregating/tray_id.csv",
                     col_types = "cccc")) %>%
  mutate(
    area = raw_area / space,
    dpg = as.numeric(dpg)
  ) %>%
  filter(phenotype != "X")


#### plot ####
rosette_growth <- ggplot(
  data = rosette_area,
  aes(
    x = dpg,
    y = area
  )
) +
  geom_line(aes(
    group = interaction(tray, plant),
    colour = genotype
  ),
  size = 0.2
  ) +
  theme_leo() +
  coord_cartesian(clip = "off") +
  scale_colour_manual(values = pal_ostwald_disc, na.value = rgb(0, 0, 0, alpha = 0.15)) +
  labs(
    x = "DPG",
    y = "Rosette area [scaled px]"
  ) +
  theme(legend.position = "top",
        plot.margin = unit(c(0.1, 0.11, 0.1, 0.1), "cm"))

pdf("rosette_growth.pdf", width = onecol, height = onecol)
rosette_growth
dev.off()

rosette_hue <- ggplot(
  data = rosette_area,
  aes(
    x = dpg,
    y = mean_hue
  )
) +
  geom_line(aes(
    group = interaction(tray, plant),
    colour = genotype
  ),
  size = 0.2
  ) +
  theme_leo() +
  coord_cartesian(clip = "off") +
  scale_colour_manual(values = pal_ostwald_disc, na.value = rgb(0, 0, 0, alpha = 0.15)) +
  labs(
    x = "DPG",
    y = "Rosette mean hue"
  ) +
  theme(legend.position = "top",
        plot.margin = unit(c(0.1, 0.11, 0.1, 0.1), "cm"))

pdf("rosette_hue.pdf", width = onecol, height = onecol)
rosette_hue
dev.off()

rosette_pca_data <- rosette_area %>%
  select(tray, plant, dpg, area, mean_hue, phenotype) %>%
  complete(dpg, nesting(tray, plant, phenotype)) %>%
  group_by(tray, plant, phenotype) %>%
  fill(mean_hue, area) %>%
  drop_na() %>%
  filter(dpg > 6) %>%
  pivot_wider(id_cols = c(tray, plant, phenotype),
              names_from = c(dpg),
              values_from = c(area, mean_hue)) %>%
  filter(!(tray == 4 & plant == 10)) %>%
  unite("ID", tray, plant, phenotype, sep = "_")

pca_result <- prcomp(rosette_pca_data[, -c(1, 2, 3)], center = TRUE, scale. = TRUE)
rownames(pca_result$x) <- rosette_pca_data$ID

gg_pca <- as_tibble(pca_result$x, rownames = NA) %>%
  rownames_to_column(var = "ID") %>%
  separate(ID, sep = "_", into = c("tray", "plant", "phenotype"), remove = FALSE)

pca_plot <- ggplot(gg_pca, aes(x = PC1, y = PC2, fill = phenotype)) +
  geom_hline(yintercept = 0, linetype = 1, size = 0.2) +
  geom_vline(xintercept = 0, linetype = 1, size = 0.2) +
  stat_ellipse(geom = "polygon", alpha = 0.5, type = "t", level = 0.4) +
  stat_ellipse(aes(fill = NULL), colour = "black", linetype = 2, size = 0.2) +
  geom_point(
    size = 3,
    stroke = 0.4,
    alpha = 0.75,
    shape = 21
  ) +
  geom_text_repel(
    data = filter(gg_pca, ID == "2_0_L" |
                    ID == "2_2_S" |
                    ID == "2_6_M" |
                    ID == "2_12_L" |
                    ID == "2_13_S" |
                    ID == "2_14_M"),
    aes(label = phenotype),
    # min.segment.length = 0,
    family = "Helvetica",
    size = ggtext_size
  ) +
  labs(
    x = "PC 1",
    y = "PC 2"
  ) +
  scale_fill_viridis_d() +
  theme_leo() 

pdf("rosette_PCA.pdf", width = onecol, height = onecol)
pca_plot
dev.off()
