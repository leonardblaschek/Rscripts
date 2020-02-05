library(tidyverse)

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
        size = 0.25,
        lineend = "square",
        color = "black"
      ),
      axis.title = element_text(size = 6),
      axis.text.x = element_text(
        colour = "black", # flipped coords
        margin = margin(1, 1, 1, 1),
        size = 6
      ),
      axis.text.y = element_text(
        colour = "black",
        angle = 0,
        vjust = 0.5,
        hjust = 1,
        size = 6,
        margin = margin(1, 1, 1, 1)
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", size = 0.25),
      panel.spacing = unit(1.5, "mm"),
      legend.position = "none",
      legend.text = element_text(size = 6),
      legend.key.height = unit(4, "mm"),
      complete = TRUE
    )
}

py <- read_csv("/home/leonard/Documents/Uni/Phloroglucinol/GCHO_content.csv") %>%
  mutate(genotype = ordered(genotype, levels = c(
    "Col-0",
    "4cl1",
    "4cl2",
    "4cl1x4cl2",
    "ccoaomt1",
    "fah1",
    "omt1",
    "ccr1-3",
    "cad4",
    "cad5",
    "cad4xcad5",
    "Poplar",
    "Spruce"
  )),
  Species = c(rep("At", 11), "Pt", "Pa"))

pdf("GCHO_content.pdf", width = 3.4, height = 1.75)
ggplot(py, aes(x = genotype, y = mean, ymin = mean - st_dev, ymax = mean + st_dev)) +
  geom_col(aes(fill = Species),
           colour = NA) +
  scale_x_discrete(
    labels = c(
      expression(paste(italic("A. thaliana"), " WT")),
      expression(italic("4cl1")),
      expression(italic("4cl2")),
      expression(paste(italic("4cl1"), "x", italic("4cl2"))),
      expression(italic("ccoaomt1")),
      expression(italic("fah1")),
      expression(italic("omt1")),
      expression(italic("ccr1")),
      expression(italic("cad4")),
      expression(italic("cad5")),
      expression(paste(italic("cad4"), "x", italic("cad5"))),
      expression(paste(italic("Poplar"), " WT")),
      expression(paste(italic("Spruce"), " WT"))
    )
  ) +
  scale_fill_manual(values = c( "At" = "#04253a", "Pt" = "#836616", "Pa" = "#833216")) +
  geom_errorbar(width = 0.2, size = 0.25) +
  labs(y = "Lignin coniferaldehyde [%]") +
  theme_leo() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1
    )
  )
dev.off()