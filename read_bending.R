library(tidyverse)

# bending <- read_delim("/home/leonard/Documents/Uni/PhD/IRX/Bending/raw/WT A 26cm.txt",
#                       delim = ";",
#                       skip = grep("Raw Data", readLines("/home/leonard/Documents/Uni/PhD/IRX/Bending/raw/WT A 26cm.txt")),
#                       col_types = "cccccc",
#                       col_names = FALSE) %>%
#   fill(X1) %>%
#   filter(X2 != "Time" & X2 != "(s)")


bending_files <-
  list.files(
    path = "/home/leonard/Documents/Uni/PhD/IRX/Bending/raw/",
    pattern = "*.txt",
    recursive = TRUE,
    full.names = TRUE
  )

read_plus <- function(flnm) {
  read_delim(flnm,
    delim = ";",
    skip = grep("Raw Data", readLines(flnm)),
    col_types = "cccccc",
    col_names = FALSE
  ) %>%
    fill(X1) %>%
    filter(X2 != "Time" & X2 != "(s)") %>%
    mutate(
      filename = basename(flnm)
    ) %>%
    rename(
      "section" = X1,
      "time" = X2,
      "displacement" = X3,
      "force" = X4,
      "strain" = X5,
      "stress" = X6
           )
}

bending_data <- lapply(bending_files[-1], read_plus) %>%
  bind_rows()

fahA <- read_delim(bending_files[1],
                      delim = ";",
                      skip = grep("Raw Data", readLines(bending_files[1])),
                      col_types = "cccccc",
                      col_names = FALSE) %>%
  mutate(filename = basename(bending_files[1])) %>%
  fill(X1) %>%
  filter(X2 != "Displacement" & X2 != "(mm)") %>%
  rename(
    "section" = X1,
    "displacement" = X2,
    "force" = X3,
    "strain" = X4,
    "stress" = X5
  )

bending_data <- bending_data %>%
  bind_rows(fahA) %>%
  separate(filename, into = c("genotype", "replicate"), sep = " ", extra = "merge") %>%
  mutate(replicate = str_remove(replicate, "Air"),
         replicate = str_remove(replicate, ".txt"),
         replicate = str_remove(replicate, "\\d\\dcm"),
         replicate = str_remove(replicate, " "),
         replicate = as.factor(replicate),
         genotype = ordered(recode(genotype, "Fah1" = "fah1"), levels = c("WT", "fah1"))) %>%
  mutate_if(is.character, ~ as.numeric(str_replace(., fixed(","), fixed("."))))

bending_pre <- bending_data %>%
  group_by(genotype, replicate, section) %>%
  summarise(break_point = strain[stress == max(stress)])

pdf("break_point.pdf", height = onecol / 2, width = onecol / 2)
ggplot(bending_pre, aes(x = genotype, y = break_point)) +
  geom_point() +
  geom_boxplot() +
  theme_leo()
dev.off()


curves <- ggplot(
  data = bending_data,
  aes(x = strain, y = stress)
) +
  # geom_vline(xintercept = 0.65,
  #            size = 0.2,
  #            linetype = 2) +
  # annotate("text",
  #          label = "WT breaking point",
  #          family = "Helvetica",
  #          size = ggtext_size,
  #          x = 0.7,
  #          y = 8.5,
  #          hjust = 0,
  #          vjust = 0) +
  geom_line(aes(group = interaction(genotype, replicate, section)),
            size = 0.2,
            alpha = 0.15
  ) +
  geom_smooth(aes(colour = genotype),
              size = 0.4
  ) +
  scale_color_manual(values = pal_ostwald_disc[c(1, 3)]) +
  theme_leo() +
  xlim(0, 3) +
  labs(
    x = "Flexural strain [%]",
    y = "Flexural stress [MPa]"
  ) +
  facet_wrap(~genotype, nrow = 2)

pdf("bending_curves_fah.pdf", width = onecol / 2, height = onecol / 2)
curves
dev.off()