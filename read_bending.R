library(tidyverse)

#### fah1 ####
bending_files <-
  list.files(
    path = "/home/leonard/Documents/Uni/PhD/IRX/Bending/fah1_raw/",
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

bending_data <- lapply(bending_files[-1], read_plus) %>%
  bind_rows()

bending_data <- bending_data %>%
  bind_rows(fahA) %>%
  separate(filename, into = c("genotype", "replicate"), sep = " ", extra = "merge") %>%
  mutate(replicate = str_remove(replicate, "Air"),
         replicate = str_remove(replicate, ".txt"),
         replicate = str_remove(replicate, "\\d\\dcm"),
         replicate = str_remove(replicate, " "),
         replicate = as.factor(replicate),
         section = as.factor(section),
         experiment = factor("fah1", levels = c("4cl", "fah1", "cad")),
         genotype = ordered(recode(genotype, "Fah1" = "fah1"), levels = c("WT", "fah1", "cad4xcad5", "4cl1x4cl2x4cl3"))) %>%
  mutate_if(is.character, ~ as.numeric(str_replace(., fixed(","), fixed("."))))

read_stiffness <- function(flnm) {
  read_delim(flnm,
             delim = ";",
             skip = 1,
             col_names = FALSE
  ) %>%
    select(-X1) %>%
    filter(str_detect(X2, 'Fah1|WT')) %>%
    rename("name" = X2,
           "modulus" = X3,
           "max.stress" = X4,
           "AUC" = X5,
           "N.max.stress" = X6) %>%
    mutate(name = str_remove(name, " Air"),
           experiment = factor("fah1", levels = c("4cl", "fah1", "cad"))) %>%
    separate(name, into = c("genotype", "replicate", "sample"), extra = "merge") %>%
    mutate(genotype = ordered(recode(genotype, "Fah1" = "fah1"), levels = c("WT", "fah1", "cad4xcad5", "4cl1x4cl2x4cl3")),
           replicate = as.factor(replicate),
           sample = as.factor(sample)) %>%
    mutate_if(is.character, ~ as.numeric(str_replace(., fixed(","), fixed("."))))
}

stiffness_fah1 <- lapply(bending_files, read_stiffness) %>%
  bind_rows()

#### cad4xcad5 ####
bending_files <-
  list.files(
    path = "/home/leonard/Documents/Uni/PhD/IRX/Bending/cad_raw/",
    pattern = "*.txt",
    recursive = TRUE,
    full.names = TRUE
  )

read_plus <- function(flnm) {
  read_delim(flnm,
             delim = ";",
             skip = grep("Raw Data", readLines(flnm)),
             col_types = "ccccc",
             col_names = FALSE
  ) %>%
    # fill(X1) %>%
    filter(X1 != "Time" & X1 != "(s)") %>%
    mutate(
      filename = basename(flnm),
      filename = str_remove(filename, ".txt"),
      genotype = "cad4xcad5"
    ) %>%
    rename(
      "time" = X1,
      "displacement" = X2,
      "force" = X3,
      "strain" = X4,
      "stress" = X5
    ) %>%
    separate(filename, into = c("section", "replicate"), sep = "_")
}

bending_data_cad <- lapply(bending_files, read_plus) %>%
  bind_rows() %>%
  mutate(replicate = as.factor(replicate),
         section = as.factor(section),
         experiment = factor("cad", levels = c("4cl", "fah1", "cad")),
         genotype = factor(genotype, levels = c("WT", "fah1", "cad4xcad5", "4cl1x4cl2x4cl3"))) %>%
  mutate_if(is.character, ~ as.numeric(str_replace(., fixed(","), fixed("."))))

read_stiffness <- function(flnm) {
  read_delim(flnm,
             delim = ";",
             skip = 1,
             col_names = FALSE
  ) %>%
    select(-X1) %>%
    filter(str_detect(X2, 'top|mid|Base')) %>%
    rename("name" = X2,
           "modulus" = X3,
           "max.stress" = X4) %>%
    mutate(experiment = factor("cad", levels = c("4cl", "fah1", "cad")),
           genotype = factor("cad4xcad5", levels = c("WT", "fah1", "cad4xcad5", "4cl1x4cl2x4cl3"))) %>%
    separate(name, into = c("replicate", "sample"), extra = "merge") %>%
    mutate(replicate = as.factor(replicate),
           sample = as.factor(sample)) %>%
    mutate_if(is.character, ~ as.numeric(str_replace(., fixed(","), fixed("."))))
}

bending_files_stiff <- bending_files %>%
  str_subset(pattern = fixed("_1."))

stiffness_cad <- lapply(bending_files_stiff, read_stiffness) %>%
  bind_rows()

#### 4cl1x2x3 ####

bending_4cl <- read_csv("/home/leonard/Documents/Uni/PhD/IRX/Bending/bending_curve.csv") %>%
  mutate(genotype = ordered(genotype, levels = c("WT", "fah1", "cad4xcad5", "4cl1x4cl2x4cl3")),
         replicate = as.factor(replicate),
         experiment = factor("4cl", levels = c("4cl", "fah1", "cad")),
         section = factor("base"))

stiffness_4cl <- read_csv("/home/leonard/Documents/Uni/PhD/IRX/Bending/stiffness_4cl.csv") %>%
  mutate(genotype = factor(genotype, levels = c("WT", "fah1", "cad4xcad5", "4cl1x4cl2x4cl3")),
         experiment = factor("4cl", levels = c("4cl", "fah1", "cad")),
         sample = as.factor(sample))

#### figures ####
bending_data_complete <- bind_rows(bending_data, bending_data_cad) %>%
  bind_rows(bending_4cl)

bending_pre <- bending_data_complete %>%
  group_by(genotype, replicate, section) %>%
  summarise(break_point = strain[stress == max(stress)],
            strength = max(stress[strain < 3]))

stiffness_complete <- stiffness_fah1 %>%
  bind_rows(stiffness_cad, stiffness_4cl) %>%
  select(-AUC, -N.max.stress) %>%
  pivot_longer(cols = c(modulus, max.stress), names_to = "variable", values_to = "value")

pdf("break_point.pdf", height = onecol / 2, width = onecol / 2)
ggplot(bending_pre, aes(x = genotype, y = break_point)) +
  geom_quasirandom(size = 0.75,
                   alpha = 0.5,
                   shape = 16) +
  geom_violin(alpha = 0.5,
              colour = "black",
              size = 0.2,
              draw_quantiles = 0.5) +
  theme_leo()
dev.off()

curves <- ggplot(
  data = bending_data_complete,
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
  scale_colour_manual(values = c("black",pal_ostwald_disc)) +
  theme_leo() +
  xlim(0, 3) +
  labs(
    x = "Strain [%]",
    y = "Stress [MPa]"
  ) +
  facet_wrap(~genotype, nrow = 1)

pdf("bending_curves_fah.pdf", width = onecol, height = onecol / 2)
curves
dev.off()

flex_violin <- function(x) {
  letters <- letter_groups(stiffness_complete %>% filter(variable == x),
                           value,
                           genotype,
                           "tukey",
                           print_position = "below",
                           print_adjust = 0.5
  )
  
  ggplot(
    data = stiffness_complete %>% filter(variable == x),
    aes(x = genotype, y = value)
  ) +
    # geom_line(aes(group = interaction(cell.type, variable, cell, measurement)),
    #           size = 0.2,
    #           alpha = 0.25) +
    geom_quasirandom(
      aes(fill = genotype),
      dodge.width = 0.6,
      shape = 21,
      # width = 0.1,
      # colour = NA,
      alpha = 0.75,
      size = 1,
      stroke = 0.2
    ) +
    # geom_boxplot(
    #   aes(group = interaction(cell.type, variable), fill = cell.type),
    #   # fill = "white",
    #   alpha = 0.5,
    #   width = 0.6,
    #   position = position_dodge2(width = 1),
    #   outlier.alpha = 0,
    #   lwd = 0.25,
    #   fatten = 1,
    #   width = 0.25
    # ) +
  geom_violin(
    # aes(group = interaction(days, cell.type),
    # fill = cell.type
    # ),
    draw_quantiles = 0.5,
    fill = "white",
    colour = "black",
    alpha = 0.85,
    width = 0.2,
    position = position_dodge(width = 0.6),
    size = 0.2,
    scale = "width"
  ) +
    geom_text(
    data = letters,
    aes(label = groups),
    size = ggtext_size,
    family = "Helvetica",
    position = position_dodge(width = 0.6)
  ) +
    scale_x_discrete(labels = c("WT", 
                                expression(italic("fah1")),
                                expression(italic("cad4x5")),
                                expression(italic("4cl1x2x3")))) +
    scale_fill_manual(values = c(pal_ostwald_disc[1], 
                                 pal_ostwald_disc[3], 
                                 pal_ostwald_disc[3], 
                                 pal_ostwald_disc[3])) +
    # scale_fill_gradientn(colours = pal_ostwald_cont, trans = "reverse")
    expand_limits(y = 0) +
    scale_y_continuous(expand = expand_scale(mult = c(0.1, 0.05))) +
    theme_leo() +
    theme(
      text = element_text(family = "Helvetica"),
      legend.position = "none",
      axis.title.x = element_blank()
    )
}


violin_strength <- flex_violin("max.stress") +
  labs(y = "Strength [MPa]")

pdf("bending_strength.pdf", width = onecol  / 2, height = onecol / 2)
violin_strength
dev.off()


violin_stiffness <- flex_violin("modulus") +
  labs(y = "Stiffness [MPa]")

pdf("bending_stiffness.pdf", width = onecol  / 2, height = onecol / 2)
violin_stiffness
dev.off()

top <- plot_grid(curves)
bottom <- plot_grid(violin_stiffness,
                    violin_strength,
                    ncol = 2,
                    labels = c("B", "C"),
                    label_fontfamily = "Helvetica",
                    label_size = 10,
                    vjust = 1,
                    hjust = 0,
                    rel_widths = c(1, 1))

pdf("IRX_figure5.pdf", width = onecol, height = onecol / 2)
plot_grid(top,
          bottom,
          ncol = 1,
          labels = c("A", ""),
          label_fontfamily = "Helvetica",
          label_size = 10,
          vjust = 1,
          hjust = 0)
dev.off()
          