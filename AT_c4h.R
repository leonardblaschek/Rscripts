library(tidyverse)
library(ggridges)
library(baseline)

#### import and tidy poplar Raman data ####
rmn_files <-
  list.files(
    path = "/data/PhD/Raman/C4H/spectra",
    pattern = "*.txt",
    recursive = TRUE,
    full.names = TRUE
  )

read_plus <- function(flnm) {
  read_tsv(flnm,
           comment = "#",
           col_names = c("wavenumber", "intensity"),
           col_types = "cc",
           # locale = locale(decimal_mark = ",")
  ) %>%
    mutate(
      flnm = str_remove(basename(flnm), fixed(".txt"))
      ) %>%
    separate(flnm, into = c("genotype", "replicate", "cell_type", "cell", "layer"), sep = "_") %>%
    mutate(
      wavenumber = as.numeric(str_replace(wavenumber, fixed(","), fixed("."))),
      wavenumber = case_when(genotype == "WT" & replicate == 1 ~ wavenumber,
                             TRUE ~ wavenumber + 30),
      intensity = as.numeric(str_replace(intensity, fixed(","), fixed("."))),
      layer = case_when(is.na(layer) ~ "full",
                        TRUE ~ layer)
    ) %>%
    pivot_wider(
      id_cols = c(genotype, replicate, cell_type, cell, layer),
      names_from = wavenumber,
      values_from = intensity
    ) %>%
    group_by(genotype, replicate, cell_type, cell, layer) %>%
    nest()
}

rmn_c4h <- lapply(rmn_files, read_plus) %>%
  bind_rows()

#### baseline correction ####
rmn_als <- function(x) {
  corrected_spectra <- baseline.als(as.matrix(x),
                                    lambda = 5,
                                    p = 0.01,
                                    maxit = 100
  )
  
  as_tibble(corrected_spectra$corrected, rownames = NA)
}

rmn_c4h_corrected <- rmn_c4h %>%
  mutate(corrected_data = map(data, rmn_als)) %>%
  select(-data) %>%
  unnest(cols = corrected_data) %>%
  pivot_longer(-c(genotype:layer),
               names_to = "wavenumber",
               values_to = "corrected.intensity"
  ) %>%
  drop_na() %>%
  mutate(wavenumber = as.numeric(wavenumber))

#### filtering and scaling ####
rmn_c4h_scaled <- rmn_c4h_corrected %>%
  filter(wavenumber > 200 & wavenumber < 1700) %>%
  group_by(genotype, replicate, cell_type, cell, layer) %>%
  mutate(
    AUC = MESS::auc(wavenumber, corrected.intensity),
    scaled = corrected.intensity / AUC,
    scaled_lig = corrected.intensity / (
      (corrected.intensity[abs(wavenumber - 1603) == min(abs(wavenumber - 1603))] / 2) + corrected.intensity[abs(wavenumber - 1334) == min(abs(wavenumber - 1334))]
    )
  )

rmn_c4h_values <- rmn_c4h_scaled %>%
  group_by(genotype, cell_type, replicate, cell) %>%
  summarise(raw_995 = corrected.intensity[abs(wavenumber - 995) == min(abs(wavenumber - 995))][1],
            scaled_995 = scaled[abs(wavenumber - 995) == min(abs(wavenumber - 995))][1])

rmn_max <- rmn_c4h_scaled %>%
  group_by(genotype, replicate, cell, cell_type) %>%
  summarise (peak_position = wavenumber[which.max(scaled)])

#### export ####
rmn_c4h_raw <- rmn_c4h_scaled %>%
  mutate(wavenumber = round(wavenumber, digits = 0)) %>%
  pivot_wider(id_cols = c(genotype, replicate, cell_type, cell, layer),
              names_from = wavenumber,
              values_from = corrected.intensity)

rmn_c4h_raw <- rmn_c4h_raw %>%
  select(genotype, replicate, cell_type, cell, layer, order(as.numeric(colnames(.))))

rmn_c4h_auc <- rmn_c4h_scaled %>%
  mutate(wavenumber = round(wavenumber, digits = 0)) %>%
  pivot_wider(id_cols = c(genotype, replicate, cell_type, cell, layer),
              names_from = wavenumber,
              values_from = corrected.intensity)

rmn_c4h_auc <- rmn_c4h_auc %>%
  select(genotype, replicate, cell_type, cell, layer, order(as.numeric(colnames(.))))

write_excel_csv(rmn_c4h_values, "Raman_c4h_995.csv")
write_excel_csv(rmn_c4h_raw, "Raman_c4h_raw.csv")
write_excel_csv(rmn_c4h_auc, "Raman_c4h_scaled.csv")

#### plot ####

c4h_celltypes <- ggplot(rmn_c4h_scaled,
                        # aes(x = wavenumber,
                        #     y = mean.intensity,
                        #     colour = cell_type)
                        ) +
  geom_vline(xintercept = 995,
             linetype = 2,
             size = 0.2) +
  geom_ridgeline(
    aes(x = wavenumber,
        y = reorder(cell_type, -corrected.intensity),
        height = corrected.intensity,
        fill = reorder(cell_type, -corrected.intensity),
        group = interaction(genotype, cell_type, replicate, cell)
        ),
    min_height = -200,
    scale = 0.0002,
    size = 0.2,
    alpha = 0.75
  ) +
  # geom_line(size = 0.2,
  #           alpha = 0.75) +
  theme_leo() +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  # scale_colour_brewer(palette = "Dark2") +
  scale_fill_viridis_d() +
  scale_x_continuous(trans = "reverse") +
  labs(x = "Wavenumber",
       y ="Intensity (baseline corrected)") +
  facet_wrap(~genotype)

c4h_celltypes_scaled <- ggplot(rmn_c4h_scaled,
                        # aes(x = wavenumber,
                        #     y = mean.intensity,
                        #     colour = cell_type)
) +
  geom_vline(xintercept = 995,
             linetype = 2,
             size = 0.2) +
  geom_ridgeline(
    aes(x = wavenumber,
        y = reorder(cell_type, -corrected.intensity),
        height = scaled,
        fill = reorder(cell_type, -corrected.intensity),
        group = interaction(genotype, cell_type, replicate, cell)
        ),
    min_height = -0.02,
    scale = 200,
    size = 0.2,
    alpha = 0.75
  ) +
  # geom_line(size = 0.2,
  #           alpha = 0.75) +
  theme_leo() +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  # scale_colour_brewer(palette = "Dark2") +
  scale_fill_viridis_d() +
  scale_x_continuous(trans = "reverse") +
  labs(x = "Wavenumber",
       y ="Intensity (corrected and scaled to AUC)") +
  facet_wrap(~genotype)

layers <- ggplot(rmn_c4h_scaled %>% filter(cell_type == "IF"),
       aes(x = wavenumber,
           y = corrected.intensity,
           colour = layer)) +
  geom_vline(xintercept = 995,
             linetype = 2,
             size = 0.2) +
  geom_line(aes(group = interaction(layer, replicate, cell)),
            size = 0.2,
            alpha = 0.75) +
  theme_leo() +
  theme(legend.position = "bottom") +
  scale_x_continuous(trans = "reverse") +
  scale_colour_manual(values = pal_ostwald_disc) +
  facet_wrap(~genotype)

rmn_box <- function(y) {
  y <- enquo(y)
  letters <- letter_groups(
    rmn_c4h_values %>% filter(cell_type != "PA"),
    !!y,
    genotype,
    "tukey",
    cell_type,
    print_position = "above"
  )
  
  ggplot(rmn_c4h_values %>% filter(cell_type != "PA"),
                     aes(x = genotype,
                         y = !!y)) +
  geom_quasirandom(
    shape = 16,
    alpha = 0.5,
    groupOnX = TRUE,
    width = 0.2) +
  # stat_summary(
  #   fun = mean,
  #   fun.min = mean,
  #   fun.max = mean,
  #   geom = "crossbar",
  #   width = 0.5,
  #   size = 0.5,
  #   colour = pal_ostwald_disc[3],
  #   fatten = 0
  # ) +
  geom_boxplot(
    fill = "white",
    colour = "black",
    alpha = 0.85,
    outlier.alpha = 0,
    width = 0.2,
    size = 0.2,
  ) +
    geom_text(data = letters,
              aes(label = Letters),
              family = "Helvetica",
              size = ggtext_size) +
  # geom_violin(
  #   draw_quantiles = 0.5,
  #   fill = "white",
  #   colour = "black",
  #   alpha = 0.85,
  #   width = 0.2,
  #   size = 0.2,
  #   scale = "width") +
  theme_leo() +
  theme(axis.title.x = element_blank()) +
  facet_wrap(~cell_type, ncol = 5)
}

raw <- rmn_box(raw_995) + labs(y = "Corrected intensity")
scaled <- rmn_box(scaled_995) + labs(y = "AUC-scaled intensity")

patch <-raw + scaled


pdf("rmn_c4h.pdf", width = twocol, height = onecol)
# c4h_celltypes
# layers
patch +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(plot.tag = element_text(size = 10, face = "bold", family = "Helvetica")) &
  plot_layout(ncol = 1)
dev.off()