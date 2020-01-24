library(baseline)
library(interactions)
library(effects)
library(grid)
library(png)
library(magick)
library(agricolae)
library(ggthemes)
library(showtext)
library(zoo)
library(cowplot)
library(kableExtra)
library(tidyverse)
library(ggrepel)
library(ggbeeswarm)
library(ggridges)
library(sf)
library(piecewiseSEM)
library(tukeygrps)
library(scales)


pal_ostwald_cont <- c(
  "#155DA7",
  "#0A75B9",
  # "#0C89C9",
  "#FED32F",
  # "#F8A63A",
  "#EF663A",
  "#ED4137"
)
pal_ostwald_cont2 <- c(
  "#275d95",
  "#286e9b",
  "#e8c245",
  "#d37456",
  "#d25952"
)
pal_ostwald_disc <- c(
  "#275d95",
  "#e8c245",
  "#d25952"
)
pal_ostwald_disc2 <- c(
  "#8fab1d",
  "#2d7d73",
  "#1d566f"
)
pal_ostwald_disc3 <- c(
  "#275d95",
  "#275d95",
  "#e8c245",
  "#d25952",
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
pal_ylgnbu <- c(
  "#c7e9b4",
  "#7fcdbb",
  "#41b6c4",
  "#1d91c0",
  "#225ea8",
  "#253494",
  "#081d58"
)


#### import CooperHewitt ####
font_add(
  "Helvetica",
  regular = "/OFL_fonts/IBMPlexSans-Light.otf",
  italic = "/prop_fonts/01. Helvetica   [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
  bold = "/OFL_fonts/IBMPlexSans-SemiBold.otf",
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
      panel.border = element_blank(),
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

scale_max <- function(x) {
  x / max(x)
}

ggtext_size <- 6 / (14 / 5)
twocol <- 18 / 2.54
onehalfcol <- 14 / 2.54
onecol <- 9 / 2.54

#### WT spectrum ####
rmn_files <-
  list.files(
    path = "/home/leonard/Documents/Uni/PhD/Raman/measurements/",
    pattern = "*.CSV",
    recursive = TRUE,
    full.names = TRUE
  )

read_plus <- function(flnm) {
  read_csv(flnm,
    comment = "#",
    col_names = FALSE,
    skip = 0,
    cols_only(
      X1 = col_number(),
      X2 = col_number()
    )
  ) %>%
    mutate(filename = flnm)
}

rmn_data <- lapply(rmn_files, read_plus) %>%
  bind_rows()

# rmn_data <- rmn_data %>%
#   select(c(1:3))

rmn_data$filename <- basename(rmn_data$filename)

# rmn_data$filename <-
#   str_replace(rmn_data$filename, fixed("fah 1"), "fah1")

rmn_data <- rmn_data %>%
  distinct(X1, filename, .keep_all = TRUE) %>%
  separate(
    filename,
    into = c("cell.type", "replicate"),
    sep = "[XhFi]",
    extra = "merge"
  ) %>%
  rename("wavenumber" = X1, "intensity" = X2) %>%
  mutate(
    cell.type = recode(cell.type,
      "Proc-S" = "Vascular bundle",
      "Proc-pt" = "Pith",
      "Proc-I" = "Interfascicular fibre",
      "Proc-Ep" = "Epidermis"
    ),
    replicate = str_remove(replicate, ".CSV")
  )

rmn_data$wavenumber <- round(rmn_data$wavenumber, digits = 0)

rmn_data_wide <- rmn_data %>%
  group_by(cell.type, replicate) %>%
  mutate(group_id = row_number()) %>%
  spread(wavenumber, intensity) %>%
  select(-group_id) %>%
  summarise_all(funs(sum(., na.rm = TRUE)))

rmn_data_mat <- as.matrix(rmn_data_wide[, -c(1:4)])
rownames(rmn_data_mat) <- paste(rmn_data_wide$cell.type,
  rmn_data_wide$replicate,
  sep = "_"
)

rmn_data_correction <- baseline.als(rmn_data_mat,
  lambda = 5,
  p = 0.01,
  maxit = 100
)
rmn_data_all <- as_tibble(rmn_data_correction$corrected, rownames = NA)

rmn_data_all <- rmn_data_all %>%
  rownames_to_column() %>%
  separate(rowname, sep = "_", into = c("cell.type", "replicate"), remove = TRUE) %>%
  gather(key = "wavenumber", value = "corrected.intensity", -c(cell.type, replicate)) %>%
  mutate(wavenumber = str_replace(wavenumber, pattern = fixed("X."), replacement = "-")) %>%
  mutate(wavenumber = str_replace(wavenumber, pattern = fixed("X"), replacement = "")) %>%
  mutate(wavenumber = as.numeric(wavenumber)) %>%
  group_by(cell.type, replicate) %>%
  inner_join(., rmn_data, by = c("cell.type", "replicate", "wavenumber")) %>%
  filter(wavenumber %in% c(300:1800)) %>%
  mutate("total_auc" = MESS::auc(wavenumber, corrected.intensity))

# rmn_data_all$genotype <- ordered(rmn_data_all$genotype,
#                                          levels = c(
#                                            "Col-0",
#                                            "4cl1",
#                                            "4cl2",
#                                            "4cl1x4cl2",
#                                            "ccoaomt1",
#                                            "fah1",
#                                            "omt1",
#                                            "ccr1-3",
#                                            "ccr1xfah1",
#                                            "cad4",
#                                            "cad5",
#                                            "cad4xcad5"
#                                          ))
write_csv(rmn_data_all, "raman_spectra_cell_types.csv")

rmn_data_all_avg <- rmn_data_all %>%
  group_by(cell.type, wavenumber) %>%
  summarise(
    intensity = mean(corrected.intensity),
    auc = mean(total_auc),
    auc_sd = sd(total_auc)
  ) %>%
  ungroup() %>%
  mutate(
    cell.type = as.factor(cell.type),
    fill_col = case_when(
      wavenumber %in% c(1550:1640) ~ "Lignin",
      wavenumber %in% c(1080:1115) ~ "Cellulose",
      TRUE ~ "none"
    )
  )

spectra <- ggplot(data = rmn_data_all_avg %>% filter(cell.type == "Vascular bundle")) +
  scale_fill_manual(values = c("Cellulose" = "#275d95", "Lignin" = "#d25952", "none" = "#0000ff00")) +
  geom_ridgeline_gradient(aes(
    x = wavenumber,
    y = 0,
    height = intensity,
    group = 1,
    fill = fill_col
  ),
  min_height = -40,
  size = 0.2
  ) +
  geom_label_repel(
    data = filter(rmn_data_all_avg %>% filter(cell.type == "Vascular bundle" & wavenumber %in% c(1660, 1620, 1340))),
    aes(
      x = wavenumber,
      label = c("S-units", "Aldehydes", "G-units"),
      y = intensity
    ),
    min.segment.length = 0,
    segment.size = 0.2,
    nudge_x = -100,
    nudge_y = 200,
    label.size = NA,
    family = "Helvetica",
    size = ggtext_size,
    fill = NA
  ) +
  annotate(geom = "text", 
           x = 1099, y = 500, hjust = 0, label = "Cellulose", family = "Helvetica", 
           size = ggtext_size, colour = "#275d95", fontface = "bold") +
  annotate(geom = "text", 
           x = 1593, y = 1986, hjust = 0, label = "Lignin", family = "Helvetica", 
           size = ggtext_size, colour = "#d25952", fontface = "bold") +
  # geom_text(
  #   data = filter(rmn_data_all_avg, wavenumber == 1699), aes(
  #     x = wavenumber,
  #     y = cell.type,
  #     label = cell.type
  #   ),
  #   hjust = 0,
  #   nudge_y = -0.1,
  #   family = "Helvetica",
  #   size = ggtext_size
  # ) +
  # geom_text(
  #   data = filter(rmn_data_all_avg, wavenumber == 403), aes(
  #     x = wavenumber,
  #     y = cell.type,
  #     label = paste("Total area under curve = ",
  #                   round(auc / 1000, digits = 0),
  #                   " AU ± ",
  #                   round(auc_sd / 1000, digits = 0))
  #   ),
  #   hjust = 1,
  #   nudge_y = -0.1,
  #   family = "Helvetica",
  #   size = ggtext_size
  # ) +
  # scale_fill_manual(values = c("Cellulose" = "#8fab1d", "Lignin" = "#2d7d73"), na.translate = F) +
  # scale_fill_viridis_d() +
  scale_x_reverse(expand = expand_scale(mult = 0)) +
  labs(
    x = expression(paste("Wavenumber [cm"^-1, "]")),
    y = "Intensity"
  ) +
  theme_leo() +
  theme(
    plot.margin = unit(c(2, 0, 2, 2), "mm"),
    # axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
pdf("WT_spectra.pdf", width = 4, height = 1.5)
spectra
dev.off()

#### Monomer/DHP spectra ####

rmn_model_files <-
  list.files(
    path = "/home/leonard/Documents/Uni/PhD/Raman/Zoltan/",
    pattern = ".csv",
    recursive = FALSE,
    full.names = TRUE
  )

read_plus <- function(flnm) {
  read_csv(flnm,
           comment = "#",
           col_names = FALSE,
           skip = 5,
           cols_only(
             X1 = col_number(),
             X2 = col_number()
           )
  ) %>%
    mutate(filename = flnm)
}

rmn_models <- lapply(rmn_model_files, read_plus) %>%
  bind_rows()

# rmn_models_wide <- rmn_models %>%
#   group_by(filename) %>%
#   rename("wavenumber" = "X1",
#          "intensity" = "X2") %>%
#   mutate(group_id = row_number()) %>%
#   spread(wavenumber, intensity) %>%
#   select(-group_id) %>%
#   summarise_all(funs(sum(., na.rm = TRUE)))

rmn_models_nest <- rmn_models %>%
  group_by(filename) %>%
  rename("wavenumber" = "X1",
         "intensity" = "X2") %>%
  nest()

pvt <- function(df) pivot_wider(data = df, names_from = wavenumber, values_from = intensity)

bl_corr <- function(df) {
  correction <- baseline(as.matrix(df),
                         method = "als",
                         lambda = 5,
                         p = 0.01,
                         maxit = 100
  )
  as_tibble((getCorrected(correction)))
}

rmn_models_corrected <- rmn_models_nest %>%
  mutate(wide = map(data, pvt),
         correction = map(wide, bl_corr)) %>%
  select(-wide, -data) %>%
  unnest_wider(correction) %>%
  pivot_longer(cols = -filename, names_to = "wavenumber", values_to = "intensity") %>%
  drop_na() %>%
  ungroup()

rmn_models <- rmn_models_corrected %>%
  mutate(filename = str_remove(basename(filename), ".\\.csv$")) %>%
  separate(filename, into = c("state", "unit", "gamma"), remove = FALSE) %>%
  # rename(
  #   "wavenumber" = "X1",
  #   "intensity" = "X2"
  # ) %>%
  mutate(
    wavenumber = round(as.numeric(wavenumber), digits = 0),
    unit = recode(unit,
                  "g" = "G",
                  "s" = "S"
    ),
    state = recode(state,
                   "mono" = "Monomeric",
                   "poly" = "Polymerised")
  ) %>%
  filter(wavenumber %in% c(300:1700)) %>%
  group_by(state, unit, gamma) %>%
  mutate(
    AUC = MESS::auc(wavenumber, intensity),
    # AUC_1600_prop = MESS::auc(wavenumber, intensity, from = 1580, to = 1615) / AUC,
    # scaled = intensity
    scaled = intensity / AUC
  ) %>%
  ungroup() %>%
  mutate(scaled = scale_max(scaled))

rmn_models <- rmn_models %>%
  ungroup() %>%
  mutate(filename = ordered(filename,
                            levels =
                              rev(c(
                                "poly-g-oh",
                                "mono-g-oh",
                                "poly-g-cho",
                                "mono-g-cho",
                                "mono-g-cooh",
                                "mono-s-oh",
                                "poly-s-cho",
                                "mono-s-cho",
                                "mono-s-cooh",
                                "mono-h-cooh",
                                "mono-c-cooh",
                                "mono-gsat-oh",
                                "mono-5h-cooh"
                              ))
  ), gamma = ordered(gamma, levels = rev(c("oh", "cho", "cooh"))))


spectra_plot <- ggplot(filter(rmn_models, 
                              filename %in% c("mono-g-oh", "mono-g-cho"))) +
  geom_text_repel(
    data = tibble(
      wavenumber = c(1656, 1620, 1600, 1334, 1276),
      unit = rep("G", len = 5),
      gamma = rep("CHO", 5),
      scaled = rep(1, 5)
    ),
    aes(
      label = wavenumber,
      x = wavenumber,
      y = scaled
    ),
    angle = 90,
    hjust = 0.5,
    direction = "x",
    min.segment.length = 0,
    segment.size = 0.4,
    nudge_y = -0.25,
    segment.colour = "grey85",
    family = "Helvetica",
    size = ggtext_size
  ) +
  annotate("segment",
           x = c(1656, 1620, 1600, 1334, 1276),
           xend = c(1656, 1620, 1600, 1334, 1276),
           y = 1,
           yend = Inf,
           colour = "grey85",
           size = 0.4
  ) +
  geom_ridgeline(aes(
    x = wavenumber,
    y = gamma,
    height = rollmean(scaled, 3, na.pad = TRUE),
    fill = state,
    colour = state,
  ),
  alpha = 0.25,
  size = 0.2,
  # colour = NA,
  # fill = NA,
  scale = 2.5,
  min_height = -0.05
  ) +
  # geom_ridgeline(
  #   data = tibble(
  #     wavenumber = rep(300:1700, 3),
  #     gamma = rep(c("cooh", "oh", "cooh"), each = 1401),
  #     unit = rep(c("G", "S", "S"), each = 1401),
  #     state = "poly",
  #     scaled = 0
  #   ),
  #   aes(
  #     x = wavenumber,
  #     y = gamma,
  #     height = rollmean(scaled, 3, na.pad = TRUE),
  #     # fill = state,
  #     colour = state
  #   ),
  #   alpha = 0.5,
  #   size = 0.2,
  #   # colour = NA,
  #   fill = NA,
  #   linetype = 2,
  #   scale = 2.5,
  #   min_height = -0.002
  # ) +
  scale_fill_manual(values = c("Polymerised" = "#d73027", "Monomeric" = "#4575b4"), name = "State") +
  scale_colour_manual(values = c("Polymerised" = "#d73027", "Monomeric" = "#4575b4"), name = "State") +
  scale_x_reverse(
    limits = c(1700, 300),
    # sec.axis = sec_axis(~., breaks = unique(peak_pres$peak))
  ) +
  scale_y_discrete(
    labels = c(
      expression(bold("G")[CHO]),
      expression(bold("G")[CHOH])
    )
  ) +
  labs(x = expression(paste("Wavenumber [cm"^-1, "]"))) +
  theme_leo() +
  theme(
    # panel.spacing = unit(5, "mm"),
    legend.position = c(0.75, 1),
    legend.justification = c(0, 1),
    legend.key.height = unit(2, "mm"),
    legend.key.width = unit(2, "mm"),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 8),
    axis.text.x.top = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
  ) 

pdf("models_pres_mono.pdf", width = 4, height = 2.5)
spectra_plot
dev.off()

#### Signal treatment ####
rmn_files <-
  list.files(
    path = "/home/leonard/Documents/Uni/PhD/IRX/RAMAN/nuoendagula/",
    pattern = "*.txt",
    recursive = TRUE,
    full.names = TRUE
  )

read_plus <- function(flnm) {
  read_tsv(flnm,
           comment = "#",
           col_names = FALSE,
           skip = 25,
           cols_only(
             X1 = col_number(),
             X2 = col_number()
           )
  ) %>%
    mutate(
      filename = flnm,
      exp_time = readLines(flnm) %>%
        str_split("\t") %>%
        grep("Exposure Time: .* s", ., value = T) %>%
        str_remove_all("#Exposure Time: | s") %>%
        as.numeric(),
      pwr = readLines(flnm) %>%
        str_split("\t") %>%
        grep("Excitation Power: .* mW", ., value = T) %>%
        str_remove_all("#Excitation Power: | mW") %>%
        as.numeric(),
      nd_fltr = 1 - (readLines(flnm) %>%
                       str_split("\t") %>%
                       grep("% \\([0-9]{3}/255\\)", ., value = T) %>%
                       str_remove_all("#ND Filter : | % \\([0-9]{3}/255\\)") %>%
                       as.numeric() / 100)
    )
}

rmn_data <- lapply(rmn_files, read_plus) %>%
  bind_rows()

# rmn_data <- rmn_data %>%
#   select(c(1:3))

rmn_data$filename <- basename(rmn_data$filename)

rmn_data$filename <-
  str_replace(rmn_data$filename, fixed("fah 1"), "fah1")

rmn_data <- rmn_data %>%
  distinct(X1, filename, .keep_all = TRUE) %>%
  separate(
    filename,
    into = c("genotype", "replicate", "cell.type", "technical"),
    sep = "\\s*#|-|\\s",
    extra = "merge"
  ) %>%
  rename("wavenumber" = X1, "intensity" = X2) %>%
  mutate(
    genotype = recode(
      genotype,
      "4cl1＆2" = "4cl1x4cl2",
      "cad4＆5" = "cad4xcad5",
      "ccoaomt" = "ccoaomt1",
      "ccr1" = "ccr1-3",
      "col.o" = "Col-0",
      "Col.0" = "Col-0",
      "ccr1＆fah1" = "ccr1xfah1"
    ),
    cell.type = recode(cell.type,
                       "px" = "PX",
                       "SX" = "SMX",
                       "MX" = "PMX"
    ),
    replicate = str_extract(
      replicate,
      "\\d"
    ),
    technical = str_extract(
      technical,
      "(\\d)+"
    ),
    wavenumber = round(wavenumber, digits = 0)
  ) %>%
  filter(cell.type != "？？") %>%
  group_by(genotype, cell.type, replicate, technical) %>%
  mutate(intensity = intensity / (exp_time * pwr * nd_fltr)) %>%
  select(-exp_time, -pwr, -nd_fltr)

rmn_data_wide <- rmn_data %>%
  group_by(genotype, cell.type, replicate, technical) %>%
  mutate(group_id = row_number()) %>%
  spread(wavenumber, intensity) %>%
  select(-group_id) %>%
  summarise_all(funs(sum(., na.rm = TRUE)))

rmn_data_mat <- as.matrix(rmn_data_wide[, -c(1:4)])
rownames(rmn_data_mat) <- paste(rmn_data_wide$genotype,
                                rmn_data_wide$cell.type,
                                rmn_data_wide$replicate,
                                rmn_data_wide$technical,
                                sep = "_"
)

rmn_data_correction <- baseline.als(rmn_data_mat,
                                    lambda = 5,
                                    p = 0.01,
                                    maxit = 100
)
rmn_data_corrected <- as_tibble(rmn_data_correction$corrected, rownames = NA)

rmn_data_corrected <- rmn_data_corrected %>%
  rownames_to_column() %>%
  separate(rowname, sep = "_", into = c("genotype", "cell.type", "replicate", "technical"), remove = TRUE) %>%
  gather(key = "wavenumber", value = "corrected.intensity", -c(genotype, cell.type, replicate, technical)) %>%
  mutate(wavenumber = str_replace(wavenumber, pattern = fixed("X."), replacement = "-")) %>%
  mutate(wavenumber = str_replace(wavenumber, pattern = fixed("X"), replacement = "")) %>%
  mutate(wavenumber = as.numeric(wavenumber)) %>%
  group_by(genotype, cell.type, replicate, technical) %>%
  inner_join(., rmn_data, by = c("genotype", "cell.type", "replicate", "technical", "wavenumber")) %>%
  ungroup() %>%
  filter(genotype != "ccr1xfah1" &
           wavenumber %in% c(300:1700) &
           cell.type %in% c("IF", "PMX")) %>%
  anti_join(tibble(genotype = "cad4xcad5",
                   cell.type = "PMX",
                   replicate = "4",
                   technical = "2")) %>%
  mutate(
    genotype = ordered(genotype,
                       levels = c(
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
                         "cad4xcad5"
                       )
    ),
    cell.type = recode(cell.type, "PMX" = "VB")
  )

# write_csv(rmn_data_corrected, "rmn_data_corrected.csv")
# rmn_data_corrected <- read_csv("rmn_data_corrected.csv", col_types = "ccccnnn") %>%
#   mutate(
#     genotype = ordered(genotype,
#       levels = c(
#         "Col-0",
#         "4cl1",
#         "4cl2",
#         "4cl1x4cl2",
#         "ccoaomt1",
#         "fah1",
#         "omt1",
#         "ccr1-3",
#         "cad4",
#         "cad5",
#         "cad4xcad5"
#       )
#     )
#   )

# lig_peak <- rmn_data_corrected %>%
#   group_by(genotype, cell.type, replicate, technical) %>%
#   filter(wavenumber %in% c(1580:1620)) %>%
#   mutate(
#     peak.pos = wavenumber[which.max(corrected.intensity)],
#     peak.height = corrected.intensity[which(wavenumber == peak.pos)]
#   ) %>%
#   filter(wavenumber == 1580) %>%
#   select(genotype:technical, peak.pos, peak.height)

rmn_data_corrected <- rmn_data_corrected %>%
  # left_join(lig_peak) %>%
  group_by(genotype, cell.type, replicate, technical) %>%
  mutate(
    AUC = MESS::auc(wavenumber, corrected.intensity),
    AUC_1603 = MESS::auc(wavenumber, corrected.intensity, from = 1580, to = 1640),
    AUC_1334 = MESS::auc(wavenumber, corrected.intensity, from = 1320, to = 1358),
    AUC_1099 = MESS::auc(wavenumber, corrected.intensity, from = 1081, to = 1111),
    AUC_381 = MESS::auc(wavenumber, corrected.intensity, from = 376, to = 398),
    scaled = corrected.intensity / AUC,
    scaled_lig = corrected.intensity / (
      corrected.intensity[wavenumber == 1603] + 2 * corrected.intensity[wavenumber == 1334]
    ),
    scaled_max = corrected.intensity / corrected.intensity[wavenumber == 1603],
    scaled_cellu = corrected.intensity / corrected.intensity[wavenumber == 381]
  )

rmn_data_pre <- rmn_data_corrected %>%
  group_by(genotype, cell.type, replicate, wavenumber) %>%
  summarise(
    samples.pre = length(corrected.intensity),
    mean.raw = mean(corrected.intensity, na.rm = TRUE),
    mean.scaled = mean(scaled, na.rm = TRUE),
    mean.scaled_lig = mean(scaled_lig, na.rm = TRUE),
    mean.scaled_max = mean(scaled_max, na.rm = TRUE),
    mean.scaled_cellu = mean(scaled_cellu, na.rm = TRUE)
    # mean.lig_peak = mean(lig_peak, na.rm = TRUE)
  )

rmn_data_avg <- rmn_data_pre %>%
  group_by(genotype, cell.type, wavenumber) %>%
  summarise(
    plants = paste("plants: ", length(mean.scaled)),
    bundles = ifelse(length(mean.scaled) > 1,
                     paste("cells/plant: ", min(samples.pre), "—", max(samples.pre), sep = ""),
                     paste("cells/plant: ", min(samples.pre))
    ),
    mean.raw = mean(mean.raw, na.rm = TRUE),
    mean.scaled = mean(mean.scaled, na.rm = TRUE),
    mean.scaled_lig = mean(mean.scaled_lig, na.rm = TRUE),
    mean.scaled_max = mean(mean.scaled_max, na.rm = TRUE),
    mean.scaled_cellu = mean(mean.scaled_cellu, na.rm = TRUE)
    # mean.lig_peak = mean(mean.lig_peak, na.rm = TRUE)
  )

spectra_plot <- ggplot(filter(rmn_data_pre,
                              genotype == "Col-0" &
                                cell.type == "IF")) +
  # geom_text_repel(
  #   data = tibble(
  #     wavenumber = c(1656, 1620, 1600, 1334, 1276),
  #     unit = rep("G", len = 5),
  #     gamma = rep("CHO", 5),
  #     scaled = rep(1, 5)
  #   ),
  #   aes(
  #     label = wavenumber,
  #     x = wavenumber,
  #     y = scaled
  #   ),
  #   angle = 90,
  #   hjust = 0.5,
  #   direction = "x",
  #   min.segment.length = 0,
  #   segment.size = 0.4,
  #   nudge_y = -0.25,
  #   segment.colour = "grey85",
  #   family = "Helvetica",
  #   size = ggtext_size
  # ) +
  # annotate("segment",
  #          x = c(1656, 1620, 1600, 1334, 1276),
  #          xend = c(1656, 1620, 1600, 1334, 1276),
  #          y = 1,
  #          yend = Inf,
  #          colour = "grey85",
  #          size = 0.4
  # ) +
  geom_ridgeline(aes(
    x = wavenumber,
    height = rollmean(mean.raw, 3, na.pad = TRUE),
    fill = cell.type,
    y = 0,
    group = replicate
  ),
  alpha = 0.25,
  size = 0.2,
  # colour = NA,
  # fill = NA,
  scale = 2.5,
  min_height = -15
  ) +
  # geom_ridgeline(
  #   data = tibble(
  #     wavenumber = rep(300:1700, 3),
  #     gamma = rep(c("cooh", "oh", "cooh"), each = 1401),
  #     unit = rep(c("G", "S", "S"), each = 1401),
  #     state = "poly",
  #     scaled = 0
  #   ),
  #   aes(
  #     x = wavenumber,
  #     y = gamma,
#     height = rollmean(scaled, 3, na.pad = TRUE),
#     # fill = state,
#     colour = state
#   ),
#   alpha = 0.5,
#   size = 0.2,
#   # colour = NA,
#   fill = NA,
#   linetype = 2,
#   scale = 2.5,
#   min_height = -0.002
# ) +
scale_fill_manual(values =  "#d73027") +
# scale_fill_brewer(palette = "Reds", direction = -1) +
  scale_x_reverse(
    limits = c(1700, 300),
    # sec.axis = sec_axis(~., breaks = unique(peak_pres$peak))
  ) +
  # scale_y_discrete(
  #   labels = c(
  #     expression(bold("G")[CHO]),
  #     expression(bold("G")[CHOH])
  #   )
  # ) +
  labs(x = expression(paste("Wavenumber [cm"^-1, "]"))) +
  theme_leo() +
  theme(
    # panel.spacing = unit(5, "mm"),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    strip.text = element_text(size = 8),
    axis.text.x.top = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
  ) 

pdf("singal_treat_before_IF.pdf", width = 4, height = 2.5)
spectra_plot
dev.off()

#### correlations ####

selected_peaks <- read_csv("/home/leonard/Documents/Uni/PhD/Raman/raman_peak_detection_2.csv")

py <- read.csv("/home/leonard/Documents/Uni/PhD/Raman/pyrolysis_2018_amount.csv")

py$Coniferaldehyde <- py$X22
py$Vanillin <- py$X18
py$Sinapaldehyde <- py$X32
py$Syringaldehyde <- py$X29
py$GOH <- ifelse(
  py$variable == "mean",
  py$X17 + py$X19 + py$X23,
  sqrt(py$X17^2 + py$X19^2 + py$X23^2)
)
py$`S-OH` <- py$X33
py$`S-CHO` <- py$X32
py$`S-CHOs` <- ifelse(
  py$variable == "mean",
  py$X32 + py$X30 + py$X29,
  sqrt(py$X32^2 + py$X30^2 + py$X29^2)
)
py$`SnoCHO` <- ifelse(
  py$variable == "mean",
  py$X26 + py$X27 + py$X31 + py$X28 + py$X33,
  sqrt(py$X27^2 + py$X31^2 + py$X28^2 + py$X26^2 + py$X33^2)
)
py$GCHO <- ifelse(
  py$variable == "mean",
  py$X18 + py$X22,
  sqrt(py$X18^2 + py$X22^2)
)
py$Aldehydes <-
  ifelse(
    py$variable == "mean",
    py$Coniferaldehyde + py$Sinapaldehyde + py$Vanillin + py$Syringaldehyde,
    sqrt(
      py$Coniferaldehyde^2 + py$Sinapaldehyde^2 + py$Vanillin^2 + py$Syringaldehyde^2
    )
  )
py$H <-
  ifelse(py$variable == "mean", rowSums(py[, 7:9]), sqrt(rowSums(py[, 7:9]^
                                                                   2)))
py$G <-
  ifelse(py$variable == "mean", rowSums(py[, 7:19]), sqrt(rowSums(py[, 7:19]^
                                                                    2)))
py$S <-
  ifelse(py$variable == "mean", rowSums(py[, 21:28]), sqrt(rowSums(py[, 21:28]^
                                                                     2)))
py$SnoOH <-
  ifelse(py$variable == "mean", rowSums(py[, 21:27]), sqrt(rowSums(py[, 21:27]^
                                                                     2)))
py$Lignin <-
  ifelse(py$variable == "mean", rowSums(py[, 3:28]), sqrt(rowSums(py[, 3:28]^
                                                                    2)))
py$`GCHO/G` <- py$GCHO / py$G
py$`GOH/G` <- py$GOH / py$G
py$`GOH/lig` <- py$GOH / py$Lignin
py$`GCHO/GOH` <- py$Coniferaldehyde / py$GOH
py$`GCHO/lig` <- py$GCHO / py$Lignin
py$`Aldehydes/lig` <- py$Aldehydes / py$Lignin
py$`S/G` <- py$S / py$G
py$`H/G` <- py$H / py$G
py$`H/S` <- py$H / py$S

py <- py %>%
  select(-starts_with("X")) %>%
  pivot_longer(-c(genotype, variable), names_to = "unit", values_to = "value") %>%
  pivot_wider(id_cols = c("genotype", "unit"), names_from = variable, values_from = value)

py_rmn <- rmn_data_corrected %>%
  filter(wavenumber %in% selected_peaks$wavenumber) %>%
  group_by(cell.type, genotype, wavenumber, replicate, technical) %>%
  pivot_wider(
    id_cols = c(genotype, cell.type, replicate, technical), names_from = wavenumber,
    values_from = scaled
  ) %>%
  group_by(cell.type, genotype, replicate, technical) %>%
  transmute(
    "CLB" = (`1603` / 2) + `1334`,
    "1334/1603" = `1334` / `1603`,
    "1664/1603" = `1664` / `1603`,
    "1664/lig" = `1664` / ((`1603` / 2) + `1334`),
    "1334/1276" = `1334` / `1276`,
    "1625/1603" = `1625` / `1603`,
    "1625/1664" = `1625` / `1664`,
    "1625/lig" = `1625` / ((`1603` / 2) + `1334`),
    "1625" = `1625`,
    "1334" = `1334`,
    "1603" = `1603`
  ) %>%
  pivot_longer(cols = -c(genotype, cell.type, replicate, technical), names_to = "wavenumber", values_to = "scaled") %>%
  group_by(genotype, wavenumber, cell.type) %>%
  summarise(
    mean.scaled = mean(scaled),
    sd.scaled = sd(scaled) / mean.scaled
  ) %>%
  group_by(genotype, wavenumber) %>%
  summarise(
    total.mean.scaled = mean(mean.scaled),
    total.sd.scaled = total.mean.scaled * sqrt(sum(sd.scaled^2))
  )

py_rmn_corr <- left_join(py, py_rmn) %>%
  mutate(WT = case_when(
    genotype == "Col-0" ~ "yes",
    TRUE ~ "no"
  ))

py_rmn_corr_mat <- py_rmn_corr %>%
  # filter(variable == "mean") %>%
  group_by(unit, wavenumber) %>%
  nest() %>%
  mutate(
    cor = map(data, ~ cor.test(.x$mean, .x$total.mean.scaled)),
    tidied = map(cor, tidy)
  ) %>%
  unnest(tidied)

broom::glance(lm(mean ~ total.mean.scaled, data = filter(py_rmn_corr, unit == "Lignin" &
                                                           wavenumber == "CLB")))

# broom::tidy(lm(Lignin ~ `1603` + `1334`, data = filter(py_rmn_corr, variable == "mean")))





py_plot <- function(x, y, z = "bork", col = "white") {
  tidycor <- broom::tidy(cor.test(
    as.matrix(filter(py_rmn_corr, unit == x &
                       wavenumber == y &
                       genotype != z)[, "mean"]),
    as.matrix(filter(py_rmn_corr, unit == x &
                       wavenumber == y &
                       genotype != z)[, "total.mean.scaled"])
  ))
  print(tidycor)
  
  regr <- broom::glance(lm(mean ~ total.mean.scaled, data = filter(py_rmn_corr, unit == x &
                                                                     wavenumber == y &
                                                                     genotype != z)))
  
  ggplot(
    data = filter(py_rmn_corr, unit == x & wavenumber == y & genotype != z),
    aes(x = mean, y = total.mean.scaled)
  ) +
    geom_rect(
      data = tidycor,
      xmax = Inf,
      ymax = Inf,
      xmin = -Inf,
      ymin = -Inf,
      aes(
        fill = estimate,
        x = 0,
        y = 0
      )
    ) +
    geom_errorbar(aes(
      ymin = total.mean.scaled - total.sd.scaled,
      ymax = total.mean.scaled + total.sd.scaled
    ),
    size = 0.2,
    width = 0,
    alpha = 0.75
    ) +
    # geom_errorbarh(aes(xmin = mean - sd,
    #                   xmax = mean + sd),
    #               size = 0.2,
    #               height = 0) +
    geom_point(aes(colour = WT),
               size = 1,
               alpha = 0.75,
               shape = 16
    ) +
    geom_smooth(
      method = "lm",
      fullrange = T,
      size = 0.4,
      se = F,
      linetype = 1,
      colour = col
    ) +
    scale_colour_manual(values = c("yes" = "white", no = "black")) +
    scale_fill_distiller(
      palette = "RdYlBu",
      na.value = "white",
      name = "Pearson's r",
      limits = c(-1, 1)
    ) +
    scale_x_continuous(limits = c(0, NA)) +
    scale_y_continuous(limits = c(0, NA)) +
    annotate("text",
             label = paste("R² = ", round(regr$adj.r.squared, digits = 2)),
             x = 0,
             y = 0,
             size = ggtext_size,
             family = "Helvetica",
             hjust = 0,
             vjust = 0,
             colour = col
    ) +
    theme_leo() +
    theme(
      axis.text.y = element_text(
        angle = 90,
        hjust = 0.5
      ),
      plot.margin = unit(c(3, 2, 3, 3), "mm")
    )
  # geom_label_repel(aes(label = genotype),
  #                  size = 1,
  #                  segment.size = 0.2)
}

py_plot_lig <- py_plot("Lignin", "CLB") +
  labs(
    x = expression(paste("Lignin * total pyrolysate"^-1)),
    y = expression(paste("CLB [AU * AUC"^-1, "]"))
  ) +
  scale_y_continuous(breaks = c(0, 0.003, 0.006), labels = c("0", "0.003", "0.006")) +
  scale_x_continuous(labels = c("0", "0.05", "0.1", "0.15"))

py_plot_sg <- py_plot("S/G", "1334/1276") +
  labs(
    x = expression(paste(bold("S")[R], " * ",bold("G")[R]^{
      -1
    })),
    y = expression(paste("1334 cm"^-1, "/ 1276 cm"^-1))
  ) +
  scale_y_continuous(breaks = c(0, 1, 2, 3)) +
  scale_x_continuous(labels = c("0", "0.1", "0.2", "0.3", "0.4"))

py_plot_cho <- py_plot("GCHO/lig", "1625/lig") +
  labs(
    x = expression(paste(bold("G")[CHO], " * lignin"^-1)),
    y = expression(paste("1625 cm"^-1, "/ CLB"))
  ) +
  # scale_y_continuous(breaks = c(0, 0.4, 0.8), labels = c("0", "0.4", "0.8")) +
  scale_x_continuous(labels = c("0", "0.1", "0.2", "0.3", "0.4")) +
  geom_point(
    data = filter(py_rmn_corr, unit == "GCHO/lig" & wavenumber == "1625/lig" & genotype == "cad4xcad5"),
    aes(colour = WT),
    size = 1,
    alpha = 0.75,
    shape = 16
  ) +
  geom_errorbar(
    data = filter(py_rmn_corr, unit == "GCHO/lig" & wavenumber == "1625/lig" & genotype == "cad4xcad5"),
    aes(
      ymin = total.mean.scaled - total.sd.scaled,
      ymax = total.mean.scaled + total.sd.scaled
    ),
    size = 0.2,
    width = 0,
    alpha = 0.75
  ) +
  geom_smooth(
    data = filter(py_rmn_corr, unit == "GCHO/lig" & wavenumber == "1625/lig" & genotype != "cad4xcad5"),
    method = "lm",
    colour = "black",
    fullrange = T,
    size = 0.4,
    se = F,
    linetype = 2
  )
# annotate("text",
#          label = "R² = −0.12",
#          x = 0,
#          y = 0.1,
#          size = ggtext_size,
#          family = "Helvetica",
#          hjust = 0,
#          vjust = 0,
#             colour = "blue")

py_plot_goh <- py_plot("GOH/lig", "1664/lig") +
  labs(
    x = expression(paste(bold("G")[CHOH], " * lignin"^-1)),
    y = expression(paste("1664 cm"^-1, "/CLB"))
  )
# scale_y_continuous(breaks = c(0, 1, 2, 3))

py_plot_snooh <- py_plot("SnoOH", "1334") +
  labs(
    x = expression(paste("S-lignin " - " S-OH")),
    y = expression(paste("1334 cm"^-1))
  ) +
  scale_y_continuous(breaks = c(0, 0.0015, 0.003), labels = c("0", "0.0015", "0.003")) +
  scale_x_continuous(breaks = c(0, 0.01, 0.02), labels = c("0", "0.01", "0.02"))

py_plot_soh <- py_plot("S-OH", "1334") +
  labs(
    x = expression(paste("S-OH")),
    y = expression(paste("1334 cm"^-1))
  ) +
  scale_y_continuous(breaks = c(0, 0.0015, 0.003), labels = c("0", "0.0015", "0.003")) +
  scale_x_continuous(breaks = c(0, 0.005, 0.01), labels = c("0", "0.005", "0.01"))

py_plot_scho <- py_plot("S-CHO", "1334") +
  labs(
    x = expression(paste("S-CHO")),
    y = expression(paste("1334 cm"^-1))
  )
# scale_y_continuous(breaks = c(0, 0.0015, 0.003), labels = c("0", "0.0015", "0.003")) +
# scale_x_continuous(breaks = c(0, 0.005, 0.01), labels = c("0", "0.005", "0.01"))

py_plot_snocho <- py_plot("SnoCHO", "1334") +
  labs(
    x = expression(paste("S" - "S-CHO")),
    y = expression(paste("1334 cm"^-1))
  )

pdf("Raman_correlations_pres.pdf", width = 5, height = 1.3)
plot_grid(py_plot_lig,
          py_plot_sg,
          py_plot_goh,
          py_plot_cho,
          nrow = 1)
dev.off()

#### Vessel type Raman ####

rmn_irx <- read.csv("file:///home/leonard/Documents/Uni/PhD/IRX/raman_IRX_revisited.csv") %>%
  select(genotype:Perim., Circ., Round, Height, Width) %>%
  filter(object == 2) %>%
  select(-object)

rmn_nb <- read.csv("file:///home/leonard/Documents/Uni/PhD/IRX/raman_IRX.csv") %>%
  select(genotype:n_f)
rmn_irx <- left_join(rmn_irx, rmn_nb)
rmn_irx$replicate <- as.character(rmn_irx$replicate)
rmn_irx$technical <- as.character(rmn_irx$technical)

raman_convex <- read.csv("/home/leonard/Documents/Uni/PhD/IRX/raman_irx_convex.csv")

raman_convex$File <-
  str_replace(raman_convex$File, fixed("fah 1"), "fah1")

raman_convex <- raman_convex %>%
  select(-X) %>%
  separate(
    File,
    into = c("genotype", "replicate", "cell.type", "technical"),
    sep = "\\s*#|-|\\s",
    extra = "merge"
  ) %>%
  mutate(
    genotype = recode(
      genotype,
      "4cl1＆2" = "4cl1x4cl2",
      "cad4＆5" = "cad4xcad5",
      "ccoaomt" = "ccoaomt1",
      "ccr1" = "ccr1-3",
      "col.o" = "Col-0",
      "Col.0" = "Col-0",
      "ccr1＆fah1" = "ccr1xfah1"
    ),
    cell.type = recode(cell.type,
                       "px" = "PX",
                       "SX" = "SMX",
                       "MX" = "PMX"
    ),
    replicate = str_extract(
      replicate,
      "\\d"
    ),
    technical = str_extract(
      technical,
      "(\\d)+"
    )
  ) %>%
  filter(cell.type != "？？")

rmn_irx <- left_join(rmn_irx, raman_convex) %>%
  # group_by(genotype, replicate, technical, cell.type) %>%
  mutate(
    Perim. = Perim. * (5.9 / 3.2), # correct wrong scale while measuring
    Area = Area * (5.9 / 3.2)^2,
    ConvexArea = ConvexArea * (5.9 / 3.2)^2
  )


rmn_data_corrected <- read_csv("rmn_data_corrected.csv", col_types = "ccccnnn") %>%
  mutate(
    genotype = ordered(genotype,
      levels = c(
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
        "cad4xcad5"
      )
    )
  )

rmn_data_corrected <- rmn_data_corrected %>%
  # left_join(lig_peak) %>%
  group_by(genotype, cell.type, replicate, technical) %>%
  mutate(
    AUC = MESS::auc(wavenumber, corrected.intensity),
    scaled = corrected.intensity / AUC,
    scaled_lig = corrected.intensity / (
      (corrected.intensity[wavenumber == 1603] / 2) + corrected.intensity[wavenumber == 1334]
    )
  )

rmn_data_pre <- rmn_data_corrected %>%
  group_by(genotype, cell.type, replicate, wavenumber) %>%
  summarise(
    samples.pre = length(corrected.intensity),
    mean.scaled = mean(scaled, na.rm = TRUE),
    mean.scaled_lig = mean(scaled_lig, na.rm = TRUE),
    # mean.lig_peak = mean(lig_peak, na.rm = TRUE)
  )

rmn_data_avg <- rmn_data_pre %>%
  group_by(genotype, cell.type, wavenumber) %>%
  summarise(
    mean.scaled = mean(mean.scaled, na.rm = TRUE),
    mean.scaled_lig = mean(mean.scaled_lig, na.rm = TRUE),
    # mean.lig_peak = mean(mean.lig_peak, na.rm = TRUE)
  )

rmn_sum <- rmn_data_corrected %>%
  group_by(genotype, cell.type, replicate, technical) %>%
  summarise(
    "total_lignin" = (scaled[wavenumber == 1603] / 2) + scaled[wavenumber == 1334],
    "GOH_rel" = scaled_lig[wavenumber == 1664],
    "GOH.G" = scaled_lig[wavenumber == 1664] / scaled_lig[wavenumber == 1603],
    # "GOH_tot" = scaled[wavenumber == 1664],
    # "G_tot" = scaled[wavenumber == 1603],
    "G_rel" = scaled_lig[wavenumber == 1603],
    "S.G" = scaled[wavenumber == 1334] / scaled[wavenumber == 1276],
    "GCHO.lig" = scaled_lig[wavenumber == 1625],
    "S_rel" = scaled_lig[wavenumber == 1334],
    # "S_tot" = scaled[wavenumber == 1334],
    "cellulose" = scaled[wavenumber == 1099] + scaled[wavenumber == 381],
    "AUC" = mean(AUC)
  )

wiesner.data <-
  read.csv("/home/leonard/Documents/Uni/Phloroglucinol/measurements_revisited.csv",
           skip = 2
  )

wiesner.data$genotype <- recode(wiesner.data$genotype, cad4x5 = "cad4xcad5")

# set cell types according to measurement order
wiesner.data[1:50 + rep(seq(0, (nrow(wiesner.data) - 50), by = 300), each = 50), 4] <-
  "IF"
wiesner.data[51:100 + rep(seq(0, (nrow(wiesner.data) - 50), by = 300), each = 50), 4] <-
  "MX"
wiesner.data[101:150 + rep(seq(0, (nrow(wiesner.data) - 50), by = 300), each = 50), 4] <-
  "XF"
wiesner.data[151:200 + rep(seq(0, (nrow(wiesner.data) - 50), by = 300), each = 50), 4] <-
  "PX"
wiesner.data[201:250 + rep(seq(0, (nrow(wiesner.data) - 50), by = 300), each = 50), 4] <-
  "LP"
wiesner.data[251:300 + rep(seq(0, (nrow(wiesner.data) - 50), by = 300), each = 50), 4] <-
  "PH"

# import SMX measurements
wiesner.data.SMX <- read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol/measurements_SMX.csv",
                             skip = 2
)
wiesner.data.SMX$replicate <- as.factor(wiesner.data.SMX$replicate)

wiesner.data <- rbind(select(wiesner.data, -technical), select(wiesner.data.SMX, -technical))

wiesner.data$genotype <-
  ordered(
    wiesner.data$genotype,
    levels = c(
      "col-0",
      "4cl1",
      "4cl2",
      "4cl1x2",
      "ccoaomt1",
      "fah1",
      "omt1",
      "ccr1-3",
      "ccr1xfah1",
      "cad4",
      "cad5",
      "cad4xcad5"
    )
  )

# calculate the correct hue on the 360 point circular scale
wiesner.data$hue <- ((wiesner.data$h.stained + 128) / 255 * 360)

wiesner.data$replicate <-
  as.factor(as.character(wiesner.data$replicate))

# calculate stained - unstained diff. and adjust for bleaching by subtracting the diff. for the unlignified phloem
wiesner.data$diff <-
  wiesner.data$OD.stained - wiesner.data$OD.unstained
# summary(wiesner.data)

wiesner.data.bg <- wiesner.data %>%
  filter(cell.type == "PH") %>%
  select(1:2, 9) %>%
  ungroup() %>%
  group_by(genotype, replicate) %>%
  summarise(OD.bg = mean(diff, na.rm = TRUE))

wiesner.data <-
  full_join(
    wiesner.data,
    wiesner.data.bg,
    all = TRUE,
    by = c("genotype", "replicate")
  )
wiesner.data$diff.adj <- wiesner.data$diff - wiesner.data$OD.bg
wiesner.data <- subset(wiesner.data, cell.type != "PH")

# average per replicate (for boxplots)
wiesner.data.pre <- wiesner.data %>%
  mutate(
    cell.type = recode(cell.type, "MX" = "PMX"),
    genotype = recode(genotype,
                      "col-0" = "Col-0",
                      "4cl1x2" = "4cl1x4cl2"
    )
  ) %>%
  group_by(genotype, cell.type, replicate) %>%
  summarise(
    mean.hue1 = mean(hue, na.rm = TRUE),
    SD.hue1 = sd(hue, na.rm = TRUE),
    mean.OD1 = mean(diff.adj, na.rm = TRUE),
    SD.OD1 = sd(diff.adj, na.rm = TRUE)
  ) %>%
  select(genotype, cell.type, replicate, mean.OD1)

rmn_irx <- read.csv("file:///home/leonard/Documents/Uni/PhD/IRX/raman_IRX_revisited.csv") %>%
  select(genotype:Perim., Circ., Round, Height, Width) %>%
  filter(object == 2) %>%
  select(-object)

raman.nb <- read.csv("file:///home/leonard/Documents/Uni/PhD/IRX/raman_IRX.csv") %>%
  select(genotype:n_f)
rmn_irx <- left_join(rmn_irx, raman.nb)
rmn_irx$replicate <- as.character(rmn_irx$replicate)
rmn_irx$technical <- as.character(rmn_irx$technical)

raman_convex <- read.csv("/home/leonard/Documents/Uni/PhD/IRX/raman_irx_convex.csv")

raman_convex$File <-
  str_replace(raman_convex$File, fixed("fah 1"), "fah1")

raman_convex <- raman_convex %>%
  select(-X) %>%
  separate(
    File,
    into = c("genotype", "replicate", "cell.type", "technical"),
    sep = "\\s*#|-|\\s",
    extra = "merge"
  ) %>%
  mutate(
    genotype = recode(
      genotype,
      "4cl1＆2" = "4cl1x4cl2",
      "cad4＆5" = "cad4xcad5",
      "ccoaomt" = "ccoaomt1",
      "ccr1" = "ccr1-3",
      "col.o" = "Col-0",
      "Col.0" = "Col-0",
      "ccr1＆fah1" = "ccr1xfah1"
    ),
    cell.type = recode(cell.type,
                       "px" = "PX",
                       "SX" = "SMX",
                       "MX" = "PMX"
    ),
    replicate = str_extract(
      replicate,
      "\\d"
    ),
    technical = str_extract(
      technical,
      "(\\d)+"
    )
  ) %>%
  filter(cell.type != "？？")

rmn_irx <- left_join(rmn_irx, raman_convex)

# write_csv(rmn_irx, "raman_irx_shape.csv")


raman.heights <- read.csv("file:///home/leonard/Documents/Uni/Master/Summer project 16/phenotyping/phenotyping.csv") %>%
  select(genotype, replicate, height) %>%
  mutate(genotype = recode(genotype, "col-0" = "Col-0")) %>%
  rename("Plant.height" = "height") %>%
  rbind(read.csv("file:///home/leonard/Documents/Uni/PhD/IRX/heights_haris.csv"))
raman.heights$replicate <- as.character(raman.heights$replicate)


rmn_irx <- left_join(rmn_irx, rmn_sum)
rmn_irx.switch <- left_join(read_tsv("/home/leonard/Documents/Uni/PhD/IRX/PMX_which_SMX.tsv",
                                     col_types = "cccc"
), rmn_irx) %>%
  mutate(
    technical = str_c(technical, "switched", sep = "_"),
    cell.type = "SMX"
  )
rmn_irx <- rmn_irx %>%
  anti_join(read_tsv("/home/leonard/Documents/Uni/PhD/IRX/PMX_which_SMX.tsv",
                     col_types = "cccc"
  )) %>%
  bind_rows(rmn_irx.switch)
rmn_irx <- left_join(rmn_irx, raman.heights)
rmn_irx <- left_join(rmn_irx, wiesner.data.pre) %>%
  mutate(convexity = Area / ConvexArea) 

rmn_irx_long <- rmn_irx %>%
  pivot_longer(-c(genotype:cell.type), names_to = "variable", values_to = "value") %>%
  ungroup() %>%
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
    "cad4xcad5"
  )),
  cell.type = ordered(cell.type, levels = c("PX", "PMX", "SMX")))



rmn_violin <- function(x) {
  letters <- letter_groups(rmn_irx_long %>% filter(variable == x & genotype %in% c("Col-0", "4cl1x4cl2", "cad4xcad5", "fah1")),
                           value,
                           genotype,
                           "tukey",
                           cell.type,
                           print_position = "below",
                           print_adjust = 1
  ) %>%
    ungroup() %>%
    mutate(cell.type = ordered(cell.type, levels = c("PX", "PMX", "SMX")))
  
  ggplot(
    data = rmn_irx_long %>% filter(variable == x & genotype %in% c("Col-0", "4cl1x4cl2", "cad4xcad5", "fah1")),
    aes(x = genotype, y = value)
  ) +
    geom_quasirandom(
      aes(fill = cell.type),
      dodge.width = 0.6,
      shape = 21,
      # width = 0.1,
      # colour = NA,
      alpha = 0.75,
      size = 2,
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
  geom_violin(aes(
    group = genotype,
    # fill = cell.type
  ),
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
    scale_x_discrete(
      labels = c(
        "Col-0",
        expression(paste(italic("4cl1"), "x", italic("4cl2"))),
        expression(italic("fah1")),
        expression(paste(italic("cad4"), "x", italic("cad5")))
      )
    ) +
    scale_fill_manual(values = pal_flame_disc) +
    expand_limits(y = 0) +
    theme_leo() +
    theme(
      text = element_text(family = "Helvetica"),
      legend.position = "none",
      axis.ticks = element_blank()
      # axis.text.y = element_blank()
    ) +
    facet_wrap(~ cell.type)
}

rmn_lig <- rmn_violin("total_lignin") +
  labs(
    x = "",
    y = expression(paste("Relative lignin"))
  ) +
  theme(axis.text.x = element_blank())

rmn_s.g <- rmn_violin("S.G") +
  labs(
    x = "",
    y = expression(paste("Relative S/G-ratio"))
  ) +
  theme(strip.text = element_blank(),
        axis.text.x = element_blank())

rmn_goh <- rmn_violin("GOH.G") +
  labs(
    x = "",
    y = expression(paste("Relative G"[CHOH]))
  ) +
  theme(strip.text = element_blank())

pdf("raman_overview_pres.pdf", height = 4, width = 6)
plot_grid(rmn_lig,
          rmn_s.g,
          rmn_goh,
          axis = "b",
          align = "hv",
          nrow = 3)
dev.off()

rmn_violin_wiesner <- function(x) {
  letters <- letter_groups(wiesner.data.pre %>% filter(cell.type %in% c("PX", "PMX", "SMX") & genotype %in% c("Col-0", "4cl1x4cl2", "cad4xcad5", "fah1")),
                           mean.OD1,
                           genotype,
                           "tukey",
                           cell.type,
                           print_position = "below",
                           print_adjust = 1
  ) %>%
    ungroup() %>%
    mutate(cell.type = ordered(cell.type, levels = c("PX", "PMX", "SMX")))
  
  ggplot(
    data = wiesner.data.pre %>% 
      filter(cell.type %in% c("PX", "PMX", "SMX") & 
               genotype %in% c("Col-0", "4cl1x4cl2", "cad4xcad5", "fah1")) %>%
      ungroup() %>%
      mutate(cell.type = ordered(cell.type, levels = c("PX", "PMX", "SMX"))),
    aes(x = genotype, y = mean.OD1)
  ) +
    geom_quasirandom(
      aes(fill = cell.type),
      dodge.width = 0.6,
      shape = 21,
      # width = 0.1,
      # colour = NA,
      alpha = 0.75,
      size = 2,
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
  geom_violin(aes(
    group = genotype,
    # fill = cell.type
  ),
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
    scale_x_discrete(
      labels = c(
        "Col-0",
        expression(paste(italic("4cl1"), "x", italic("4cl2"))),
        expression(italic("fah1")),
        expression(paste(italic("cad4"), "x", italic("cad5")))
      )
    ) +
    scale_fill_manual(values = pal_flame_disc) +
    expand_limits(y = 0) +
    theme_leo() +
    theme(
      text = element_text(family = "Helvetica"),
      legend.position = "none",
      axis.ticks = element_blank()
    ) +
    facet_wrap(~ cell.type)
}

rmn_wiesner <- rmn_violin_wiesner("mean.OD1") +
  labs(
    x = "",
    y = expression(paste("Coniferaldehyde"))
  ) 

pdf("wiesner_vessels_pres.pdf", height = 1.4, width = 6)
rmn_wiesner
dev.off()


#### Wiesner poplar ####

pop_irx <- read_csv("/home/leonard/Documents/Uni/PhD/IRX/Poplar/poplar_irx.csv")

#### create reference at the heigth of the cambium for distance calculation ####
pop_irx_ref <- pop_irx %>%
  filter(object == "ref") %>%
  select(c(
    "genotype",
    "replicate",
    "technical",
    "X",
    "Y",
    "Length"
  ))
pop_irx_ref$ref.y1 <- pop_irx_ref$Y
pop_irx_ref$ref.y2 <- pop_irx_ref$Y
pop_irx_ref$ref.x1 <- pop_irx_ref$X - (pop_irx_ref$Length / 2)
pop_irx_ref$ref.x2 <- pop_irx_ref$X + (pop_irx_ref$Length / 2)

#### merge reference into data frame ####
pop_irx <- pop_irx %>%
  filter(object != "ref") %>%
  full_join(select(pop_irx_ref, -X, -Y, -Length),
            by = c("genotype", "replicate", "technical")
  ) %>%
  mutate(genotype = ordered(genotype, levels = c("WT", "c4h", "ccr1")))

#### calculate difference of the centre of each vessel to the cambium ####
pop_irx$Distance <-
  apply(
    pop_irx[, c("X", "Y", "ref.x1", "ref.x2", "ref.y1", "ref.y2")],
    1,
    function(x) {
      a <- c(x[1], x[2])
      b <- c(x[3], x[5])
      c <- c(x[4], x[6])
      v1 <- b - c
      v2 <- a - b
      m <- cbind(v1, v2)
      d <- abs(det(m)) / sqrt(sum(v1 * v1))
      d
    }
  )

pop_irx <- pop_irx %>%
  mutate(bin = cut(
    Distance,
    breaks = c(-Inf, 50, 100, Inf),
    labels = c("I", "II", "III")
  ))

pop_irx <- pop_irx %>%
  mutate(
    convexity = Area / ConvexArea,
    Mean = Mean / 255,
    type = recode(type, "P" = "Primary",
                  "S" = "Secondary"),
    genotype = recode(genotype,
                      "c4h" = "C4H-RNAi",
                      "ccr1" = "CCR1-RNAi")
  )

pop_dist <- ggplot(pop_irx, aes(x = Distance, y = Mean)) +
  geom_point(aes(fill = type),
             shape = 21,
             size = 2,
             alpha = 0.75,
             stroke = 0.2
  ) +
  scale_fill_manual(values = pal_flame_disc[c(2,3)], name = "Vessel type") +
  theme_leo() +
  labs(x = "Distance from cambium [µm]",
       y = "Coniferaldehyde in vessels") +
  geom_smooth(
    method = "lm",
    span = 1,
    size = 0.2,
    colour = "black",
    linetype = 2
  ) +
    theme(axis.ticks = element_blank(),
          legend.position = c(0.9, 0.9)) +
  facet_wrap(~ genotype, nrow = 3)

pdf("pop_dist_pres.pdf", width = 6, height = 3)
pop_dist
dev.off()

