library(tidyverse)
library(broom)
library(Hmisc)
library(PerformanceAnalytics)
library(showtext)
library(cocor)

#### import Helvetica Neue ####
font_add("Helvetica", regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf",
         bolditalic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-BdIt.otf")
showtext_auto()

#### Arabidopsis main figure ####
#load data
fw_data <- read_csv("/home/leonard/Documents/Uni/Phloroglucinol/foodweb.csv") %>%
  mutate(
    cell.type = recode(cell.type, "V" = "MX"),
    adj.cell.type = recode(adj.cell.type, "V" = "MX"),
    od = ODx255 / 255 # Pixel values of the measured image were multiplied by 255 for better visualisation
  )
 
#average by individual plant
fw_avg <- fw_data %>%
  group_by(genotype, cell.type, adj.cell.type, replicate) %>%
  summarise(od = mean(od, na.rm = TRUE))

fw_spread <- fw_avg %>%
  unite(cell.wall, cell.type, adj.cell.type) %>%
  spread(cell.wall, od)

#data frame of occuring combinations of adjacent cell types
fw_comb <- read_csv("/home/leonard/Documents/Uni/Phloroglucinol/foodweb_comb_at.csv")

#create function for calculation of delta-r, significance and relative impact
corr.comp <- function(data, x, y) {
  coop <- cor.test(data[[paste0(x, "_", y)]], data[[paste0(y, "_", y)]])
  ind <- cor.test(data[[paste0(x, "_", x)]], data[[paste0(y, "_", y)]])
  dep <- cor.test(data[[paste0(x, "_", x)]], data[[paste0(x, "_", y)]])
  tibble(estimate = coop$estimate - ind$estimate,
         p.value = cocor.dep.groups.overlap(r.jk=coop$estimate, 
                           r.jh=ind$estimate, 
                           r.kh=dep$estimate, 
                           n=length(data[[paste0(x, "_", y)]]), 
                           alternative="two.sided", 
                           alpha=0.05, 
                           conf.level=0.95, 
                           null.value=0)@williams1959[["p.value"]],
         rel_effect = mean(data[[paste0(x, "_", y)]], na.rm = T) / mean(data[[paste0(x, "_", x)]], na.rm = T)
  )
}

#map corr.comp over the occurring combinations of adjacent cell types
fw_corrs <- pmap(fw_comb, ~ corr.comp(fw_spread, .x, .y)) %>%
  map_df(as_tibble) %>%
  bind_cols(fw_comb, .) %>%
  mutate(sig = case_when(p.value < 0.05 ~ "yes",
                         TRUE ~"no"))

#plot colour coded relatie impact for use in the network
rel_effects <- ggplot(fw_corrs) +
  geom_point(aes(x = cell.type,
                 y = adj.cell.type,
                 fill = rel_effect,
                 colour = sig),
             shape = 21,
             size = 15,
             stroke = 2) +
  geom_text(aes(label = round(estimate, digits = 2),
                x = cell.type,
                y = adj.cell.type)) +
  scale_colour_manual(values = c("yes" = "darkgreen", "no" = "white")) +
  scale_fill_distiller(palette = "RdBu",
                       limits = c(0.7, 1.3),
                       direction = 1) +
  # theme_leo() +
  theme(legend.position = "bottom")

pdf("fw_reworked_colours.pdf")
rel_effects
dev.off()

#### Arabidopsis supplemental figures ####

fw_matrix <- ungroup(fw_spread) %>%
  select(-genotype, -replicate) %>%
  as.matrix(.)

pdf("corrplots.pdf", width = 12, height = 12)
chart.Correlation(fw_matrix, histogram=TRUE, pch=21)
dev.off()

fw_impact <- fw_spread %>%
  ungroup() %>%
  mutate(genotype = recode(genotype, "col-0" = "WT")) %>%
  group_by(genotype, replicate) %>%
  summarise(
    "CB -> XF" = XF_CB / XF_XF,
    "MX -> XF" = XF_MX / XF_XF,
    "IF -> XF" = XF_IF / XF_XF,
    "PA -> PX" = PX_PA / PX_PX,
    "MX -> PX" = PX_MX / PX_PX,
    "LP -> IF" = IF_LP / IF_IF,
    "CB -> IF" = IF_CB / IF_IF,
    "XF -> IF" = IF_XF / IF_IF,
    "IF -> LP" = LP_IF / LP_LP,
    "CB -> MX" = MX_CB / MX_MX,
    "PA -> MX" = MX_PA / MX_MX,
    "XF -> MX" = MX_XF / MX_MX,
    "PX -> MX" = MX_PX / MX_MX
  ) %>%
  gather(-genotype, -replicate, key = "relationship", value = "od")

impact <-
  ggplot(data = fw_impact, aes(x = ordered(
    genotype,
    levels = rev(c(
      "WT",
      "4cl1",
      "4cl2",
      "4cl1x4cl2",
      "ccoaomt1",
      "fah1",
      "omt1",
      "ccr1",
      "cad4",
      "cad5",
      "cad4xcad5"
    ))
  ),
  y = od,
  fill = od)) + 
  geom_hline(
    yintercept = 1, 
    linetype = 2, 
    size = 0.5
  ) +
  geom_jitter(size = 3, shape = 21, width = 0.05, alpha = 0.75) + 
  scale_fill_gradientn(
    colours = c("#2166ac", "#67a9cf", "#d1e5f0", "white", "#fddbc7", "#ef8a62", "#b2182b"), 
    limits = c(0, 2),
    na.value = "#b2182b"
  ) +
  scale_x_discrete(
    labels = c(
      "WT",
      expression(italic("4cl1")),
      expression(italic("4cl2")),
      expression(paste(italic("4cl1"), "x", italic("4cl2"))),
      expression(italic("ccoaomt1")),
      expression(italic("fah1")),
      expression(italic("omt1")),
      expression(italic("ccr1")),
      expression(italic("cad4")),
      expression(italic("cad5")),
      expression(paste(italic("cad4"), "x", italic("cad5")))
    )
  ) +
  scale_y_continuous(limits = c(0, 4.5)) +
  theme_minimal() + 
  labs(y = "Relative Impact") +
  theme(
    text = element_text(size = 14, family = "Helvetica"),
    axis.ticks = element_line(
      size = 0.25,
      lineend = "square",
      color = "black"
    ),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      colour = "black"
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      fill = NA, 
      color = "black", 
      size = 0.25
    ),
    panel.spacing.x = unit(1.5, "mm"),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    strip.text = element_text(
      vjust = 0.1,
      hjust = 0,
      face = "italic"
    )
  ) +
  facet_wrap( ~ relationship, nrow = 1) +
  coord_flip()

pdf("fw_impact.pdf", width = 15, height = 4)
impact
dev.off()

#### Poplar main figures ####
#load data
poplar <-
  read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol/poplar_foodweb.csv")
poplar <- poplar[, -16]

# calculate distance to the cambium reference line 
# 5.9 is the number of pixels per µm
# X and Y are already measured according to scale, the ref values are given in pixels (see Fiji macro)
poplar$Distance <-
  apply(poplar[, c("X", "Y", "ref.x1", "ref.x2", "ref.y1", "ref.y2")],
        1 ,
        function(x) {
          a <- c(x[1], x[2])
          b <- c((x[3] / 5.9), (x[5] / 5.9))
          c <- c((x[4] / 5.9), (x[6] / 5.9))
          v1 <- b - c
          v2 <- a - b
          m <- cbind(v1, v2)
          d <- abs(det(m)) / sqrt(sum(v1 * v1))
          d
        })

poplar$diff <- poplar$OD.stained - poplar$OD.unstained
poplar_bg <-
  subset(
    poplar,
    cell.type == "CB" & adj.cell.type == "CB",
    select = c("genotype", "replicate", "technical", "diff")
  )
colnames(poplar_bg)[4] <- "OD.bg"

poplar_bg <- poplar_bg %>%
  group_by(genotype, replicate, technical) %>%
  mutate(OD.bg = mean(OD.bg, na.rm = TRUE))

poplar <-
  full_join(poplar,
            unique(poplar_bg),
            by = c("genotype", "replicate", "technical"))
poplar$diff.adj <- poplar$diff - poplar$OD.bg

#bin the measurements by distance from the cambium
poplar_bin <- poplar %>%
  mutate(bin = cut(
    Distance,
    breaks = c(-Inf, 50, 100, Inf),
    labels = c("I", "II", "III")
  ))

poplar_bin_count <- poplar_bin %>%
  group_by(genotype, replicate, cell.type, adj.cell.type, bin) %>%
  summarise(count = n())

poplar_bin_pre <- poplar_bin %>%
  filter(cell.type %in% c("V", "F", "R"), adj.cell.type %in% c("V", "F", "R")) %>%
  group_by(genotype, bin, cell.type, adj.cell.type, replicate) %>%
  summarise(od = mean(diff.adj))

poplar_bin_spread <- poplar_bin_pre %>%
  unite(cell.wall, cell.type, adj.cell.type) %>%
  spread(cell.wall, od) %>%
  group_by(bin) %>%
  group_split()

#data frame of occuring combinations of adjacent cell types
fw_comb <- read_csv("/home/leonard/Documents/Uni/Phloroglucinol/foodweb_comb_pt.csv")

#map corr.comp over the occurring combinations of adjacent cell types in each bin
fw_corrs_1 <- pmap(fw_comb, ~ corr.comp(poplar_bin_spread[[1]], .x, .y)) %>%
  map_df(as_tibble) %>%
  bind_cols(fw_comb, .) %>%
  mutate(sig = case_when(p.value < 0.05 ~ "yes",
                         TRUE ~"no"),
         bin = "1")

fw_corrs_2 <- pmap(fw_comb, ~ corr.comp(poplar_bin_spread[[2]], .x, .y)) %>%
  map_df(as_tibble) %>%
  bind_cols(fw_comb, .) %>%
  mutate(sig = case_when(p.value < 0.05 ~ "yes",
                         TRUE ~"no"),
         bin = "2")

fw_corrs_3 <- pmap(fw_comb, ~ corr.comp(poplar_bin_spread[[3]], .x, .y)) %>%
  map_df(as_tibble) %>%
  bind_cols(fw_comb, .) %>%
  mutate(sig = case_when(p.value < 0.05 ~ "yes",
                         TRUE ~"no"),
         bin = "3")

#bind data frames together and plot colour coded relative impacts
fw_corrs_pt <- rbind(fw_corrs_1, fw_corrs_2, fw_corrs_3)

rel_effects <- ggplot(fw_corrs_pt) +
  geom_point(aes(x = cell.type,
                 y = adj.cell.type,
                 fill = rel_effect,
                 colour = sig),
             shape = 21,
             size = 15,
             stroke = 2) +
  geom_text(aes(label = round(estimate, digits = 2),
                x = cell.type,
                y = adj.cell.type)) +
  scale_colour_manual(values = c("yes" = "darkgreen", "no" = "white")) +
  scale_fill_distiller(palette = "RdBu",
                       limits = c(0, 2),
                       direction = 1, 
                       na.value = "#2166ac") +
  # scale_fill_gradientn(colours = pal_rdbu, limits = c(0.7, 1.3)) +
  # theme_leo() +
  theme(legend.position = "bottom") +
  facet_wrap(~ bin, ncol = 3)

pdf("fw_reworked_colours_pt.pdf")
rel_effects
dev.off()

#### Poplar supplemental figures ####

poplar.matrix.I <- poplar_bin_spread[[1]] %>%
  select(-genotype, -replicate, -bin) %>%
  as.matrix()

poplar.matrix.II <- poplar_bin_spread[[2]] %>%
  select(-genotype, -replicate, -bin) %>%
  as.matrix()

poplar.matrix.III <- poplar_bin_spread[[3]] %>%
  select(-genotype, -replicate, -bin) %>%
  as.matrix()

pdf("corrplots_pt.pdf", width = 6, height = 6)
chart.Correlation(poplar.matrix.I, histogram=TRUE, pch=21)
chart.Correlation(poplar.matrix.II, histogram=TRUE, pch=21)
chart.Correlation(poplar.matrix.III, histogram=TRUE, pch=21)
dev.off()

poplar_plot <- function(x) {
fw_impact_pt <- poplar_bin_spread[[x]] %>%
  ungroup() %>%
  group_by(genotype, replicate) %>%
  summarise(
    "V -> F" = F_V / F_F,
    "V -> R" = R_V / R_R,
    "R -> F" = F_R / F_F,
    "R -> V" = V_R / V_V,
    "F -> R" = R_F / R_R,
    "F -> V" = V_F / V_V,
  ) %>%
  gather(-genotype, -replicate, key = "relationship", value = "od")

print(filter(fw_impact_pt, od > 10))

  ggplot(data = fw_impact_pt, aes(x = ordered(
    genotype,
    levels = rev(c(
      "WT",
      "c4h",
      "ccr"
    ))
  ),
  y = od,
  fill = od)) + 
  geom_hline(
    yintercept = 1, 
    linetype = 2, 
    size = 0.5
  ) +
  geom_jitter(size = 5, shape = 21, width = 0.05, alpha = 0.75) + 
  scale_fill_gradientn(
    colours = c("#2166ac", "#67a9cf", "#d1e5f0", "white", "#fddbc7", "#ef8a62", "#b2182b"), 
    limits = c(0, 2),
    na.value = "#b2182b"
  ) +
  scale_x_discrete(
    labels = rev(c(
      "WT",
      expression(paste(italic("C4H"), "-RNAi")),
      expression(paste(italic("CCR"), "-RNAi"))
    ))
  ) +
  scale_y_continuous(limits = c(0, 10)) +
  theme_minimal() + 
  labs(y = "Relative Impact") +
  theme(
    text = element_text(size = 14, family = "Helvetica"),
    axis.ticks = element_line(
      size = 0.25,
      lineend = "square",
      color = "black"
    ),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      colour = "black"
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      fill = NA, 
      color = "black", 
      size = 0.25
    ),
    # panel.spacing.x = unit(1.5, "mm"),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    strip.text = element_text(
      vjust = 0.1,
      hjust = 0,
      face = "italic"
    )
  ) +
  facet_wrap( ~ relationship, nrow = 3) +
  coord_flip()
}

pdf("fw_impact_pt.pdf", width = 6, height = 6)
poplar_plot(1)
poplar_plot(2)
poplar_plot(3)
dev.off()
