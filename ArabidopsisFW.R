library(dplyr)
library(broom)
library(purrr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(PerformanceAnalytics)
library(showtext)

#### import Helvetica Neue ####
font_add("Helvetica", regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf",
         bolditalic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-BdIt.otf")
showtext_auto()

#### import data ####
fw.data <-
  read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol/foodweb.csv")

fw.data <-
  mutate(
    fw.data,
    cell.type = recode(cell.type, "V" = "MX"),
    adj.cell.type = recode(adj.cell.type, "V" = "MX")
  )

fw.data$od <- fw.data$ODx255 / 255 # Pixel values of the measured image were multiplied by 255 for better visualisation

fw.avg <- fw.data %>%
  group_by(genotype, cell.type, adj.cell.type, replicate) %>%
  summarise(od = mean(od, na.rm = TRUE))

fw.spread <- fw.avg %>%
  unite(cell.wall, cell.type, adj.cell.type) %>%
  spread(cell.wall, od)

#### calculate relative influence ####
fw.loading <- fw.spread %>%
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

#### calculate linear regressions for the network ####
fw.lm <- ungroup(fw.spread) %>%
  select(-genotype, -replicate) %>%
  do("MX -> PX" = lm(PX_MX ~ MX_MX, data = .),
     "MX -> XF" = lm(XF_MX ~ MX_MX, data = .),
     "PX -> MX" = lm(MX_PX ~ PX_PX, data = .),
     "XF -> MX" = lm(MX_XF ~ XF_XF, data = .),
     "XF -> IF" = lm(IF_XF ~ XF_XF, data = .),
     "IF -> XF" = lm(XF_IF ~ IF_IF, data = .),
     "IF -> LP" = lm(LP_IF ~ IF_IF, data = .),
     "LP -> IF" = lm(IF_LP ~ LP_LP, data = .)
  )
fw.lm <- fw.lm %>%
  gather(key ="relationship", value = "regression") 

fw.lm.tidy <- full_join(
  map_dfr(fw.lm$regression, glance, .id = "relationship"),
  subset(map_dfr(fw.lm$regression, tidy, .id = "relationship"), term != "(Intercept)"),
  by = "relationship")
fw.lm.tidy$relationship <- fw.lm$relationship

#### export regression values for network ####
write.csv(fw.lm.tidy, file = "FW_regressions.csv")

fw.matrix <- ungroup(fw.spread) %>%
  select(-genotype, -replicate) %>%
  as.matrix(.)

#### calculate pearson correlations ####
fw.corr <- rcorr(fw.matrix)

#### plot pearson corelations (Fig. S7) ####
pdf("corrplots.pdf", width = 12, height = 12)
chart.Correlation(fw.matrix, histogram=TRUE, pch=21)
dev.off()

influence_ttest <- fw.spread %>%
  group_by(genotype) %>%
  summarise(
    "CB -> XF" = t.test(XF_CB, XF_XF)$p.value,
    "MX -> XF" = t.test(XF_MX, XF_XF)$p.value,
    "IF -> XF" = t.test(XF_IF, XF_XF)$p.value,
    "PA -> PX" = t.test(PX_PA, PX_PX)$p.value,
    "MX -> PX" = t.test(PX_MX, PX_PX)$p.value,
    "LP -> IF" = t.test(IF_LP, IF_IF)$p.value,
    "CB -> IF" = t.test(IF_CB, IF_IF)$p.value,
    "XF -> IF" = t.test(IF_XF, IF_IF)$p.value,
    "IF -> LP" = t.test(LP_IF, LP_LP)$p.value,
    "CB -> MX" = t.test(MX_CB, MX_MX)$p.value,
    "PA -> MX" = t.test(MX_PA, MX_MX)$p.value,
    "XF -> MX" = t.test(MX_XF, MX_MX)$p.value,
    "PX -> MX" = t.test(MX_PX, MX_MX)$p.value
  ) %>%
  pivot_longer(-genotype, names_to = "relationship", values_to = "p_value") %>%
  filter(p_value < 0.05)

#### plot relative influence (Fig. 8b) ####
interaction <-
  ggplot(data = fw.loading, aes(x = ordered(
    genotype,
    levels = rev(c(
      "col-0",
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
    linetype = 1, 
    size = 0.5
  ) +
  geom_jitter(size = 3, shape = 21, width = 0.05, alpha = 0.75) + 
  scale_fill_gradientn(
    colours = c("#2166ac", "#67a9cf", "#d1e5f0", "white", "#fddbc7", "#ef8a62", "#b2182b"), 
    limits = c(0, 2),
    na.value = "#b2182b"
  ) +
  scale_x_discrete(
    labels = rev(c(
      "Col-0",
      expression(italic("4cl1")),
      expression(italic("4cl2")),
      expression(italic("4cl1x4cl2")),
      expression(italic("ccoaomt1")),
      expression(italic("fah1")),
      expression(italic("omt1")),
      expression(italic("ccr1")),
      expression(italic("cad4")),
      expression(italic("cad5")),
      expression(italic("cad4xcad5"))
    ))
  ) +
  scale_y_continuous(limits = c(0, 4.5)) +
  theme_minimal() + 
  labs(y = "Relative Influence") +
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

pdf("fw_loading.pdf", width = 15, height = 4)
interaction
dev.off()
