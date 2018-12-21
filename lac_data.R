library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(showtext)
library(reshape2)
library(cowplot)


# import Helvetica Neue
font_add("Helvetica", regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf")
showtext_auto()

lac.data <- read.csv("file:///home/leonard/Dropbox/Review/lac_data.csv")

lac.short <- subset(lac.data, select = c(1,2,8,10:24))
lac.short <- melt(lac.short, c("Kingdom", "Species", "Temperature.optimum", "PI"))
lac.short <-
  lac.short %>% separate(variable, c("var", "substrate"), extra = "merge", remove = FALSE)
lac.short$value <- as.numeric(as.character(lac.short$value))
lac.short$var <- as.factor(as.character(lac.short$var))

lac.short <- subset(lac.short, substrate != "GOH" & substrate != "LDOPA" & substrate != "PYGL" & substrate != "GUA")
lac.short$Kingdom <- ordered(lac.short$Kingdom, levels = c("Fungi", "Prokaryota", "Plantae"))
lac.short$Temperature.optimum <- as.numeric(as.character(lac.short$Temperature.optimum))


lac.fig.km <- ggplot(data = subset(lac.short, var == "Km"), aes(x = Kingdom, y = value, fill = Kingdom)) +
  geom_violin(draw_quantiles = 0.5, alpha = 0.25) +
  geom_jitter(shape = 21, width = 0.05, stroke = 0.1, size = 2, alpha = 0.75) +
  # geom_boxplot(fill = NA, outlier.alpha = 0, width = 0.25, colour = "black") +
  scale_y_log10(limits = c(0.5, 10000)) +
  labs(y = "Km [µM]") +
  scale_fill_brewer(palette = "Set1") +
  theme_base() +
  theme(text = element_text( family = "Helvetica"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.background = element_blank(),
        strip.text = element_blank(),
        plot.margin = unit(c(0,1,-1,1), "mm"),
        legend.position = "none") +
  facet_wrap(~ substrate)

lac.fig.ph <- ggplot(data = subset(lac.short, var == "pH"), aes(x = Kingdom, y = value, fill = Kingdom)) +
  geom_violin(draw_quantiles = 0.5, alpha = 0.25) +
  geom_jitter(shape = 21, width = 0.05, stroke = 0.1, size = 2, alpha = 0.75) +
  # geom_boxplot(fill = NA, outlier.alpha = 0, width = 0.25, colour = "black") +
  labs(y = "pH optimum") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 10), breaks = c(0, 3, 6, 9)) +
  theme_base() +
  theme(text = element_text( family = "Helvetica"),
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.margin = unit(c(0,1,0,1), "mm")
        ) +
  facet_wrap(~ substrate, strip.position="bottom")

lac.fig.temp <- ggplot(data = lac.short, aes(x = Kingdom, y = Temperature.optimum, fill = Kingdom)) +
  geom_violin(draw_quantiles = 0.5, alpha = 0.25) +
  geom_jitter(shape = 21, width = 0.05, stroke = 0.1, size = 2, alpha = 0.5) +
  # geom_boxplot(fill = NA, outlier.alpha = 0, width = 0.25, colour = "black") +
  labs(y = "Temperature optimum [°C]") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 95)) +
  theme_base() +
  theme(text = element_text( family = "Helvetica"),
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank()
        # plot.margin = unit(c(0,0,0,0), "cm")
  ) 

lac.fig.pi <- ggplot(data = lac.short, aes(x = Kingdom, y = PI, fill = Kingdom)) +
  geom_violin(draw_quantiles = 0.5, alpha = 0.25) +
  geom_jitter(shape = 21, width = 0.05, stroke = 0.1, size = 2, alpha = 0.5) +
  # geom_boxplot(fill = NA, outlier.alpha = 0, width = 0.25, colour = "black") +
  labs(y = "Isoelectric point") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 10.5), breaks = c(0, 3, 6, 9)) +
  theme_base() +
  theme(text = element_text( family = "Helvetica"),
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank()
        # plot.margin = unit(c(0,0,0,0), "cm")
  )

pdf("lac_fig.pdf", width = 7, height = 5)
kinetics <- plot_grid(lac.fig.km,
          lac.fig.ph,
          ncol = 1,
          align = "v",
          rel_heights = c(1, 1.4)
          )
plot_grid(kinetics,
          lac.fig.temp,
          lac.fig.pi,
          ncol = 3,
          align = "h",
          labels = c("(A)", "(B)", "(C)"),
          label_fontfamily = "Helvetica",
          label_x = 0,
          hjust = 0,
          # axis = "b",
          rel_widths = c(1, 0.5, 0.5)
)
dev.off()