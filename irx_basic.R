library(reshape2)
library(ggplot2)
library(ggthemes)
library(showtext)
library(dplyr)
library(cowplot)

font_add("Helvetica", regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf",
         bolditalic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-BdIt.otf")
showtext_auto()

irx.basic <- subset(irx.data, select= c(1,2,4,6, 7,10,11))
irx.basic <- melt(irx.basic, c("genotype", "replicate", "object", "Distance"))

irx_basic <- ggplot(data = irx.basic, aes(x = genotype, y = value, fill = Distance)) +
  geom_violin(draw_quantiles = 0.5, adjust = 1.5) +
  geom_jitter(shape = 21, width = 0.1, stroke = 0.1) +
  scale_fill_viridis_c() +
  theme_base() +
  theme(text = element_text(family = "Helvetica")) +
  facet_grid(object ~ variable, scales = "free_x") +
  coord_flip()

pdf("irx_basic.pdf")
irx_basic
dev.off()