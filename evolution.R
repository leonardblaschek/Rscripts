library(ggplot2)
library(reshape2)
library(ggthemes)
library(showtext)
library(broom)

# import Helvetica Neue
font_add("Helvetica", regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf")
showtext_auto()

PO <-
  read.csv("file:///home/leonard/Dropbox/Review/gene_families.csv")
PO$species <- ordered(PO$species, levels = PO$species)
PO$taxon <- ordered(PO$taxon, levels = unique(PO$taxon))
PO$number <- row(PO)
PO <- PO[-8,]

PO <- subset(PO, select = c(1:7))
PO <- melt(PO, c("species", "taxon", "time", "genome"))

lac.lm <- PO %>%
  group_by(variable) %>%
  filter(genome < 800) %>%
  do(lin.fit = lm(value ~ genome, data = .))

lac.lm.coef <- tidy(lac.lm, lin.fit)

ppo.plot <- ggplot(data = PO, aes(x = genome, y = value)) +
  geom_point(aes(size = time, fill = variable), colour = "black", shape = 21, stroke = 0.1, alpha = 0.75) +
  scale_size_area(max_size = 10, name = "Divergence [mya]") +
  scale_x_continuous(limits = c(0, 800)) +
  scale_fill_few(guide = "none") +
  scale_colour_few( name = "") +
  geom_smooth(aes(colour = variable), method = "lm", se = FALSE, linetype = 2) +
  guides(colour = guide_legend(order = 0),
         size = guide_legend(order = 1)) +
  labs(x = "Genome size [Mb]",
       y = "Number of paralogues") +
  theme_base() +
  theme(text = element_text(family = "Helvetica"),
        legend.position = c(0.3, 0.8),
        legend.box = "horizontal",
        plot.background = element_blank())


pdf("gene_fams.pdf", width = 6, height = 6)
ppo.plot
dev.off()