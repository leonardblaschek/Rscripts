library(tidyr)
library(dplyr)
library(factoextra)
library(ggfortify)
library(showtext)
library(cowplot)
library(ggthemes)

#### import Helvetica Neue ####
font_add("Helvetica", regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf")
showtext_auto()

#### import poplar measurements ####
poplar <-
  read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol/poplar_foodweb.csv")
poplar <- poplar[, -16]

poplar$cell.type <- recode(poplar$cell.type,
                           "F" = "Fibre",
                           "V" = "Vessel",
                           "R" = "Ray",
                           "CB" = "Cambium"
)
poplar$adj.cell.type <- recode(poplar$adj.cell.type,
                               "F" = "Fibre",
                               "V" = "Vessel",
                               "R" = "Ray",
                               "CB" = "Cambium",
                               "PA" = "Parenchyma"
)

#### calculate distance to the cambium reference line 
# 5.9 is the number of pixels per µm
# X and Y are already measured according to scale, the ref values are given in pixels (see Fiji macro) ####
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

#### subtract background ####
poplar$diff <- poplar$OD.stained - poplar$OD.unstained
poplar.bg <-
  subset(
    poplar,
    cell.type == "Cambium" & adj.cell.type == "Cambium",
    select = c("genotype", "replicate", "technical", "diff")
  )
colnames(poplar.bg)[4] <- "OD.bg"

poplar.bg <- poplar.bg %>%
  group_by(genotype, replicate, technical) %>%
  mutate(OD.bg = mean(OD.bg, na.rm = TRUE))

poplar <-
  full_join(poplar,
            unique(poplar.bg),
            by = c("genotype", "replicate", "technical"))
poplar$diff.adj <- poplar$diff - poplar$OD.bg

#### bin measurements by distance to the cambium ####
poplar.bin <- poplar %>%
  mutate(bin = cut(
    Distance,
    breaks = c(-Inf, 50, 100, Inf),
    labels = c("I", "II", "III")
  ))

poplar.bin.pre <- poplar.bin %>%
  filter(cell.type %in% c("Vessel", "Fibre", "Ray")) %>%
  group_by(genotype, bin, cell.type, replicate) %>%
  summarise(od = mean(diff.adj))

pca.poplar <- poplar.bin.pre %>%
  unite(cellbin, cell.type, bin, sep = " ") %>%
  spread(genotype, od)
  
#### make PCA ####
log.poplar <- log(pca.poplar[, -(1:2)])
celltypes.poplar <- as.character(pca.poplar[, 1])
pca.poplar.post <- prcomp(log.poplar, center = TRUE, scale. = TRUE)

gg.pca <- data.frame(pca.poplar.post$x, cellbin = pca.poplar$cellbin)

colors = c("Fibre I" = "#7DD4E4",
           "Fibre II" = "#0490AA",
           "Fibre III" = "#025666",
           "Ray I" = "#FFF0AA",
           "Ray II" = "#AA9539",
           "Ray III" = "#554600",
           "Vessel I" = "#B97CAC",
           "Vessel II" = "#7C296A",
           "Vessel III" = "#3E0030")
shapes = c("Fibre I" = 21,
           "Fibre II" = 21,
           "Fibre III" = 21,
           "Ray I" = 22,
           "Ray II" = 22,
           "Ray III" = 22,
           "Vessel I" = 23,
           "Vessel II" = 23,
           "Vessel III" = 23)

#### plot PCA ####
pca <- ggplot(gg.pca, aes(x = PC1, y = PC2, fill = cellbin)) + 
  geom_hline(yintercept = 0, linetype = 1) +
  geom_vline(xintercept = 0, linetype = 1) +
  stat_ellipse(geom = "polygon", alpha = 0.5, type = "t", level = 0.4) +
  stat_ellipse(aes(fill = NULL), colour = "black", linetype = 2) +
  geom_point(aes(shape = cellbin),
             size = 4,
             stroke = 0.5,
             alpha = 1) +
  labs(x = "PC 1 (77.2%)",
       y = "PC 2 (12.1%)") +
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  scale_colour_viridis_d() +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica"),
    axis.ticks = element_line(size = 0.5, lineend = "square", color = "black"),
    axis.title = element_text(size = 20),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.text.x = element_text(
      size = 20,
      colour = "black",
      angle = 0,
      vjust = 0.5,
      hjust = 0.5
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.5),
    legend.position = c(0.1, 0.8),
    legend.title = element_blank(),
    legend.text = element_text(size = 15, colour = "black"),
    plot.margin = unit(c(2,2,2,2), "mm")
  ) 

pdf("PCA_poplar.pdf")
pca
dev.off()

#### plot PCA variables and loadings ####
top <- plot_grid(
  fviz_pca_var(pca.poplar.post,
               repel = TRUE,
               title = "",
               font.family = "Helvetica",
               col.var = "black",
               alpha.var = 0.4,
               ggtheme = theme_few()
  ),
  fviz_eig(pca.poplar.post, 
           main = "", 
           font.family = "Helvetica",
           ggtheme = theme_few(),
           barcolor = NA,
           barfill = "#1d91c0"
           
  ),
  labels = c("(a)", "(b)"),
  label_fontfamily = "Helvetica",
  rel_widths = c(2,1),
  nrow = 1,
  align = "h"
)

mid <- plot_grid(
  fviz_contrib(pca.poplar.post, 
               choice = "var", 
               axes = 1, 
               top = 12, 
               title = "", 
               font.family = "Helvetica",
               ggtheme = theme_few(),
               color = NA,
               fill = "#1d91c0"
  ),
  fviz_contrib(pca.poplar.post, 
               choice = "var", 
               axes = 2, 
               top = 12, 
               title = "", 
               font.family = "Helvetica",
               ggtheme = theme_few(),
               color = NA,
               fill = "#1d91c0"
  ),
  labels = c("(c)", "(d)"),
  label_fontfamily = "Helvetica"
)

pdf("supp_poplar.pdf")
plot_grid(
top,
mid,
nrow = 2,
rel_heights = c(2,1)
)
dev.off()