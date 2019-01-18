library(vegan)
library(agricolae)
library(sysfonts)
library(showtext)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(reshape2)
library(Hmisc)
library(proxy)

#### import Helvetica Neue ####
font_add("Helvetica", regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf")
showtext_auto()

color.fun <- colorRampPalette(brewer.pal(11, "RdBu"))

#### import poplar measurements ####
poplar <-
  read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol/poplar_foodweb.csv")
poplar <- poplar[, -16]
poplar$replicate <- factor(poplar$replicate)
poplar$technical <- factor(poplar$technical)
poplar$number <- row(poplar)
poplar$cell.type <- plyr::revalue(poplar$cell.type,
                                  c(
                                    "F" = "Fibre",
                                    "V" = "Vessel",
                                    "R" = "Ray",
                                    "CB" = "Cambium"
                                  ))
poplar$adj.cell.type <- plyr::revalue(
  poplar$adj.cell.type ,
  c(
    "F" = "Fibre",
    "V" = "Vessel",
    "R" = "Ray",
    "CB" = "Cambium",
    "PA" = "Parenchyma"
  )
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
  merge(poplar,
        unique(poplar.bg),
        by = c("genotype", "replicate", "technical"))

poplar$diff.adj <- poplar$diff - poplar$OD.bg
poplar$genotype <-
  ordered(poplar$genotype, levels = c("WT", "c4h", "ccr"))


#### calculate the correct hue on the 360 point circular scale ####
poplar$hue <- ((poplar$H.stained + 128) / 255 * 360)

poplar <- poplar %>%
  group_by(genotype, replicate, technical) %>%
  mutate(rel.dist = Distance / max(Distance, na.rm = TRUE))

poplar <- poplar %>%
  group_by(genotype, replicate, technical, cell.type, adj.cell.type) %>%
  mutate(rel.od = diff.adj / max(diff.adj, na.rm = TRUE))

poplar <-
  subset(
    poplar,
    cell.type != "Cambium" &
      adj.cell.type != "Parenchyma" & adj.cell.type != "Cambium"
  )
poplar$cell.type <- as.character(poplar$cell.type)
poplar$adj.cell.type <- as.character(poplar$adj.cell.type)
# poplar <- subset(poplar, cell.type == adj.cell.type) # uncomment to select only self anjacent cell walls
poplar$cell.type <- as.factor(as.character(poplar$cell.type))

#### bin the measurements by distance from the cambium ####
poplar.bin <- poplar %>%
  mutate(bin = cut(
    Distance,
    breaks = c(-Inf, 50, 100, Inf),
    labels = c("I", "II", "III")
  ))

poplar.bin.pre <- poplar.bin %>%
  group_by(genotype, bin, cell.type, replicate) %>%
  summarise(
    mean.od = mean(diff.adj, na.rm = TRUE),
    sd.od = sd(diff.adj, na.rm = TRUE),
    mean.hue = mean(hue, na.rm = TRUE),
    sd.hue = sd(hue, na.rm = TRUE)
  )

poplar.bin.avg <- poplar.bin.pre %>%
  group_by(genotype, bin, cell.type) %>%
  summarise(
    mean.od = mean(mean.od, na.rm = TRUE),
    sd.od = sd(sd.od, na.rm = TRUE),
    mean.hue = mean(mean.hue, na.rm = TRUE),
    sd.hue = sd(sd.hue, na.rm = TRUE)
  )

poplar.bin$bin <-
  ordered(poplar.bin$bin, levels = c("I", "II", "III"))

#### calculate differences from the WT and range normalise ####
poplar.cluster <- poplar.bin.avg
poplar.cluster.WT <- subset(poplar.cluster, genotype == "WT" & bin == "I")
colnames(poplar.cluster.WT)[7] <- "hue.WT"
colnames(poplar.cluster.WT)[4] <- "OD.WT"
poplar.cluster <- merge(poplar.cluster, poplar.cluster.WT[, c(3, 4, 7)], all = TRUE)
poplar.cluster$mean.od <- poplar.cluster$mean.od - poplar.cluster$OD.WT
poplar.cluster$mean.od <- poplar.cluster$mean.od / max(abs(poplar.cluster$mean.od))
poplar.cluster$mean.hue <- poplar.cluster$mean.hue - poplar.cluster$hue.WT
poplar.cluster$mean.hue <- poplar.cluster$mean.hue / max(abs(poplar.cluster$mean.hue))

phlog.matrix <- subset(poplar.cluster, select = c(1,2,3,4))
phlog.matrix <- dcast(phlog.matrix, genotype + bin ~ cell.type)
phlog.matrix.rownames <- paste(phlog.matrix[,1], phlog.matrix[,2])
phlog.matrix <- data.matrix(phlog.matrix[, 3:5])
rownames(phlog.matrix) <- phlog.matrix.rownames

#### calculate euclidean distances ####
d <- vegdist(t(phlog.matrix), method = "euclidean")
h <- hclust(d, method = "average")
h <- reorder(h, d, agglo.fun = "max")

d2 <- vegdist(phlog.matrix, method = "euclidean")
h2 <- hclust(d2, method = "average")
h2 <- reorder(h2, d2, "max")

#### plot heatmap ####
htmp <- function () {
  par(family = "Helvetica")
  heatmap.2(
    phlog.matrix,
    na.rm = TRUE,
    trace = "none",
    Colv = as.dendrogram(h),
    Rowv = as.dendrogram(h2),
    dendrogram = "both",
    key = TRUE,
    key.title = NA,
    key.xlab = NA,
    key.ylab = NA,
    margins = c(2, 20),
    col = color.fun(15),
    cexCol = 1.5,
    cexRow = 1.5,
    srtCol = 0,
    adjCol = c(0.5, 0.5),
    adjRow = c(0, 0.5),
    offsetRow = 0,
    rowsep = c(1:11),
    colsep = c(1:4),
    sepwidth = c(0.01,0.02),
    denscol = "black"
  )}

pdf(file = "heatmap_pop_OD.pdf", height = 4, width = 10)
htmp()
dev.off()