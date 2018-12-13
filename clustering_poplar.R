library(vegan)
library(agricolae)
library(sysfonts)
library(showtext)
library(ggplot2)
library(ggthemes)
library(gplots)
library(RColorBrewer)
library(reshape2)
library(Hmisc)
library(proxy)

font_add("Helvetica", regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf")
showtext_auto()

color.fun <- colorRampPalette(brewer.pal(11, "RdBu"))

poplar.cluster <- poplar.bin.avg
poplar.cluster.WT <- subset(poplar.cluster, genotype == "WT" & bin == "I")
colnames(poplar.cluster.WT)[7] <- "hue.WT"
colnames(poplar.cluster.WT)[4] <- "OD.WT"
poplar.cluster <- merge(poplar.cluster, poplar.cluster.WT[, c(3, 4, 7)], all = TRUE)
poplar.cluster$mean.od <- poplar.cluster$mean.od - poplar.cluster$OD.WT
poplar.cluster$mean.od <- poplar.cluster$mean.od / max(abs(poplar.cluster$mean.od))
poplar.cluster$mean.hue <- poplar.cluster$mean.hue - poplar.cluster$hue.WT
poplar.cluster$mean.hue <- poplar.cluster$mean.hue / max(abs(poplar.cluster$mean.hue))

# SUPERCLUSTER
# poplar.cluster.super <- melt(poplar.cluster[, c(1,2,3,5)], c("genotype", "cell.type"))
# poplar.cluster.super <- dcast(poplar.cluster.super, formula = cell.type ~ genotype + variable)
# 
# phlog.matrix.super <- data.matrix(poplar.cluster.super)
# rownames(phlog.matrix.super) <- poplar.cluster.super[, 1]
# phlog.matrix.super <- phlog.matrix.super[, -1]
# 
# d = vegdist(phlog.matrix.super, na.rm = TRUE, method = "euclidean")
# h = hclust(d, method = "average")
# pdf("supercluster.pdf")
# plot(h, hang = -1)
# dev.off()

# CLUSTER BY CELL TYPE
# Set cell type for plots
# CT <- "IF"


phlog.matrix <- subset(poplar.cluster, select = c(1,2,3,4))
phlog.matrix <- dcast(phlog.matrix, genotype + bin ~ cell.type)
phlog.matrix.rownames <- paste(phlog.matrix[,1], phlog.matrix[,2])
phlog.matrix <- data.matrix(phlog.matrix[, 3:5])
rownames(phlog.matrix) <- phlog.matrix.rownames

pdf("vegdist_absorbance_pop.pdf")
d <- vegdist(t(phlog.matrix), method = "euclidean")
h <- hclust(d, method = "average")
h <- reorder(h, d, agglo.fun = "max")

d2 <- vegdist(phlog.matrix, method = "euclidean")
h2 <- hclust(d2, method = "average")
h2 <- reorder(h2, d2, "max")

plot(h)
dev.off()

htmp <- function () {
  par(family = "Helvetica")
  heatmap.2(
    phlog.matrix,
    na.rm = TRUE,
    trace = "none",
    Colv = as.dendrogram(h),
    Rowv = as.dendrogram(h2),
    # Colv = TRUE,
    dendrogram = "both",
    key = TRUE,
    # density.info = "none",
    key.title = NA,
    key.xlab = NA,
    key.ylab = NA,
    margins = c(2, 20),
    col = color.fun(15),
    # breaks = seq(-1, 1, length.out=51),
    # symbreak = FALSE,
    # symkey = FALSE,
    cexCol = 1.5,
    cexRow = 1.5,
    srtCol = 0,
    adjCol = c(0.5, 0.5),
    adjRow = c(0, 0.5),
    offsetRow = 0,
    # lmat=rbind(c(5,4), c(3,2), c(0,1)),
    # lhei=c(2,4,0.2),
    rowsep = c(1:11),
    colsep = c(1:4),
    sepwidth = c(0.01,0.02),
    denscol = "black"
  )}

pdf(file = "heatmap_pop_OD.pdf", height = 4, width = 10)
htmp()
dev.off()


# plot grid
# pdf("heatmap_grid.pdf", height = 5, width = 7.5)
# plot_grid(
#           htmp,
#           labels = c('', 'L'), 
#           ncol=2, 
#           nrow = 1, 
#           label_fontfamily = "Helvetica",
#           rel_widths = c(1, 1),
#           scale = c(1,0.9))
# dev.off()
# system("gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/default -dNOPAUSE -dQUIET -dBATCH -dDetectDuplicateImages -dCompressFonts=true -r150 -sOutputFile=heat_grid.pdf heatmap_grid.pdf")
# 
