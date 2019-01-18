library(vegan)
library(agricolae)
library(sysfonts)
library(showtext)
library(ggplot2)
library(ggthemes)
library(plyr)
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


#### import measurements ####
phlog.monol <-
  read.csv("/home/leonard/Documents/Uni/Phloroglucinol/measurements_revisited.csv",
           skip = 2)

#### calculate pixel values from OD ####
phlog.monol$genotype <-
  ordered(
    phlog.monol$genotype,
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
      "cad4x5"
    )
  )


#### set cell types according to measurement order ####
phlog.monol[1:50 + rep(seq(0, (nrow(phlog.monol) - 50), by = 300), each = 50), 4] <-
  "IF"
phlog.monol[51:100 + rep(seq(0, (nrow(phlog.monol) - 50), by = 300), each = 50), 4] <-
  "MX"
phlog.monol[101:150 + rep(seq(0, (nrow(phlog.monol) - 50), by = 300), each = 50), 4] <-
  "XF"
phlog.monol[151:200 + rep(seq(0, (nrow(phlog.monol) - 50), by = 300), each = 50), 4] <-
  "PX"
phlog.monol[201:250 + rep(seq(0, (nrow(phlog.monol) - 50), by = 300), each = 50), 4] <-
  "LP"
phlog.monol[251:300 + rep(seq(0, (nrow(phlog.monol) - 50), by = 300), each = 50), 4] <-
  "PH"


#### calculate the correct hue on the 360 point circular scale ####
phlog.monol$hue <- ((phlog.monol$h.stained + 128) / 255 * 360)

phlog.monol$replicate <-
  as.factor(as.character(phlog.monol$replicate))


#### calculate stained - unstained diff. and adjust for bleaching by subtracting the diff. for the unlignified phloem ####
phlog.monol$diff <-
  phlog.monol$OD.stained - phlog.monol$OD.unstained
phlog.monol.bg <-
  ddply(
    subset(phlog.monol, cell.type == "PH", select = c(1, 2, 3, 4, 10)),
    c("genotype", "replicate", "technical"),
    summarise,
    OD.bg = mean(diff, na.rm = TRUE)
  )
phlog.monol.bg$cell.type <- NULL
phlog.monol <-
  merge(
    phlog.monol,
    phlog.monol.bg,
    all = TRUE,
    by = c("genotype", "replicate", "technical")
  )
phlog.monol$diff.adj <- phlog.monol$diff - phlog.monol$OD.bg
phlog.monol <- subset(phlog.monol, cell.type != "PH")


#### average per replicate (for boxplots) ####
phlog.monol.pre <-
  ddply(
    phlog.monol,
    c("genotype", "cell.type", "replicate"),
    summarise,
    mean.hue1 = mean(hue, na.rm = TRUE),
    SD.hue1 = sd(hue, na.rm = TRUE),
    mean.OD1 = mean(diff.adj, na.rm = TRUE),
    SD.OD1 = sd(diff.adj, na.rm = TRUE)
  )


#### average per genotype (for barplots) ####
phlog.monol.avg <-
  ddply(
    phlog.monol.pre,
    c("genotype", "cell.type"),
    summarise,
    mean.hue2 = mean(mean.hue1, na.rm = TRUE),
    SD.hue2 = sd(mean.hue1, na.rm = TRUE),
    mean.OD2 = mean(mean.OD1, na.rm = TRUE),
    SD.OD2 = sd(mean.OD1, na.rm = TRUE)
  )

#### calculate differences from the WT and range normalise ####
phlog.monol.avg.WT <- subset(phlog.monol.avg, genotype == "col-0")
colnames(phlog.monol.avg.WT)[3] <- "hue.WT"
colnames(phlog.monol.avg.WT)[5] <- "OD.WT"
phlog.monol.avg <- merge(phlog.monol.avg, phlog.monol.avg.WT[, c(2, 3, 5)], all = TRUE)
phlog.monol.avg$mean.OD2 <- phlog.monol.avg$mean.OD2 - phlog.monol.avg$OD.WT
phlog.monol.avg$mean.OD2 <- phlog.monol.avg$mean.OD2 / max(abs(phlog.monol.avg$mean.OD2))
phlog.monol.avg$mean.hue2 <- phlog.monol.avg$mean.hue2 - phlog.monol.avg$hue.WT
phlog.monol.avg$mean.hue2 <- phlog.monol.avg$mean.hue2 / max(abs(phlog.monol.avg$mean.hue2))

phlog.matrix <- subset(phlog.monol.avg, select= c(1,2,5))
phlog.matrix <- dcast(phlog.matrix, genotype ~ cell.type)
phlog.matrix.rownames <- phlog.matrix[,1]
phlog.matrix <- data.matrix(phlog.matrix[, 2:6])
rownames(phlog.matrix) <- phlog.matrix.rownames

#### calculate euclidean distances ####
d <- vegdist(t(phlog.matrix), method = "euclidean")
h <- hclust(d, method = "average")
h <- reorder(h, d, "sum")

d2 <- vegdist(phlog.matrix, method = "euclidean")
h2 <- hclust(d2, method = "average")
h2 <- reorder(h2, d2, "min")

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
    margins = c(2, 10),
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

pdf(file = "heatmap_phlog_OD.pdf", height = 4, width = 10)
htmp()
dev.off()