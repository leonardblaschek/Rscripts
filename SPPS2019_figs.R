library(showtext)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(gghighlight)


#### import Helvetica Neue ####
font_add("Helvetica",
  regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
  italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
  bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf"
)
showtext_auto()


#### import measurements ####
phlog.monol <-
  read.csv("/home/leonard/Documents/Uni/Phloroglucinol/measurements_revisited.csv",
    skip = 2
  )

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

phlog.monol$cell.type <- factor(phlog.monol$cell.type)

#### calculate the correct hue on the 360 point circular scale ####
phlog.monol$hue <- ((phlog.monol$h.stained + 128) / 255 * 360)

phlog.monol$replicate <-
  as.factor(as.character(phlog.monol$replicate))


#### calculate stained - unstained diff. and adjust for bleaching by subtracting the diff. for the unlignified phloem ####
phlog.monol$diff <-
  phlog.monol$OD.stained - phlog.monol$OD.unstained
phlog.monol.bg <- phlog.monol %>%
  filter(cell.type == "PH") %>%
  select(1:4, 10) %>%
  group_by(genotype, replicate, technical) %>%
  summarise(OD.bg = mean(diff, na.rm = TRUE))

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
phlog.monol.pre <- phlog.monol %>%
  group_by(genotype, cell.type, replicate) %>%
  summarise(
    mean.hue1 = mean(hue, na.rm = TRUE),
    SD.hue1 = sd(hue, na.rm = TRUE),
    mean.OD1 = mean(diff.adj, na.rm = TRUE),
    SD.OD1 = sd(diff.adj, na.rm = TRUE)
  )


#### average per genotype (for barplots) ####
phlog.monol.avg <- phlog.monol.pre %>%
  group_by(genotype, cell.type) %>%
  summarise(
    mean.hue2 = mean(mean.hue1, na.rm = TRUE),
    SD.hue2 = sd(mean.hue1, na.rm = TRUE),
    mean.OD2 = mean(mean.OD1, na.rm = TRUE),
    SD.OD2 = sd(mean.OD1, na.rm = TRUE)
  )

#### plot absorbance ####
ab_plot <- function(ct, genotype1, genotype2, genotype3) {
  ggplot(filter(phlog.monol.avg, cell.type == ct)) +
    geom_col(aes(x = genotype, y = mean.OD2), fill = "#04253A") +
    gghighlight(genotype %in% c(genotype1, genotype2, genotype3)) +
    geom_jitter(
      data = filter(
        phlog.monol.pre,
        cell.type == ct &
          genotype %in% c(genotype1, genotype2, genotype3)
      ),
      aes(x = genotype, y = mean.OD1),
      colour = "#04253A",
      fill = "#FFCC3D",
      shape = 21,
      size = 2,
      width = 0.2
    ) +
    labs(y = "Absorbance") +
    scale_y_continuous(expand = c(0,0.01,0,0.01), limits = c(-0.005, 0.5)) +
    scale_x_discrete(
      labels = c(
        "Col-0",
        "",
        "",
        expression(paste(italic("4cl1"), "x", italic("4cl2"))),
        "",
        "",
        "",
        expression(italic("ccr1")),
        "",
        "",
        "",
        ""
      )
    ) +
    theme_minimal() +
    theme(text = element_text(family = "Helvetica", size = 14),
          axis.text = element_text(colour = "#04253A"),
          axis.title.y = element_text(size = 12, colour = "#04253A"),
          axis.title.x = element_blank(),
          panel.grid = element_blank())
}
pdf("IF_col_4cl.pdf", width = 3, height = 3)
ab_plot("IF", "col-0", "ccr1-3", "4cl1x2")
dev.off()
pdf("MX_col_4cl.pdf", width = 3, height = 3)
ab_plot("MX", "col-0", "ccr1-3", "4cl1x2")
dev.off()
pdf("PX_col_4cl.pdf", width = 3, height = 3)
ab_plot("PX", "col-0", "ccr1-3", "4cl1x2")
dev.off()
pdf("XF_col_4cl.pdf", width = 3, height = 3)
ab_plot("XF", "col-0", "ccr1-3", "4cl1x2")
dev.off()

