library(showtext)
library(ggplot2)
library(ggthemes)
library(dplyr)
# library(ggthemr)


#### import Helvetica Neue ####
font_add("Helvetica",
         regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf")
showtext_auto()

#### generating plot theme ####
theme_leo <- function(base_size = 12,
                      base_family = "Helvetica"){
  theme_minimal(base_size = base_size,
                base_family = base_family) %+replace%
    theme(
      # axis.line = element_line(
      #   size = 0.25,
      #   lineend = "square",
      #   color = "black"
      # ),
      strip.text = element_text(hjust = 0, face = "italic"),
      axis.ticks = element_line(
        size = 0.25,
        lineend = "square",
        color = "black"
      ),
      axis.title = element_blank(),
      axis.text.x = element_text(colour = "black", # flipped coords
                                 margin = margin(r = 1)),
      axis.text.y = element_text(
        colour = "black",
        angle = 0,
        vjust = 0.5,
        hjust = 1,
        margin = margin(r = 1)
      ),
      # axis.text.y = element_text(colour = "black",
      #                            margin = margin(r = 1)),
      # axis.text.x = element_text(
      #   colour = "black",
      #   angle = 90,
      #   vjust = 0.5,
      #   hjust = 1,
      #   margin = margin(t = 1)
      # ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", size = 0.25),
      panel.spacing = unit(1.5, "mm"),
      legend.position = "bottom",
      legend.text = element_text(size = rel(0.8)),
      legend.key.height = unit(4, "mm"),
      legend.key.width = unit(30, "mm"),
      # plot.margin = unit(c(0, 0, 0, 0), "cm"),   
      
      complete = TRUE
    )
}

#### import measurements ####
phlog.monol <-
  read.csv("/home/leonard/Documents/Uni/Phloroglucinol/measurements_revisited.csv",
           skip = 2)

phlog.monol$genotype <- recode(phlog.monol$genotype, cad4x5 = "cad4xcad5")

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

#### import SMX measurements ####
phlog.monol.SMX <- read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol/measurements_SMX.csv",
                            skip = 2)
phlog.monol.SMX$replicate <- as.factor(phlog.monol.SMX$replicate)

phlog.monol <- full_join(select(phlog.monol, -technical), select(phlog.monol.SMX, -technical))

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
      "cad4xcad5"
    )
  )

#### calculate the correct hue on the 360 point circular scale ####
phlog.monol$hue <- ((phlog.monol$h.stained + 128) / 255 * 360)

phlog.monol$replicate <-
  as.factor(as.character(phlog.monol$replicate))

#### calculate stained - unstained diff. and adjust for bleaching by subtracting the diff. for the unlignified phloem ####
phlog.monol$diff <-
  phlog.monol$OD.stained - phlog.monol$OD.unstained
phlog.monol.bg <- phlog.monol %>%
  filter(cell.type == "PH") %>%
  select(1:3, 9) %>%
  group_by(genotype, replicate) %>%
  summarise(OD.bg = mean(diff, na.rm = TRUE))

phlog.monol.bg$cell.type <- NULL
phlog.monol <-
  merge(
    phlog.monol,
    phlog.monol.bg,
    all = TRUE,
    by = c("genotype", "replicate")
  )
phlog.monol$diff.adj <- phlog.monol$diff - phlog.monol$OD.bg
phlog.monol <- subset(phlog.monol, cell.type != "PH")

#### average per replicate (for boxplots) ####
phlog.monol.pre <-  phlog.monol %>%
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
# ggthemr("dust")
wiesner_smx <- ggplot(phlog.monol.avg, aes(x = cell.type, y = mean.OD2)) +
  theme_leo() +
  scale_y_continuous(limits = c(-0.01,0.55), expand = c(0.0,0.0)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.85),
           width = 0.75) +
  facet_wrap(~ genotype, nrow = 6) +
  geom_jitter(
    data = phlog.monol.pre,
    aes(x = cell.type, y = mean.OD1),
    shape = 21,
    width = 0.1,
    stroke = 0.1,
    alpha = 0.95,
    fill = "white"
  ) +
  # scale_x_discrete(
  #   labels = c(
  #     "Col-0",
  #     expression(italic("4cl1")),
  #     expression(italic("4cl2")),
  #     expression(paste(italic("4cl1"), "x", italic("4cl2"))),
  #     expression(italic("ccoaomt1")),
  #     expression(italic("fah1")),
  #     expression(italic("omt1")),
  #     expression(italic("ccr1")),
  #     expression(paste(italic("ccr1"), "x", italic("fah1"))),
  #     expression(italic("cad4")),
  #     expression(italic("cad5")),
  #     expression(paste(italic("cad4"), "x", italic("cad5")))
  #   )) +
  coord_flip()

pdf("wiesner_smx.pdf", height = 9, width = 3)
wiesner_smx
dev.off()
