library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(showtext)

font_add("Helvetica",
         regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf")
showtext_auto()

comp <-
  read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol/chiang_comp.csv")
sacc <-
  read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol/chiang_sacc.csv")
wood <-
  read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol/chiang_wood.csv")
moe <-
  read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol/chiang_moe.csv")

comp.avg <- ddply(comp, c("Line.ID"), summarise,
                  mean.aldehydes = mean(Aldehydes),
                  mean.sg = mean(S.G.Ratio))

sacc.avg <- ddply(
  sacc,
  c("Line.ID"),
  summarise,
  mean.glu = mean(Glu.Unpret),
  mean.glu.pret = mean(Glu.Pret),
  mean.xyl = mean(Xyl.Unpret),
  mean.xyl.pret = mean(Xyl.Pret)
)

wood.avg <- ddply(wood, c("Line.ID", "Target.Gene.s."), summarise,
                  mean.lignin = mean(Lignin.))

chiang.data <- merge(wood.avg, sacc.avg, all = TRUE)
chiang.data <- merge(chiang.data, comp.avg, all = TRUE)
chiang.data <- merge(chiang.data, moe, all = TRUE)

r2.xyl <-
  round(summary(lm((mean.aldehydes / mean.lignin) ~ mean.xyl.pret,
                   data = subset(chiang.data, mean.aldehydes > 0)
  ))[["adj.r.squared"]], 2)
r2.moe <-
  round(summary(lm((mean.aldehydes / mean.lignin) ~ MOE,
             data = subset(chiang.data, mean.aldehydes > 0)
  ))[["adj.r.squared"]], 2)
r2.sg <-
  round(summary(lm(mean.sg ~ mean.glu.pret,
                   data = subset(chiang.data, mean.aldehydes > 0)
  ))[["adj.r.squared"]], 2)


chiang.scatter <-
  ggplot(data = chiang.data, aes(x = (mean.aldehydes / mean.lignin) * 100, y = mean.xyl.pret)) +
  xlim(0.01, 2.5) +
  annotate(
    "text",
    family = "Helvetica",
    label = expression(R ^ 2 == 0.81),
    x = 2.2,
    y = 0.175,
    size = 8,
    colour = "#04253a"
  ) +
  theme_few() +
  theme(
    text = element_text(family = "Helvetica", colour = "#04253a"),
    axis.ticks = element_line(colour = "#04253a"),
    panel.border = element_rect(colour = "#04253a"),
    axis.text = element_text(size = 16, colour = "#04253a"),
    axis.title = element_text(size = 20)
  ) +
  labs(x = "Lignin Aldehydes [%]",
       y = "Xylose Release") +
  geom_point(
    size = 6,
    shape = 21,
    alpha = 0.75,
    fill = "#04253a",
    colour = "#04253a"
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linetype = 2,
    colour = "#cc2236ff",
    size = 2
  )

pdf("aldehydes_glu.pdf")
chiang.scatter
dev.off()

chiang.scatter <-
  ggplot(data = chiang.data, aes(x = (mean.aldehydes / mean.lignin) * 100, y = MOE)) +
  xlim(0.01, 2.5) +
  annotate(
    "text",
    family = "Helvetica",
    label = expression(R ^ 2 == 0.67),
    x = 2.2,
    y = 4400,
    size = 8,
    colour = "#04253a"
  ) +
  theme_few() +
  theme(
    text = element_text(family = "Helvetica", colour = "#04253a"),
    axis.ticks = element_line(colour = "#04253a"),
    panel.border = element_rect(colour = "#04253a"),
    axis.text = element_text(size = 16, colour = "#04253a"),
    axis.title = element_text(size = 20)
  ) +
  labs(x = "Lignin Aldehydes [%]",
       y = "Modulus of Elasticity") +
  geom_point(
    size = 6,
    shape = 21,
    alpha = 0.75,
    fill = "#04253a",
    colour = "#04253a"
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linetype = 2,
    colour = "#cc2236ff",
    size = 2
  )

pdf("aldehydes_MOE.pdf")
chiang.scatter
dev.off()
