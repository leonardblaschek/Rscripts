library(dplyr)
library(PerformanceAnalytics)

phlog.irx <-
  read.csv("/home/leonard/Documents/Uni/Phloroglucinol/measurements_revisited.csv",
           skip = 2)

phlog.irx$genotype <- recode(phlog.irx$genotype, cad4x5 = "cad4xcad5")

#### set cell types according to measurement order ####
phlog.irx[1:50 + rep(seq(0, (nrow(phlog.irx) - 50), by = 300), each = 50), 4] <-
  "IF"
phlog.irx[51:100 + rep(seq(0, (nrow(phlog.irx) - 50), by = 300), each = 50), 4] <-
  "MX"
phlog.irx[101:150 + rep(seq(0, (nrow(phlog.irx) - 50), by = 300), each = 50), 4] <-
  "XF"
phlog.irx[151:200 + rep(seq(0, (nrow(phlog.irx) - 50), by = 300), each = 50), 4] <-
  "PX"
phlog.irx[201:250 + rep(seq(0, (nrow(phlog.irx) - 50), by = 300), each = 50), 4] <-
  "LP"
phlog.irx[251:300 + rep(seq(0, (nrow(phlog.irx) - 50), by = 300), each = 50), 4] <-
  "PH"

#### import SMX measurements ####
phlog.irx.SMX <- read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol/measurements_SMX.csv",
                            skip = 2)
phlog.irx.SMX$replicate <- as.factor(phlog.irx.SMX$replicate)

phlog.irx <- full_join(select(phlog.irx, -technical), select(phlog.irx.SMX, -technical))

phlog.irx$genotype <-
  ordered(
    phlog.irx$genotype,
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
phlog.irx$hue <- ((phlog.irx$h.stained + 128) / 255 * 360)

phlog.irx$replicate <-
  as.factor(as.character(phlog.irx$replicate))

#### calculate stained - unstained diff. and adjust for bleaching by subtracting the diff. for the unlignified phloem ####
phlog.irx$diff <-
  phlog.irx$OD.stained - phlog.irx$OD.unstained
phlog.irx.bg <- phlog.irx %>%
  filter(cell.type == "PH") %>%
  select(1:3, 9) %>%
  group_by(genotype, replicate) %>%
  summarise(OD.bg = mean(diff, na.rm = TRUE))

phlog.irx.bg$cell.type <- NULL
phlog.irx <-
  merge(
    phlog.irx,
    phlog.irx.bg,
    all = TRUE,
    by = c("genotype", "replicate")
  )
phlog.irx$diff.adj <- phlog.irx$diff - phlog.irx$OD.bg
phlog.irx <- subset(phlog.irx, cell.type != "PH")

#### average per replicate (for boxplots) ####
phlog.irx.pre <-  phlog.irx %>%
  mutate(cell.type = recode(cell.type, "MX" = "PMX"),
         genotype = recode(genotype, "col-0" = "Col-0",
                           "4cl1x2" = "4cl1x4cl2")) %>%
  group_by(genotype, cell.type, replicate) %>%
  summarise(
    mean.hue1 = mean(hue, na.rm = TRUE),
    SD.hue1 = sd(hue, na.rm = TRUE),
    mean.OD1 = mean(diff.adj, na.rm = TRUE),
    SD.OD1 = sd(diff.adj, na.rm = TRUE)
  ) %>%
  select(genotype, cell.type, replicate, mean.OD1)
  
#### AVERAGES ####
raman.spread <- raman.data.plot %>%
  ungroup() %>%
  mutate(cell.type = dplyr::recode(cell.type, "MX" = "PMX")) %>%
  group_by(genotype, cell.type, replicate) %>%
  select(-value.scaled) %>%
  spread(variable, value)


raman.irx <- read.csv("file:///home/leonard/Dropbox/raman_IRX.csv") %>%
  select(genotype:Perim., Circ., Round)
raman.irx$replicate <- as.character(raman.irx$replicate)
raman.irx$technical <- as.character(raman.irx$technical)

raman.heights <- read.csv("file:///home/leonard/Documents/Uni/Master/Summer project 16/phenotyping/phenotyping.csv") %>%
  select(genotype, replicate, height) %>%
  mutate(genotype = recode(genotype, "col-0" = "Col-0")) %>%
  rename("Plant.height" = "height") %>%
  rbind(read.csv("file:///home/leonard/Documents/Uni/PhD/IRX/heights_haris.csv"))

raman.heights$replicate <- as.character(raman.heights$replicate)

#### linear regression with lignin content ####
raman.irx <- left_join(raman.irx, raman.spread)
raman.irx <- left_join(raman.irx, raman.heights)
raman.irx <- left_join(raman.irx, phlog.irx.pre)
write.csv(raman.irx, file = "SEM_data.csv")

lignin_plot <- ggplot(data = irx.corr.lm, aes(x = Circ., y = intensity1599)) +
  geom_point(shape = 21, stroke = 0.2, size = 2, aes(fill = cell.type)) +
  geom_smooth(method = "lm", aes(colour = cell.type), se = FALSE)

#### basic multiple linear regression ####
irx.corr.lm <- inner_join(irx.spread, raman.spread) %>%
  ungroup() %>%
  filter(cell.type == "PMX") %>%
  select(-genotype, -cell.type)

model.null <- lm(Circ. ~ 1, data = irx.corr.lm)
model.full <- lm(Circ. ~ Perim. + Distance + intensity1119 + intensity1599, data = irx.corr.lm)

lm.best <- step(
  model.null,
  scope =
    list(upper = model.full),
  direction = "both",
  data = irx.corr.lm
)

irx.corr.lm$pred <- predict(lm.best)

lm_plot <- ggplot(data = irx.corr.lm, aes(x = Circ., y = pred)) +
  ggtitle(summary(lm.best)$call) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme(plot.title = element_text(family = "Helvetica",size = 10))

#### mixed glm model ####
irx.corr.glm <- inner_join(irx.spread, raman.spread) %>%
  ungroup() %>%
  select(-genotype)

glm.null <- glm(Circ. ~ 1, data = irx.corr.glm)
glm.full <- glm(Circ. ~ Area + Perim. + Distance + intensity1119 + intensity381 + intensity1599 + intensity1621 + intensity1662 + cell.type, data = irx.corr.glm)

step(
  glm.null,
  scope =
    list(upper = glm.full),
  direction = "both",
  data = irx.corr.glm
)

glm.best <- glm(formula = Circ. ~ cell.type + Area + Perim. + intensity381 + 
                  Distance + LigBy1119, data = irx.corr.glm)


glm.best.noratio <- glm(formula = Circ. ~ cell.type + Area + Perim. + intensity381 + 
                          Distance + intensity1119 + intensity1662, data = irx.corr.glm)

irx.corr.glm$pred <- predict(glm.best)
irx.corr.glm$pred_noratio <- predict(glm.best.noratio)

glm_plot <- ggplot(data = irx.corr.glm, aes(x = Circ., y = pred)) +
  ggtitle(summary(glm.best)$call) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme(plot.title = element_text(family = "Helvetica",size = 10))

glm_plot_noratio <- ggplot(data = irx.corr.glm, aes(x = Circ., y = pred_noratio)) +
  ggtitle(summary(glm.best.noratio)$call) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme(plot.title = element_text(family = "Helvetica",size = 10))

pdf("models.pdf")
chart.Correlation(irx.corr.lm)
lignin_plot
lm_plot
glm_plot
glm_plot_noratio
dev.off()

#### RAW DATA ####
raman.spread <- raman.data.corrected %>%
  ungroup() %>%
  mutate(cell.type = dplyr::recode(cell.type, "MX" = "PMX")) %>%
  group_by(genotype, cell.type, wavenumber, replicate, technical) %>%
  filter(wavenumber %in% c(381, 1119, 1599, 1621, 1662)) %>%
  select(-intensity) %>%
  rename("intensity" = "corrected.intensity") %>%
  gather(variable, value, -(genotype:wavenumber)) %>%
  unite(WN_var, variable, wavenumber, sep = "") %>%
  spread(WN_var, value) %>%
  mutate(LigBy1119 = intensity1599/intensity1119,
         LigBy381 = intensity1599/intensity381) %>%
  ungroup() %>%
  select(-replicate, -technical) 

irx.spread <- irx.data %>%
  rename(cell.type = object) %>%
  ungroup() %>%
  select(-technical, -replicate, -number)


irx.corr.lm <- inner_join(raman.spread, irx.spread, by= c("genotype", "cell.type")) %>%
  ungroup() %>%
  select(-genotype)

model.null <- glm(Circ. ~ 1, data = irx.corr.lm)
model.full <- glm(Circ. ~ Perim. + Distance + intensity1119 + intensity1599 + cell.type, data = irx.corr.lm)

lm.best <- step(
  model.null,
  scope =
    list(upper = model.full),
  direction = "both",
  data = irx.corr.lm
)

irx.corr.lm$pred <- predict(lm.best)

lm_plot <- ggplot(data = irx.corr.lm, aes(x = Circ., y = pred)) +
  ggtitle(summary(lm.best)$call) +
  geom_point(shape = 21, aes(fill = cell.type)) +
  geom_smooth(method = "lm", aes(colour = cell.type)) +
  theme(plot.title = element_text(family = "Helvetica",size = 10))

pdf("raw_model.pdf")
lm_plot
dev.off()
