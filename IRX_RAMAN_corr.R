library(car)
library(dplyr)
library(PerformanceAnalytics)
library(lme4)


raman.spread <- raman.data.corrected %>%
  ungroup() %>%
  mutate(cell.type = dplyr::recode(cell.type, "MX" = "PMX")) %>%
  group_by(genotype, cell.type, wavenumber) %>%
  filter(wavenumber %in% c(381, 1119, 1599, 1621, 1662)) %>%
  summarise(
    intensity = mean(corrected.intensity, na.rm = TRUE)
  ) %>%
  gather(variable, value, -(genotype:wavenumber)) %>%
  unite(WN_var, variable, wavenumber, sep = "") %>%
  spread(WN_var, value) %>%
  mutate(LigBy1119 = intensity1599/intensity1119,
         LigBy381 = intensity1599/intensity381)

irx.spread <- irx.data %>%
  rename(cell.type = object) %>%
  group_by(genotype, cell.type) %>%
  select(-technical, -replicate, -number) %>%
  summarise_all(mean)


#### basic multiple linear regression ####
irx.corr.lm <- inner_join(irx.spread, raman.spread) %>%
  ungroup() %>%
  select(-genotype)

# irx.corr.lm$Circ. <- yjPower(irx.corr.glm$Circ., 1.5)

model.null <- lm(Circ. ~ 1, data = irx.corr.lm)
model.full <- lm(Circ. ~ cell.type + Area + Perim. + Distance + intensity1119 + intensity381, data = irx.corr.lm)

step(
  model.null,
  scope =
    list(upper = model.full),
  direction = "both",
  data = irx.corr.lm
)

lm.best <- lm(formula = Circ. ~ cell.type + Area + Perim. + intensity381 + 
                Distance + intensity1119, data = irx.corr.lm)

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

# irx.corr.glm$Circ. <- yjPower(irx.corr.glm$Circ., 1.5)

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

#### mixed lme4 model ####
irx.corr.lmer <- inner_join(irx.spread, raman.spread) %>%
  ungroup() %>%
  select( -genotype)

lmer.best <- lmer(Circ. ~ cell.type + Area + Perim. + intensity381 + 
                    Distance + LigBy1119 + (1 | cell.type), data = irx.corr.lmer)

irx.corr.lmer$pred <- predict(lmer.best)

lme4_plot <- ggplot(data = irx.corr.lmer, aes(x = Circ., y = pred)) +
  labs(title = summary(lmer.best)$call) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme(plot.title = element_text(family = "Helvetica",size = 10))

pdf("models.pdf")
chart.Correlation(select(irx.corr.lm, -cell.type))
lm_plot
glm_plot
glm_plot_noratio
lme4_plot
dev.off()
