library(dplyr)
library(tidyr)
library(ggplot2)
sem.data <- read.csv("file:///home/leonard/R/Output/PhD/SEM_data.csv") %>%
  select(-X, -X1599...cm....1.1119..cm....1)

summary(sem.data)

# sem.gather <- sem.data %>%
#   gather(key = "variable",
#          value = "value",
#          - genotype,
#          - replicate,
#          - cell.type,
#          - technical)
# head(sem.gather)
ggplot(sem.data, aes(x = genotype, y = Integrated.lignin.peak)) +
  geom_violin() +
  geom_jitter() +
  theme_leo() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5)) +
  facet_wrap(~ cell.type, ncol = 3)

first.model <- psem(
  lm(Circ. ~ 
        n_v +
        mean.OD1 +
        Integrated.cellulose.peak +
        Lignin.peak.position +
        Integrated.lignin.peak,
      family = gaussian,
      data = drop_na(filter(sem.data, cell.type == "PMX"))),
  glm(Plant.height ~ 
       Circ. +
       Lignin.peak.position, 
     family = poisson(),
     data = drop_na(filter(sem.data, cell.type == "PMX"))),
  data = drop_na(filter(sem.data, cell.type == "PMX"))
)
summary(first.model, .progressBar = FALSE)
# plot(first.model)


ggplot(sem.data, aes(x = Integrated.lignin.peak, y = Circ., fill = cell.type)) +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  geom_point(shape = 21, alpha = 0.75, size = 3) +
  stat_smooth(method = "loess", linetype = 2, se = FALSE, aes(colour = cell.type)) +
  theme_leo()

ggplot(sem.data, aes(x = Circ., y = Plant.height, fill = cell.type)) +
  scale_fill_viridis_d() +
  geom_point(shape = 21, alpha = 0.75, size = 3) +
  stat_smooth(method = "glm", method.args = list(family = "poisson"), colour = "black", linetype = 2, se = FALSE, aes(group = 1)) +
  theme_leo()