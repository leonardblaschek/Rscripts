library(dplyr)
library(tidyr)
library(ggplot2)

sem.data <- read.csv("file:///home/leonard/Dropbox/SEM_data.csv") %>%
  select(-replicate, -X) %>%
  filter(genotype != "ccr1xfah1") %>%
  group_by(genotype, cell.type) %>%
  summarise_all(funs(mean(., na.rm = TRUE))) %>%
  rename(Cellulose = X1119...cm....1,
         Lignin = X1599...cm....1, 
         Ratio = X1599...cm....1.1119..cm....1
         )

summary(lm(Circ. ~ Integrated.lignin.peak + 
             Integrated.cellulose.peak +
             n_v +
             n_f,
           sem.data))

#direct effect of circularity on height
rsquared(glm(Plant.height ~ Circ., family = "poisson", sem.data))
#test for possible direct effect of lignin on height
summary(lm(Plant.height ~ Circ. + Integrated.lignin.peak, sem.data))
#test for possible direct effect of cellulose on height
summary(lm(Plant.height ~ Circ. + Integrated.cellulose.peak, sem.data))

#direct effect of neighbourhood on circularity
summary(lm(Circ. ~ n_v + n_f + n_p, sem.data))
#test for possible direct effect of neighbourhood on height
summary(lm(Plant.height ~ Circ. + n_v + n_f + n_p, sem.data))

sem.ind <- read.csv("file:///home/leonard/Dropbox/raman_IRX.csv") %>%
  unite("ID", c("genotype", "replicate"))

ggplot(data = sem.ind, aes(x = ID, y = Circ.)) + 
  geom_violin() + 
  geom_point(shape = 21, aes(fill = cell.type))

summary(lm(Circ. ~ Perim. * Area, data = sem.ind))

summary(lm(Circ. ~  n_v +
             n_f +
             Perim.,
           sem.ind))
library(piecewiseSEM)
sem.data <- as.data.frame(sem.data)
betas <- summary(lm(Circ. ~ n_f + n_v,
                    data = sem.data))$coefficients[-1,1]

sem.data$composite <- rowSums(sem.data[, c("n_f", "n_v")] * betas)

head(sem.data, n = 1)
second.model <- psem(
  lm(Circ. ~ 
       composite +
       Integrated.cellulose.peak +
       Integrated.lignin.peak,
     # family = "binomial",
     data = sem.data),
  lm(Plant.height ~ Circ., data = sem.data),
  data = sem.data
)
summary(second.model, .progressBar = FALSE)
plot(second.model)