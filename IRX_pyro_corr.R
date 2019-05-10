library(dplyr)
library(broom)
library(reshape2)
library(ggplot2)
py <-
  read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol/pyrolysis_2018.csv")

py <-
  subset(py, genotype != "pal1" &
           genotype != "pal2" & genotype != "WT_cad")
py$Coniferaldehyde <- py$X22
py$Vanillin <- py$X18
py$Sinapaldehyde <- py$X32
py$Syringaldehyde <- py$X29
py$Aldehydes <-
  ifelse(
    py$variable == "mean",
    py$Coniferaldehyde + py$Sinapaldehyde + py$Vanillin + py$Syringaldehyde,
    sqrt(
      py$Coniferaldehyde ^ 2 + py$Sinapaldehyde ^ 2 + py$Vanillin ^ 2 + py$Syringaldehyde ^ 2
    )
  )
py$Lignin <-
  ifelse(py$variable == "mean", rowSums(py[, 3:28]), sqrt(rowSums(py[, 3:28] ^
                                                                    2)))

py.melt <-
  melt(
    subset(py, select = c(1, 2, 29:34)),
    id = c("genotype", "variable"),
    variable.name = "residue"
  )
py.melt$sd <- py.melt$value
py.melt <-
  merge(
    subset(py.melt, variable == "mean", select = c(1, 3, 4)),
    subset(py.melt, variable == "sd", select = c(1, 3, 5))
  )

py.melt$genotype <- recode(py.melt$genotype, "col-0" = "Col-0", "cad4x5" = "cad4xcad5", "4cl1x2" = "4cl1x4cl2")

irx.melt.avg <- irx.melt %>%
  group_by(genotype, object, variable) %>%
  dplyr::summarise(mean.value = mean(value))

irx.pyr.avg <- irx.melt.avg %>%
  full_join(., py.melt[, c(1:3)], by = "genotype")

irx.pyr.corr <- irx.pyr.avg %>%
  filter(residue == "Lignin", variable == "Circ.") %>%
  ungroup() %>%
  select(object, mean.value, value) %>%
  group_by(object) %>%
  do(lignin.fit = lm(mean.value ~ value, data = .))

irx.py.corr.coef <- tidy(irx.pyr.corr, lignin.fit)

irx.pyr.avg <- irx.pyr.avg %>%
  full_join(., subset(irx.py.corr.coef, term == "value"))

plot_df <- irx.pyr.avg %>%
  filter(residue == "Lignin", variable == "Circ.") %>%
  ungroup() %>%
  select(object, mean.value, value, p.value) %>%
  group_by(object) %>%
  do(plots = ggplot(data = ., aes(x = value, y = mean.value)) +
       geom_point() + 
       geom_smooth(method = "lm") + 
       labs(y = "Lignin content",
            x = "Circularity") +
       ggtitle(.$object)) 

pdf("irx_corr.pdf")
plot_df$plots
dev.off()
