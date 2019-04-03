library(dplyr)
library(ggplot2)
library(agricolae)

data <- read.csv("file:///home/leonard/Documents/tukey_data.csv", stringsAsFactors = TRUE)

tukey <- function(x) {
  aov1 <- aov(data = x, value ~ genotype)
  groups <- HSD.test(aov1, "genotype", alpha = 0.05)
  groups$groups[["genotype"]] <- rownames(groups$groups)
  groups$groups[["value"]] <- 0 # if your y-axis represents "value", the 0 specifies that the letters will be printed at 0
  return(groups[["groups"]])
}

letters <- data %>%
  group_by(cell.type, variable) %>%
  do(data.frame(tukey(.)))

plot <- ggplot(data = data, aes(x = genotype, y = value)) +
  geom_jitter(width = 0.2) +
  geom_text(data = letters,
            aes(label = groups),
            angle = 0,
            hjust = 0.5,
            family = "Helvetica") +
  scale_y_continuous(expand = expand_scale(mult = c(0.1,0.1))) + # this expands all y-axes by factor x to the bottome and factor y to the top
  facet_grid(cell.type ~ variable)

pdf("tukey_test.pdf")
plot
dev.off()