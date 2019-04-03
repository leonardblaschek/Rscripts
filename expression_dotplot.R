library(ggplot2)
library(ggthemes)
library(showtext)
library(dplyr)
library(tidyr)

font_add(
  "Helvetica",
  regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
  italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
  bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf",
  bolditalic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-BdIt.otf"
)
showtext_auto()

#### generating plot theme ####
theme_leo <- function(base_size = 12,
                      base_family = "Helvetica"){
  theme_minimal(base_size = base_size,
                base_family = base_family) %+replace%
    theme(
      strip.text = element_text(
        hjust = 0, 
        face = "italic"
      ),
      axis.ticks = element_line(
        size = 0.25,
        lineend = "square",
        color = "black"
      ),
      axis.title = element_blank(),
      axis.text.y = element_text(colour = "black", 
                                 margin = margin(t = 0, r = 2.5, b = 0, l = 0)),
      axis.text.x = element_text(
        colour = "black",
        angle = 0,
        vjust = 0.5,
        hjust = 0.5
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", size = 0.25),
      panel.spacing = unit(1.5, "mm"),
      legend.position = "none",
      
      complete = TRUE
    )
}

expression <- read.csv("file:///home/leonard/Dropbox/Review/figures_PO-review/expression_fig/expression.csv")
expression$PO <- ordered(expression$PO, levels = c("PPO", "LAC", "PRX"))

expression.avg <- expression %>%
  group_by(PO) %>%
  mutate(total = sum(count, na.rm = TRUE)) %>%
  group_by(PO, compartment, method) %>%
  mutate(relative = count/total)

expression.count <- expression %>%
  uncount(weights = count) %>%
  group_by(PO, compartment) %>%
  arrange(method) %>%
  mutate(number = row_number())

ex_plot <- ggplot(expression.count, aes(x = compartment, y = number, fill = PO)) +
  geom_jitter(data = filter(expression.count, method == "pred"),
              aes(x = compartment, y = number, fill = PO),
              alpha = 0) +
  geom_vline(xintercept = 1,
             size = 7.5,
             color = "grey95") +
  geom_vline(xintercept = 3,
             size = 7.5,
             color = "grey95") +
  geom_vline(xintercept = 5,
             size = 7.5,
             color = "grey95") +
  geom_vline(xintercept = 7,
             size = 7.5,
             color = "grey95") +
  geom_jitter(data = filter(expression.count, method == "pred"), aes(x = compartment, y = number, fill = PO),
              width = 0.2,
              shape = 21,
              alpha = 0.5,
              size = 1,
              stroke = 0.1) +
  geom_jitter(data = filter(expression.count, method == "exp"), aes(x = compartment, y = number, fill = PO),
              width = 0.2,
              shape = 21,
              alpha = 0.9,
              size = 2,
              stroke = 0.2) +
  theme_leo() +
  theme(axis.ticks.x = element_blank(),
        axis.text.y = element_text(hjust = 1)) +
  scale_fill_few() +
  facet_wrap(~ PO, nrow = 1) +
  coord_flip()

ex_plot_2 <- ggplot() +
  geom_jitter(data = filter(expression.count, method == "pred"), aes(x = PO, y = number, fill = PO),
              width = 0.2,
              shape = 21,
              alpha = 0.5,
              size = 1,
              stroke = 0.1) +
  geom_jitter(data = filter(expression.count, method == "exp"), aes(x = PO, y = number, fill = PO),
              width = 0.2,
              shape = 21,
              alpha = 0.9,
              size = 2,
              stroke = 0.2) +
  theme_leo() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_few() +
  facet_wrap(~ compartment, ncol = 1, strip.position = "right") +
  coord_flip()

pdf("expression_plot.pdf", height = 2, width = 6)
ex_plot
dev.off()

pdf("expression_plot_2.pdf", height = 7, width = 3)
ex_plot_2
dev.off()
