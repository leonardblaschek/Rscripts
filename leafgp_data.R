library(ggplot2)
library(ggthemes)
library(showtext)
library(dplyr)
library(tidyr)
library(agricolae)

# import Helvetica Neue
font_add("Helvetica",
         regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf")
showtext_auto()

tukey <- function(x) {
  aov1 <- aov(data = x, value ~ genotype)
  groups <- HSD.test(aov1, "genotype", alpha = 0.05)
  groups$groups[["genotype"]] <- rownames(groups$groups)
  groups$groups[["value"]] <- 0
  return(groups[["groups"]])
}

leafgp.data <-
  read.csv(
    "file:///home/leonard/Documents/Uni/PhD/Phenotyping/18-09_lac_mutants_Processed_1-10-2018_comb-raw_Tray_1/2018-10-01_LeafMeasure.csv"
  )
height.data <-
  read.csv(
    "file:///home/leonard/Documents/Uni/PhD/Phenotyping/18-09_lac_mutants/originals/front/height_measurements.csv",
    skip = 1
  )

final.pheno <- read.csv("file:///home/leonard/Documents/Uni/PhD/Phenotyping/18-09_lac_mutants/18-11-1_phenotyping_lac.csv")
final.pheno <- final.pheno %>%
  spread(., measurement, value) %>%
  group_by(genotype, replicate) %>%
  mutate(maturity.ratio = siliques.mature / (siliques.mature + siliques.green + flower.count)) %>%
  gather(., key = "measurement", value = "value", -genotype, -replicate)

final.pheno$genotype <- ordered(
  final.pheno$genotype,
  levels = rev(c("Col-0", "lac4/10/12/17", "lac5/10/12/17", "lac4/5/10/12/17")))

pheno.letters <- final.pheno %>%
  group_by(measurement) %>%
  do(data.frame(tukey(.)))
  

leafgp.data <- subset(leafgp.data, Projected_LeafArea.mm.2. > 0)
leafgp.data$genotype <-
  ifelse(
    leafgp.data$Pot_ID < 6,
    "lac5/10/12/17",
    ifelse(
      leafgp.data$Pot_ID > 5 & leafgp.data$Pot_ID < 11,
      "lac4/10/12/17",
      ifelse(
        leafgp.data$Pot_ID > 10 &
          leafgp.data$Pot_ID < 16,
        "lac4/5/10/12/17",
        "Col-0"
      )
    )
  )
A <- c(1, 6, 11, 16)
B <- c(2, 7, 12, 17)
C <- c(3, 8, 13, 18)
D <- c(4, 9, 14, 19)
E <- c(5, 10, 15, 20)
leafgp.data$replicate <- as.factor(ifelse(
  leafgp.data$Pot_ID %in% A,
  "A",
  ifelse(
    leafgp.data$Pot_ID %in% B,
    "B",
    ifelse(
      leafgp.data$Pot_ID %in% C,
      "C",
      ifelse(leafgp.data$Pot_ID %in% D,
             "D",
             "E")
    )
  )
))

leafgp.data$EXP_Date <- as.character(leafgp.data$EXP_Date)
leafgp.data$EXP_Date <-
  ifelse(leafgp.data$EXP_Date == "2018-09-07", NA, leafgp.data$EXP_Date)
leafgp.data$EXP_Date <-
  ifelse(leafgp.data$EXP_Date == "2018-09-05",
         "2018-09-07",
         leafgp.data$EXP_Date)
leafgp.data$EXP_Date <-
  ifelse(is.na(leafgp.data$EXP_Date),
         "2018-09-05",
         leafgp.data$EXP_Date)
leafgp.data$dag <-
  difftime(
    as.Date(leafgp.data$EXP_Date, format = "%Y-%m-%d"),
    as.Date("2018-08-29", format = "%Y-%m-%d"),
    units = "days"
  )
leafgp.data$genotype <- ordered(
  leafgp.data$genotype,
  levels = c("Col-0", "lac4/10/12/17", "lac5/10/12/17", "lac4/5/10/12/17")
)

height.data$genotype <- ordered(
  height.data$genotype,
  levels = c("Col-0", "lac4/10/12/17", "lac5/10/12/17", "lac4/5/10/12/17")
)
height.data$dag <-
  difftime(
    as.Date(height.data$date, format = "%Y-%m-%d"),
    as.Date("2018-08-29", format = "%Y-%m-%d"),
    units = "days"
  )
height.data$height <-
  ifelse(!is.na(height.data$pot.width),
         (height.data$height * (height.data$pot.width / 725)) / 112,
         height.data$height)

full.pheno <- height.data %>%
  full_join(
    .,
    leafgp.data %>% select(dag, EXP_Date, Projected_LeafArea.mm.2., genotype, replicate),
    by = c("dag", "genotype", "replicate", "date" = "EXP_Date")
  )

leafgp.area <-
  ggplot(data = leafgp.data,
         aes(x = dag, y = Projected_LeafArea.mm.2. / 100)) +
  geom_smooth(
    aes(fill = genotype, colour = genotype),
    method = "loess",
    formula =  y ~ x,
    span = 0.75,
    se = TRUE,
    size = 0.5,
    linetype = 2
  ) +
  geom_jitter(
    aes(group = interaction(dag, genotype), fill = genotype),
    shape = 21,
    width = 0.25,
    size = 2,
    stroke = 0.25
  ) +
  geom_vline(
    xintercept = difftime("2018-09-22", "2018-08-29", units = "days"),
    linetype = 3
  ) +
  annotate(
    "text",
    x = difftime("2018-09-21", "2018-08-29", units = "days"),
    y = 30,
    label = "Bolting",
    family = "Helvetica"
  ) +
  theme_base() +
  scale_colour_few() +
  scale_fill_few() +
  scale_x_continuous(limits = c(0, 31)) +
  labs(x = "Days after germination", y = "Projected leaf area [cm²]") +
  theme(
    text = element_text(family = "Helvetica"),
    legend.position = c(0.1, 0.8),
    legend.title = element_blank(),
    plot.background = element_blank()
  )

pdf("leafgp_plot.pdf", height = 4, width = 10)
leafgp.area
dev.off()

manual.height <-
  ggplot(data = height.data, aes(x = dag, y = height)) +
  geom_smooth(
    aes(fill = genotype, colour = genotype),
    method = "lm",
    formula = y ~ splines::bs(x, 4),
    span = 0.75,
    se = TRUE,
    size = 0.5,
    linetype = 2
  ) +
  geom_jitter(
    aes(group = interaction(dag, genotype), fill = genotype),
    shape = 21,
    width = 0.25,
    size = 2,
    stroke = 0.25
  ) +
  geom_vline(
    xintercept = difftime("2018-09-22", "2018-08-29", units = "days"),
    linetype = 3
  ) +
  annotate(
    "text",
    x = difftime("2018-09-23", "2018-08-29", units = "days"),
    y = 20,
    label = " Bolting",
    family = "Helvetica"
  ) +
  theme_base() +
  scale_colour_few() +
  scale_fill_few() +
  scale_x_continuous(limits = c(24, 60)) +
  scale_y_continuous(limits = c(0, 60)) +
  labs(x = "Days after germination", y = "Height [cm]") +
  theme(
    text = element_text(family = "Helvetica"),
    legend.position = c(0.1, 0.8),
    legend.title = element_blank(),
    plot.background = element_blank()
  )

pdf("manual_height_plot.pdf",
    height = 4,
    width = 10)
manual.height
dev.off()

leafgp.full <-
  ggplot() +
  geom_smooth(
    data = full.pheno,
    aes(
      fill = genotype,
      colour = genotype,
      x = dag,
      y = Projected_LeafArea.mm.2. / 100
    ),
    method = "lm",
    # method = "loess",
    formula = y ~ splines::bs(x, 4),
    # formula =  y ~ x,
    span = 0.75,
    se = TRUE,
    size = 0.5,
    linetype = 2
  ) +
  geom_jitter(
    data = subset(
      full.pheno,
      date != "2018-09-04" &
        date != "2018-09-11" &
        date != "2018-09-13" &
        date != "2018-09-18" &
        date != "2018-09-27"
    ),
    aes(
      group = interaction(dag, genotype),
      fill = genotype,
      x = dag,
      y = Projected_LeafArea.mm.2. / 100
    ),
    shape = 21,
    width = 0.25,
    size = 2,
    stroke = 0.25
  ) +
  geom_smooth(
    data = full.pheno,
    aes(
      fill = genotype,
      colour = genotype,
      x = dag,
      y = height
    ),
    method = "lm",
    # method = "loess",
    formula = y ~ splines::bs(x, 4),
    # formula =  y ~ x,
    span = 0.75,
    se = TRUE,
    size = 0.5,
    linetype = 1
  ) +
  geom_jitter(
    data = subset(
      full.pheno,
      date != "2018-10-01" &
        date != "2018-10-04" &
        date != "2018-10-11" &
        date != "2018-10-16"
    ),
    aes(
      group = interaction(dag, genotype),
      fill = genotype,
      x = dag,
      y = height
    ),
    shape = 22,
    width = 0.25,
    size = 2,
    stroke = 0.25
  ) +
  geom_vline(
    xintercept = difftime("2018-09-22", "2018-08-29", units = "days"),
    linetype = 3
  ) +
  annotate(
    "text",
    x = difftime("2018-09-20", "2018-08-29", units = "days"),
    y = 30,
    label = "Bolting",
    family = "Helvetica"
  ) +
  theme_base() +
  scale_colour_few(labels = c(
    " Col-0",
    expression(italic(" lac4/10/12/17")),
    expression(italic(" lac5/10/12/17")),
    expression(italic(" lac4/5/10/12/17"))
  )) +
  scale_fill_few(labels = c(
    " Col-0",
    expression(italic(" lac4/10/12/17")),
    expression(italic(" lac5/10/12/17")),
    expression(italic(" lac4/5/10/12/17"))
  )) +
  scale_x_continuous(limits = c(0, 63)) +
  scale_y_continuous(sec.axis = dup_axis(name = "Height [cm]")) +
  labs(x = "Days after germination", y = "Projected leaf area [cm²]") +
  theme(
    text = element_text(family = "Helvetica"),
    legend.position = c(0.1, 0.8),
    legend.title = element_blank(),
    legend.text.align = 0,
    plot.background = element_blank()
  )

pdf("pheno_full.pdf", height = 4, width = 10)
leafgp.full
dev.off()

final.plot <- ggplot(data = final.pheno, aes(x = genotype, y = value)) +
  geom_jitter(width = 0.1, shape = 21) +
  geom_boxplot(width = 0.1, outlier.alpha = 0) +
  theme_base() +
  theme(
    text = element_text(family = "Helvetica"),
    legend.position = c(0.1, 0.8),
    legend.title = element_blank(),
    legend.text.align = 0,
    plot.background = element_blank(),
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),
    axis.title = element_blank()
  ) +
  geom_text(data = pheno.letters, 
            aes(label = groups), 
            angle = 0,
            hjust = 0, 
            family = "Helvetica") +
  scale_x_discrete(labels = rev(c(
    "Col-0",
    expression(italic("lac4/10/12/17")),
    expression(italic("lac5/10/12/17")),
    expression(italic("lac4/5/10/12/17"))
  ))) +
  facet_wrap(. ~ measurement, ncol = 3, scales = "free_x") +
  coord_flip()

pdf("final_pheno.pdf", height = 7, width = 7)
final.plot
dev.off()