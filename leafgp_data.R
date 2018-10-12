library(ggplot2)
library(ggthemes)
library(showtext)

# import Helvetica Neue
font_add("Helvetica",
         regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf")
showtext_auto()

leafgp.data <-
  read.csv(
    "file:///home/leonard/Documents/Uni/PhD/Phenotyping/18-09_lac_mutants_Processed_1-10-2018_comb-raw_Tray_1/2018-10-01_LeafMeasure.csv"
  )
height.data <-
  read.csv(
    "file:///home/leonard/Documents/Uni/PhD/Phenotyping/18-09_lac_mutants/originals/front/height_measurements.csv",
    skip = 1
  )

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

leafgp.area <-
  ggplot(data = leafgp.data,
         aes(x = dag, y = Projected_LeafArea.mm.2.)) +
  geom_smooth(
    aes(fill = genotype, colour = genotype),
    span = 0.75,
    se = TRUE,
    size = 0.5,
    linetype = 2
  ) +
  geom_point(
    aes(group = interaction(dag, genotype), fill = genotype),
    shape = 21,
    size = 2,
    stroke = 0.25,
    position = position_dodge(width = 0.9)
  ) +
  geom_boxplot(aes(group = interaction(dag, genotype)),
               size = 0.25,
               colour = "black", 
               outlier.alpha = 0,
               position = position_dodge(width = 0.9)) +
  geom_vline(
    xintercept = difftime("2018-09-22", "2018-08-29", units = "days"),
    linetype = 3
  ) +
  annotate(
    "text",
    x = difftime("2018-09-21", "2018-08-29", units = "days"),
    y = 3000,
    label = "Bolting",
    family = "Helvetica"
  ) +
  theme_base() +
  scale_colour_few() +
  scale_fill_few() +
  scale_x_continuous(limits = c(0, 31)) +
  labs(x = "Days after germination", y = "Projected leaf area [mm²]") +
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
    span = 0.75,
    se = TRUE,
    size = 0.5,
    linetype = 2
  ) +
  geom_point(
    aes(group = interaction(dag, genotype), fill = genotype),
    shape = 21,
    size = 2,
    stroke = 0.25,
    position = position_dodge(width = 0.9)
  ) +
  geom_boxplot(aes(group = interaction(dag, genotype)),
               size = 0.25,
               colour = "black", 
               outlier.alpha = 0,
               position = position_dodge(width = 0.9)) +
  theme_base() +
  scale_colour_few() +
  scale_fill_few() +
  scale_x_continuous(limits = c(20, 60)) +
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