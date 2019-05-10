library(reshape2)
library(ggplot2)
library(ggthemes)
library(showtext)
library(dplyr)
library(cowplot)
library(agricolae)

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
        # face = "italic"
        ),
      axis.ticks = element_line(
        size = 0.25,
        lineend = "square",
        color = "black"
      ),
      axis.title = element_blank(),
      axis.text.y = element_text(colour = "black"),
      axis.text.x = element_text(
        colour = "black",
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", size = 0.25),
      panel.spacing = unit(1.5, "mm"),
      legend.position = "bottom",
      legend.text = element_text(size = rel(0.8)),
      legend.key.height = unit(4, "mm"),
      legend.key.width = unit(30, "mm"),
      # plot.margin = unit(c(0, 0, 0, 0), "cm"),   
      
      complete = TRUE
    )
}


###############################
# functions for scaling and statistics
###############################
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

tukey <- function(x) {
  aov1 <- aov(data = x, value ~ genotype)
  groups <- HSD.test(aov1, "genotype", alpha = 0.05)
  groups$groups[["genotype"]] <- rownames(groups$groups)
  groups$groups[["value"]] <- 0
  return(groups[["groups"]])
}

###############################
# read data, create descriptive data frame
###############################
irx.data <-
  read.csv("file:///home/leonard/Documents/Uni/PhD/IRX/irx_measurements.csv",
           sep = "\t")

irx.data$genotype <- recode(irx.data$genotype, "col-0" = "Col-0")

###############################
# create reference at the heigth of the cambium for distance calculation
###############################
irx.ref <-
  subset(irx.data, object == "ref", select = c(1, 2, 3, 6, 7, 27))
irx.ref$ref.y1 <- irx.ref$Y
irx.ref$ref.y2 <- irx.ref$Y
irx.ref$ref.x1 <- irx.ref$X - (irx.ref$Length / 2)
irx.ref$ref.x2 <- irx.ref$X + (irx.ref$Length / 2)

###############################
# merge reference into data frame
###############################
irx.data <-
  full_join(
    subset(irx.data, object != "ref"),
    irx.ref[, c(1, 2, 3, 7, 8, 9, 10)],
    by = c("genotype", "replicate", "technical")
  )

###############################
# calculate difference of the centre of each vessel to the cambium
###############################
irx.data$Distance <-
  apply(irx.data[, c("X", "Y", "ref.x1", "ref.x2", "ref.y1", "ref.y2")],
        1 ,
        function(x) {
          a <- c(x[1], x[2])
          b <- c(x[3], x[5])
          c <- c(x[4], x[6])
          v1 <- b - c
          v2 <- a - b
          m <- cbind(v1, v2)
          d <- abs(det(m)) / sqrt(sum(v1 * v1))
          d
        })

###############################
# add a variable to identify each vessel
###############################
irx.data <- tibble::rowid_to_column(irx.data, "number")

###############################
# factorise and order for facetting
###############################
irx.data$replicate <- as.factor(irx.data$replicate)
irx.data$technical <- as.factor(as.character(irx.data$technical))
irx.data$object <-
  ordered(irx.data$object, levels = c("PX", "PMX", "SMX"))

###############################
# create coordinate data frame
###############################
irx.roi <- irx.data

###############################
# create subsets containing the
# X and the Y coordinates of the vessels respectively
###############################
irx.roix <-
  data.frame(irx.roi %>% select(
    genotype,
    replicate,
    technical,
    object,
    number,
    contains("ROIx")
  )) %>%
  melt(
    id = c("genotype", "replicate", "technical", "object", "number"),
    value.name = "ROIx"
  ) %>%
  cbind(., colsplit(.$variable, "x", c("variable", "point")))

irx.roiy <-
  data.frame(irx.roi %>% select(
    genotype,
    replicate,
    technical,
    object,
    number,
    contains("ROIy")
  )) %>%
  melt(
    id = c("genotype", "replicate", "technical", "object", "number"),
    value.name = "ROIy"
  ) %>%
  cbind(., colsplit(.$variable, "y", c("variable", "point")))

###############################
# bind X and Y coordinate data frames together and remove superfluous rows
###############################
irx.roi <- cbind(irx.roix[, c(1:5, 7, 9)], irx.roiy[, 7])
colnames(irx.roi)[8] <- "ROIy"

irx.roi <-
  irx.roi[apply(irx.roi[c(6, 8)], 1, function(z) {
    !any(z == 0)
  }), ]
irx.roi <-
  irx.roi[apply(irx.roi[c(6, 8)], 1, function(z) {
    !any(is.na(z))
  }), ]

###############################
# clean descriptive data for each vessel
###############################
irx.data$HeightByWidth <- irx.data$Height / irx.data$Width
irx.data <-
  irx.data[, c(
    "genotype",
    "replicate",
    "technical",
    "object",
    "number",
    "Area",
    "Perim.",
    "Width",
    "Height",
    "HeightByWidth",
    "Circ.",
    "Distance"
  )]

###############################
# melt data frame for facetting
###############################
irx.melt <-
  melt(irx.data, id = c("genotype", "replicate", "technical", "object"))

irx.melt$genotype <-
  ordered(
    irx.melt$genotype,
    levels = c(
      "Col-0",
      "4cl1",
      "4cl2",
      "4cl1x4cl2",
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

irx.melt <- filter(irx.melt, genotype != "ccr1xfah1" & genotype != "cad4")
###############################
# calculate relative distance, where the innermost PX coordinate
# denotes the end (1) of the vascular bundle
###############################
irx.data <- irx.data %>%
  group_by(genotype, replicate, technical) %>%
  mutate(rel.distance = Distance / max(Distance, na.rm = TRUE))

###############################
# join descriptive data and coordinate data and order the data frame
###############################
irx.roi <-
  full_join(irx.data,
            irx.roi,
            by = c("genotype", "replicate", "technical", "object", "number"))
irx.roi <-
  irx.roi[order(
    irx.roi$genotype,
    irx.roi$replicate,
    irx.roi$technical,
    irx.roi$number,
    irx.roi$point
  ),]

###############################
# scale the vessels for plotting on continuous scales
###############################
irx.roi <- irx.roi %>%
  group_by(genotype, replicate, technical, object, number) %>%
  mutate(ROIx = (Circ. + ((
    ROIx - mean(ROIx, na.rm = TRUE)
  )) / 2500),
  ROIy = (rel.distance + ((
    ROIy - mean(ROIy, na.rm = TRUE)
  )) / 2500))

###############################
# calculate z-scores
###############################
irx.melt <- irx.melt %>%
  group_by(variable) %>%
  mutate(value.scaled = scale_this(value))

###############################
# order from innermost to outermost cell type
###############################
irx.roi$object <-
  ordered(irx.roi$object, levels = c("PX", "PMX", "SMX"))

###############################
# remove number variable for overview plot
###############################
irx.melt <- subset(irx.melt, variable != "number")

###############################
# export data
###############################
write_csv(irx.data, "Wiesner_IRX_data.csv")

###############################
# Tukey-HSD test 
###############################
irx.letters <- filter(irx.melt, variable %in% c("HeightByWidth")) %>%
  group_by(object, variable) %>%
  do(data.frame(tukey(.)))

###############################
# shape and position overview plot
###############################
irx.overview <-
  ggplot(data = filter(irx.melt, variable %in% c("HeightByWidth")), aes(x = genotype, y = value)) +
  geom_jitter(
    aes(fill = value.scaled),
    shape = 21,
    width = 0.1,
    alpha = 0.9,
    size = 2,
    stroke = 0.25
  ) +
  # geom_violin(draw_quantiles = 0.5, adjust = 1.5, fill = rgb(1,1,1,0.5)) +
  geom_boxplot(fill = rgb(1,1,1,0.5), outlier.alpha = 0) +
  # geom_label(data = irx.letters, 
  #            aes(label = groups), 
  #            angle = 90,
  #            fill = rgb(1, 1, 1, 0.75), 
  #            hjust = 1, 
  #            label.size = 0, 
  #            family = "Helvetica") +
  geom_text(data = irx.letters, 
             aes(label = groups), 
             angle = 90,
             hjust = 1, 
             family = "Helvetica") +
  scale_fill_distiller(palette = "RdBu", name = "Z-score by\nrow") +
  scale_y_continuous(expand = expand_scale(mult = c(0.28,0.05))) +
  scale_x_discrete(
    labels = c(
      "Col-0",
      expression(italic("4cl1")),
      expression(italic("4cl2")),
      expression(paste(italic("4cl1"), "x", italic("4cl2"))),
      expression(italic("ccoaomt1")),
      expression(italic("fah1")),
      expression(italic("omt1")),
      expression(italic("ccr1")),
      expression(italic("cad5")),
      expression(paste(italic("cad4"), "x", italic("cad5")))
    )
  ) +
  theme_leo() +
  theme(legend.position = "none") +
  facet_grid(variable ~ object, scales = "free_y") 
  # coord_flip()

pdf("irx_overview_4.pdf", width = 10, height = 6)
irx.overview
dev.off()