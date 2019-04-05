library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(ggthemes)
library(showtext)

#### import Helvetica Neue ####
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
      strip.text = element_text(hjust = 0, face = "italic"),
      axis.ticks = element_line(
        size = 0.25,
        lineend = "square",
        color = "black"
      ),
      axis.text.x = element_text(colour = "black", # flipped coords
                                 margin = margin(1, 1, 1, 1)),
      axis.text.y = element_text(
        colour = "black",
        angle = 0,
        vjust = 0.5,
        hjust = 1,
        margin = margin(1, 1, 1, 1)
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", size = 0.25),
      panel.spacing = unit(1.5, "mm"),
      legend.position = "bottom",
      legend.text = element_text(size = rel(0.8)),
      legend.key.height = unit(4, "mm"),
      complete = TRUE
    )
}

#### read-in function ###
read_plus <- function(flnm) {
  read_tsv(flnm,
           comment = "#",
           col_names = FALSE
  ) %>%
    mutate(filename = basename(flnm))
}

#### load data ###
pheno.files <-
  list.files(
    path = "/home/leonard/Documents/Uni/PhD/Phenotyping/2019-03-26_Laccase_mutants/Trays/",
    pattern = "*.txt",
    recursive = FALSE,
    full.names = TRUE
  )

pheno.read <- lapply(pheno.files, read_plus) %>%
  bind_rows()

pheno.data <- pheno.read %>%
  mutate(X1 = recode(X1, "HEADER_WATERSHED" = "WATERSHED HEADER", #recode for easier separation
                     "HEADER_MARKER" = "MARKER_HEADER",
                     "HEADER_SHAPES" = "SHAPES_HEADER",
                     "HEADER_HISTOGRAM" = "HISTOGRAM_HEADER")) %>%
  separate(X1, into = c("Analysis", "Information")) %>% #separate column for data frame splitting
  mutate(filename = str_remove(filename, ".txt")) %>% #clean
  separate(filename, into = c("Date", "Tray", "Img_ID"), sep = "_|\\.") %>% #separate filename into metadata
  mutate(Img_ID = str_remove(Img_ID, "JPG")) %>% #clean
  mutate(Tray = str_remove(Tray, "Tray")) %>% #clean
  group_split(Analysis) #split into separate data frames based on analysis (i.e. HISTOGRAM, SHAPES, etc.)

#### clean up individual data frames ####
pheno.hist <- pheno.data[[1]] %>%
  select_if(~sum(!is.na(.)) > 0)
colnames(pheno.hist)[3:13] <- (pheno.hist)[1, c(3:13)]
pheno.hist <- pheno.hist %>%
  filter(Information != "HEADER") %>%
  mutate_all(funs(str_remove_all(., "\\[|\\]"))) %>%
  separate_rows(c(4:13)) %>%
  select(-Analysis, -Information)

pheno.mark <- pheno.data[[2]] %>%
  select_if(~sum(!is.na(.)) > 0) %>%
  select(-Img_ID)
colnames(pheno.mark)[3:6] <- (pheno.mark)[1, c(3:6)]
pheno.mark <- pheno.mark %>%
  filter(Information != "HEADER") %>%
  select(-Analysis, -Information)

pheno.shap <- pheno.data[[3]] %>%
  select_if(~sum(!is.na(.)) > 0)
colnames(pheno.shap)[3:13] <- (pheno.shap)[1, c(3:13)]
pheno.shap <- pheno.shap %>%
  filter(Information != "HEADER") %>%
  select(-Analysis, -Information)

pheno.wtrs <- pheno.data[[4]] %>%
  select_if(~sum(!is.na(.)) > 0)
colnames(pheno.wtrs)[3] <- (pheno.wtrs)[1, 3]
pheno.wtrs <- pheno.wtrs %>%
  filter(Information != "HEADER") %>%
  select(-Analysis, -Information)

#### ID and merge dataframes ####
tray.ids <- read_csv("/home/leonard/Documents/Uni/PhD/Phenotyping/2019-03-26_Laccase_mutants/Tray_IDs.csv", col_types = "cccc")

pheno.full <- pheno.shap %>%
  left_join(pheno.mark) %>%
  left_join(pheno.wtrs) %>%
  left_join(tray.ids)

#### scale the data to marker size ####
pheno.full <- pheno.full %>%
  group_by(Date) %>%
  mutate(marker_area = max(as.numeric(marker_area)))
pheno.full$scale_factor <- as.numeric(pheno.full$marker_area) / (4^2 * pi)
pheno.full$area <- as.numeric(pheno.full$area) / pheno.full$scale_factor
pheno.full$Date <- as.Date(pheno.full$Date)

#### remove dead plants ####
# pheno.full <- pheno.full %>%
#   unite("Tray_ID", Tray, Img_ID, sep = "_", remove = FALSE) %>%
#   filter(Tray_ID != "5_5" & Tray_ID != "5_6" & Tray_ID != "5_14" & Tray_ID != "4_0" & Tray_ID != "4_4")

### plot rosette growth kinetics ###
rosette.kinetics <- ggplot(data = pheno.full, aes(x = Date, y = area, fill = Genotype)) +
  geom_jitter(aes(size = estimated_object_count), shape =21, alpha = 0.75, stroke = 0.25) +
  geom_line(aes(group = interaction(Plant, Genotype), colour = Genotype), alpha = 0.15, size = 0.25) +
  # geom_boxplot(size = 0.5, width = 0.5, position = position_dodge(width = 1), aes(group = interaction(Genotype, Date))) +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  theme_leo() +
  theme(legend.position = c(0.15, 0.8),
        legend.title = element_blank())

pdf("rosette_area.pdf")
rosette.kinetics
dev.off()

rosette.kinetics <- ggplot(data = pheno.full, aes(x = Date, y = area)) +
  geom_line(aes(group = interaction(Plant, Genotype), colour = Genotype), alpha = 0.15, size = 0.25) +
  stat_summary(aes(group = Genotype, colour = Genotype), geom = "line", fun.y = "mean") +
  scale_colour_viridis_d() +
  labs(y = expression(Rosette~area~"[mm"^"2"*"]")) +
  theme_leo() +
  theme(legend.position = c(0.1, 0.9),
        legend.title = element_blank())

pdf("rosette_area_grouped.pdf")
rosette.kinetics
dev.off()
