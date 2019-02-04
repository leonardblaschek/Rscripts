library(dplyr)
library(readr)
library(stringr)
library(tidyr)

raman.files <- list.files(path = "/home/leonard/Documents/Uni/PhD/IRX/RAMAN", 
                          pattern = "*.txt", 
                          recursive = TRUE,
                          full.names = TRUE)

read_plus <- function(flnm) {
  read_tsv(flnm, 
           comment = "#",
           col_names = FALSE,
           skip = 1) %>% 
    mutate(filename = flnm) 
}

raman.data <- lapply(raman.files, read_plus) %>% 
  bind_rows()  

raman.data <- raman.data %>%
  select(c(1:3)) %>%
  filter(!str_detect(filename, "baseline corrected")) 

raman.data$filename <- basename(raman.data$filename)

raman.data <- raman.data %>%
  separate(filename, into = c("genotype", "replicate", "cell.type","technical"), 
           sep = "\\s*#|-|\\s",
           extra = "merge")

raman.data$technical <- str_remove(raman.data$technical, ".txt")