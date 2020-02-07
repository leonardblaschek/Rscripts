library(baseline)
library(tidyverse)

rspec <-
  read_csv(
    "/home/leonard/Documents/Uni/Phloroglucinol/18-06-25_RAMAN/spectra_to_scale.csv"
  ) %>%
  pivot_longer(cols = -wavenumber, names_to = "sample", values_to = "intensity") %>%
  separate(sample, into = c("genotype", "cell.type", "replicate")) %>%
  filter(wavenumber > 600 & wavenumber < 1700)

rspec_wide <- rspec %>%
  group_by(genotype, cell.type, replicate) %>%
  mutate(group_id = row_number()) %>%
  spread(wavenumber, intensity) %>%
  select(-group_id) %>%
  summarise_all(funs(sum(., na.rm = TRUE)))

rspec_mat <- as.matrix(rspec_wide[, -c(1:3)])
rownames(rspec_mat) <- paste(rspec_wide$genotype,
                                rspec_wide$cell.type,
                                rspec_wide$replicate,
                                sep = "_"
)

rspec_correction <- baseline.als(rspec_mat,
                                    lambda = 5,
                                    p = 0.01,
                                    maxit = 100
)

rspec_corrected <- as_tibble(rspec_correction$corrected, rownames = NA) %>%
  rownames_to_column() %>%
  separate(rowname, sep = "_", into = c("genotype", "cell.type", "replicate"), remove = TRUE) %>%
  gather(key = "wavenumber", value = "corrected.intensity", -c(genotype, cell.type, replicate)) %>%
  mutate(wavenumber = as.numeric(wavenumber)) %>%
  group_by(genotype, cell.type, replicate) %>%
  mutate(
    AUC = MESS::auc(wavenumber, corrected.intensity),
    scaled = corrected.intensity / AUC
  )

pdf("check_spectra.pdf")
ggplot(rspec_corrected, 
       aes(x = wavenumber, 
           y = scaled, 
           colour = genotype, 
           group = interaction(genotype, replicate))) +
  geom_line(size = 0.2) +
  facet_wrap(~cell.type, nrow = 3)
dev.off()

pdf("check_spectra_raw.pdf")
ggplot(rspec, 
       aes(x = wavenumber, 
           y = intensity, 
           colour = genotype, 
           group = interaction(genotype, replicate))) +
  geom_line(size = 0.2) +
  facet_wrap(~cell.type, nrow = 3)
dev.off()

rspec_corrected_wide <- rspec_corrected %>%
  pivot_wider(id_cols = wavenumber, values_from = scaled, names_from = c(genotype, cell.type, replicate))

write_csv(rspec_corrected_wide, "raman_corrected_scaled.csv")