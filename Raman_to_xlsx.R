#### select path to the folder containing your spectra (.txt) files ####
spectra_path <- choose.dir(caption = "Choose folder containing your spectra (.txt) files") # works only on windows
# spectra_path <- "/home/leonard/Downloads/Manuela/" # alternatively, choose path manually

#### write function to make sure required packages are installed ####
install_and_load <- function(x) {
  for (i in x) {
    if (!require(i, character.only = TRUE)) {
      install.packages(i, dependencies = TRUE)
      require(i, character.only = TRUE)
    }
  }
}

#### load packages (and install if necessary) ####
install_and_load(c("baseline", "tidyverse", "MESS", "writexl", "zoo"))

#### list Raman spectra files ####
Raman_files <-
  list.files(
    path = spectra_path,
    pattern = "*.txt",
    recursive = TRUE,
    full.names = TRUE
  )

#### write function to read spectra ####
read_Raman <- function(flnm) {
  read_tsv(flnm,
    comment = "#",
    col_names = c("wavenumber", "intensity"),
    col_types = "cc",
  ) %>%
    mutate(
      filename = basename(flnm)
    ) %>%
    mutate(
      wavenumber = as.numeric(str_replace(wavenumber, fixed(","), fixed("."))),
      intensity = as.numeric(str_replace(intensity, fixed(","), fixed("."))),
    ) %>%
    filter(wavenumber > 300) %>%
    pivot_wider(
      id_cols = filename,
      names_from = wavenumber,
      values_from = intensity
    ) %>%
    group_by(filename) %>%
    nest()
}

#### read spectra and combine into a tibble ####
Raman_spectra <- map_dfr(Raman_files, read_Raman)

#### write baseline correction function ####
Raman_als <- function(x) {
  corrected_spectra <- baseline::baseline.als(as.matrix(x),
    lambda = 5,
    p = 0.01,
    maxit = 100
  )
  as_tibble(corrected_spectra$corrected, rownames = NA)
}

#### baseline correction ####
Raman_spectra_corrected <- Raman_spectra %>%
  mutate(corrected_data = map(data, Raman_als)) %>%
  select(-data) %>%
  unnest(cols = c(corrected_data)) %>%
  pivot_longer(-filename,
    names_to = "wavenumber",
    values_to = "corrected_intensity"
  ) %>%
  drop_na() %>%
  mutate(wavenumber = as.numeric(wavenumber))

#### find the cellulose band around 380 cm-1 ####
cellulose_peak <- Raman_spectra_corrected %>%
  filter(round(wavenumber, digits = 0) %in% c(370:390)) %>%
  group_by(filename) %>%
  mutate(
    peak_pos = wavenumber[which.max(corrected_intensity)]
  ) %>%
  select(-wavenumber, -corrected_intensity) %>%
  unique() %>%
  ungroup() %>%
  mutate(peak_pos_ref = peak_pos[1])

#### align spectra at the cellulose band around 380 cm-1 ####
Raman_spectra_aligned <- Raman_spectra_corrected %>%
  left_join(cellulose_peak) %>%
  ungroup() %>%
  mutate(wavenumber = wavenumber + (peak_pos_ref - peak_pos))

#### normalise spectra to AUC ("scaled") and to total lignin ("scaled_lig") ####
Raman_spectra_scaled <- Raman_spectra_aligned %>%
  group_by(filename) %>%
  mutate(
    AUC = MESS::auc(wavenumber, corrected_intensity),
    scaled = corrected_intensity / AUC,
    scaled_lig = corrected_intensity / (
      (corrected_intensity[abs(wavenumber - 1600) == min(abs(wavenumber - 1600))][1] / 2) + corrected_intensity[abs(wavenumber - 1330) == min(abs(wavenumber - 1330))][1]
    )
  )

#### view spectra to check their alignment ####
ggplot(
  Raman_spectra_scaled,
  aes(
    x = wavenumber,
    y = scaled,
    group = filename
  )
) +
  geom_vline(xintercept = c(380, 1330, 1600)) +
  geom_line()

#### fill missing values by linear interpolation from adjacent values ####
Raman_spectra_filled <- Raman_spectra_scaled %>%
  mutate(wavenumber = round(wavenumber, digits = 0)) %>%
  complete(wavenumber = seq(300, 1800)) %>%
  arrange(filename, wavenumber) %>%
  group_by(filename) %>%
  mutate(raw_full = zoo::na.approx(corrected_intensity, na.rm = F),
         scaled_full = zoo::na.approx(scaled, na.rm = F),
         scaled_lig_full = zoo::na.approx(scaled_lig, na.rm = F))

#### view example spectrum to check that the interpolated values are correct ####
ggplot(Raman_spectra_filled %>% filter(filename == Raman_spectra_filled$filename[1])) +
  geom_line(aes(
    x = wavenumber,
    y = scaled
  ),
  colour = "blue",
  size = 0.2
  ) +
  geom_line(aes(
    x = wavenumber,
    y = scaled_full
  ),
  colour = "red",
  size = 0.2
  )

#### collect raw (unscaled) spectra ####
export_raw <- Raman_spectra_filled %>%
  select(filename, wavenumber, raw_full) %>%
  filter(wavenumber %in% c(300:1800)) %>% 
  pivot_wider(
    id_cols = filename,
    names_from = wavenumber,
    values_from = raw_full
  ) %>%
  select(filename, order(as.numeric(colnames(.))))

#### collect spectra normalised to AUC ####
export_scaled <- Raman_spectra_filled %>%
  select(filename, wavenumber, scaled_full) %>%
  filter(wavenumber %in% c(300:1800)) %>% 
  pivot_wider(
    id_cols = filename,
    names_from = wavenumber,
    values_from = scaled_full
  ) %>%
  select(filename, order(as.numeric(colnames(.))))

#### collect spectra normalised to total lignin ####
export_scaled_lig <- Raman_spectra_filled %>%
  select(filename, wavenumber, scaled_lig_full) %>%
  filter(wavenumber %in% c(300:1800)) %>% 
  pivot_wider(
    id_cols = filename,
    names_from = wavenumber,
    values_from = scaled_lig_full
  ) %>%
  select(filename, order(as.numeric(colnames(.))))

#### export spectra to an .xlsx file with three sheets ####
write_xlsx(list(
  "Raw spectra" = export_raw,
  "AUC scaled spectra" = export_scaled,
  "Lignin scaled spectra" = export_scaled_lig
),
path = paste0(spectra_path, "Raman_spectra.xlsx"),
col_names = TRUE,
format_headers = TRUE,
use_zip64 = FALSE
)