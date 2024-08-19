## ----working_directory------------------------------------------------------------
my.wd <- "~/projects/local-files/soilspec_training"
# Or within an RStudio project
# my.wd <- getwd()


## ----setup, message=FALSE, warning=FALSE------------------------------------------
library("tidyverse")
library("asdreader")
library("opusreader2")
library("prospectr")
library("qs")
library("moments")
library("resemble")


## ----asdreader, message=FALSE, warning=FALSE--------------------------------------
# Downloading an .asd file
visnir.spectra.url <- "https://github.com/soilspectroscopy/ossl-models/raw/main/sample-data/101453MD01.asd"

visnir.spectra.path <- file.path(my.wd, "file1.asd")

download.file(url = visnir.spectra.url,
              destfile = visnir.spectra.path,
              mode = "wb")

# Reading asd file
visnir.spectra <- asdreader::get_spectra(visnir.spectra.path)

# Inspecting the file
class(visnir.spectra)
dim(visnir.spectra)
visnir.spectra[1,1:5]

# Spectral range
range(as.numeric(colnames(visnir.spectra)))


## ----opusreader2, message=FALSE, warning=FALSE------------------------------------
# Downloading an .0 file
mir.spectra.url <- "https://github.com/soilspectroscopy/ossl-models/raw/main/sample-data/235157XS01.0"

mir.spectra.path <- file.path(my.wd, "file2.0")

download.file(url = mir.spectra.url,
              destfile = mir.spectra.path,
              mode = "wb")

# Reading asd file
mir.spectra <- opusreader2::read_opus_single(dsn = mir.spectra.path)

# Inspecting the file
class(mir.spectra)
names(mir.spectra)

# Spectra is stored in file$ab$data
class(mir.spectra$ab$data)
dim(mir.spectra$ab$data)

# Spectral range
range(as.numeric(colnames(mir.spectra$ab$data)))


## ----read_csv, message=FALSE, warning=FALSE---------------------------------------
# Downloading an csv output from Neospectra
nir.spectra.url <- "https://github.com/soilspectroscopy/ossl-models/raw/main/sample-data/sample_neospectra_data.csv"

nir.spectra.path <- file.path(my.wd, "file3.csv")

download.file(url = nir.spectra.url,
              destfile = nir.spectra.path,
              mode = "wb")

# Reading csv file
nir.spectra <- readr::read_csv(nir.spectra.path)

# Inspecting the file
class(nir.spectra)
nir.spectra[1:5,1:5]

# Spectral range after removing first column
range(as.numeric(colnames(nir.spectra[,-1])))


## ----nir_reflectance--------------------------------------------------------------
# Original data
nir.spectra[1:5,1:5]

# Spectra column names
spectra.column.names <- nir.spectra %>%
  select(-sample_id) %>%
  names()

# Transforming reflectance (%) to reflectance factor (decimal, 0-1)
# Also, rounding to 5 decimal places of precision
nir.spectra.rf <- nir.spectra %>%
  mutate(across(all_of(spectra.column.names), ~round(.x/100, 5)))

nir.spectra.rf[1:5,1:5]


## ----visualization----------------------------------------------------------------
## Pivot to long format
nir.spectra.rf.long <- nir.spectra.rf %>%
  pivot_longer(all_of(spectra.column.names),
               names_to = "wavelength",
               values_to = "reflectance") %>%
  mutate(wavelength = as.numeric(wavelength),
         reflectance = as.numeric(reflectance))

head(nir.spectra.rf.long)

## Visualization
ggplot(data = nir.spectra.rf.long) +
  geom_line(aes(x = wavelength, y = reflectance,
                group = sample_id),
            alpha = 0.5, linewidth = 0.5) +
  theme_light()


## ----interpolation_columns--------------------------------------------------------
# Old columns, reversed, and as numeric
old.wavelength <- as.numeric(rev(spectra.column.names))
head(old.wavelength)

# New columns, increasing order, spaced 2 nm
new.wavelength <- seq(1350, 2550, by = 2)
head(new.wavelength)


## ----interpolation----------------------------------------------------------------
# Selecting old spectra in increasing order
# Parse to matrix (input of prospectr::resample)
# Resample
# Parse to tibble
# Bind to original sample ids
nir.spectra.rf.int <- nir.spectra.rf %>%
  select(all_of(rev(spectra.column.names))) %>%
  as.matrix() %>%
  prospectr::resample(X = .,
                      wav = old.wavelength,
                      new.wav = new.wavelength,
                      interpol = "spline") %>%
  as_tibble() %>%
  bind_cols({nir.spectra.rf %>%
      select(sample_id)}, .)

nir.spectra.rf.int[1:5,1:5]


## ----visualization_resample-------------------------------------------------------
new.spectra.column.names <- as.character(new.wavelength)

nir.spectra.rf.int %>%
  pivot_longer(all_of(new.spectra.column.names),
               names_to = "wavelength",
               values_to = "reflectance") %>%
  mutate(wavelength = as.numeric(wavelength),
         reflectance = as.numeric(reflectance)) %>%
  ggplot(data = .) +
  geom_line(aes(x = wavelength, y = reflectance, group = sample_id),
            alpha = 0.5, linewidth = 0.5) +
  theme_light()


## ----sg---------------------------------------------------------------------------
# Select spectra columns
# Parse to matrix (input of prospectr::savitzkyGolay)
# Apply preprocessing
# Parse to tibble
# Bind the spectra to id column
nir.spectra.sg <- nir.spectra.rf.int %>%
  select(all_of(new.spectra.column.names)) %>%
  as.matrix() %>%
  savitzkyGolay(X = ., p = 2, w = 11, m = 1, delta.wav = 2) %>%
  as_tibble() %>%
  bind_cols({nir.spectra.rf %>%
      select(sample_id)}, .)

nir.spectra.sg[1:5,1:5]


## ----sg_visualization-------------------------------------------------------------
nir.spectra.sg %>%
  pivot_longer(any_of(new.spectra.column.names),
               names_to = "wavelength",
               values_to = "reflectance") %>%
  mutate(wavelength = as.numeric(wavelength),
         reflectance = as.numeric(reflectance)) %>%
  ggplot(data = .) +
  geom_line(aes(x = wavelength, y = reflectance, group = sample_id),
            alpha = 0.5, linewidth = 0.5) +
  theme_light()


## ----snv--------------------------------------------------------------------------
# Select spectra columns
# Parse to matrix (input of prospectr::savitzkyGolay)
# Apply preprocessing
# Parse to tibble
# Bind the spectra to id column
nir.spectra.snv <- nir.spectra.rf.int %>%
  select(all_of(new.spectra.column.names)) %>%
  as.matrix() %>%
  prospectr::standardNormalVariate(X = .) %>%
  as_tibble() %>%
  bind_cols({nir.spectra.rf %>%
      select(sample_id)}, .)

nir.spectra.sg[1:5,1:5]


## ----snv_visualization------------------------------------------------------------
nir.spectra.snv %>%
  pivot_longer(any_of(new.spectra.column.names),
               names_to = "wavelength",
               values_to = "reflectance") %>%
  mutate(wavelength = as.numeric(wavelength),
         reflectance = as.numeric(reflectance)) %>%
  ggplot(data = .) +
  geom_line(aes(x = wavelength, y = reflectance, group = sample_id),
            alpha = 0.5, linewidth = 0.5) +
  theme_light()


## ----neospectra_read, eval=TRUE---------------------------------------------------
## Internet configuration for downloading big datasets
options(timeout = 10000)

## Reading serialized files
neospectra.soil <- qread_url("https://storage.googleapis.com/soilspec4gg-public/neospectra_soillab_v1.2.qs")
dim(neospectra.soil)
neospectra.soil[1:5,1:5]

neospectra.site <- qread_url("https://storage.googleapis.com/soilspec4gg-public/neospectra_soilsite_v1.2.qs")
dim(neospectra.site)
neospectra.site[1:5,1:5]

neospectra.nir <- qread_url("https://storage.googleapis.com/soilspec4gg-public/neospectra_nir_v1.2.qs")
dim(neospectra.nir)
neospectra.nir[1:5,1:5]


## ----neospectra_filter, eval=TRUE-------------------------------------------------
# Column names
neospectra.site %>%
  names()

neospectra.soil %>%
  names()

neospectra.nir %>%
  select(1:20) %>%
  names()

# How many distinct id.sample_local_c in the NIR file?
neospectra.nir %>%
  distinct(id.sample_local_c) %>%
  nrow()

# How many distinct countries and samples from?
# 2016 samples from USA, other from African countries
# We can use this columns to split and block datasets
neospectra.site %>%
  count(location.country_iso.3166_txt)

# Selecting relevant site data
neospectra.site <- neospectra.site %>%
  select(id.sample_local_c, location.country_iso.3166_txt)

# Selecting relevant soil data
neospectra.soil <- neospectra.soil %>%
  select(id.sample_local_c, oc_usda.c729_w.pct)

# Selecting relevant NIR data and taking average across repeats
neospectra.nir <- neospectra.nir %>%
  select(id.sample_local_c, starts_with("scan_nir")) %>%
  group_by(id.sample_local_c) %>%
  summarise_all(mean, .group = "drop")

# Renaming spectral columns headers to numeric integers
neospectra.nir %>%
  select(starts_with("scan_nir")) %>%
  names() %>%
  head()

old.names <- neospectra.nir %>%
  select(starts_with("scan_nir")) %>%
  names()

new.names <- gsub("scan_nir.|_ref", "", old.names)

neospectra.nir <- neospectra.nir %>%
  rename_with(~new.names, all_of(old.names))

spectral.column.names <- new.names

neospectra.nir[1:5,1:5]


## ----neospectra_view, eval=TRUE---------------------------------------------------
neospectra.nir %>%
  pivot_longer(any_of(spectral.column.names),
               names_to = "wavelength",
               values_to = "reflectance") %>%
  mutate(wavelength = as.numeric(wavelength),
         reflectance = as.numeric(reflectance)) %>%
  ggplot(data = .) +
  geom_line(aes(x = wavelength, y = reflectance, group = id.sample_local_c),
            alpha = 0.25, linewidth = 0.25) +
  theme_light()


## ----neospectra_snv, eval=TRUE----------------------------------------------------
neospectra.nir.snv <- neospectra.nir %>%
  select(all_of(spectral.column.names)) %>%
  as.matrix() %>%
  prospectr::standardNormalVariate(X = .) %>%
  as_tibble() %>%
  bind_cols({neospectra.nir %>%
      select(id.sample_local_c)}, .)

neospectra.nir.snv %>%
  pivot_longer(any_of(spectral.column.names),
               names_to = "wavelength",
               values_to = "reflectance") %>%
  mutate(wavelength = as.numeric(wavelength),
         reflectance = as.numeric(reflectance)) %>%
  ggplot(data = .) +
  geom_line(aes(x = wavelength, y = reflectance, group = id.sample_local_c),
            alpha = 0.25, linewidth = 0.25) +
  theme_light()


## ----neospectra_join, eval=TRUE---------------------------------------------------
# Joining data
neospectra <- left_join(neospectra.site,
                        neospectra.soil,
                        by = "id.sample_local_c") %>%
  left_join(., neospectra.nir.snv, by = "id.sample_local_c")

# Filtering out samples without SOC values
neospectra <- neospectra %>%
  filter(!is.na(oc_usda.c729_w.pct))

neospectra[1:5,1:5]

# Preparing train and test split
neospectra.train <- neospectra %>%
  filter(location.country_iso.3166_txt == "USA")

neospectra.test <- neospectra %>%
  filter(location.country_iso.3166_txt != "USA")


## ----compression, eval=TRUE-------------------------------------------------------
train.pca <- neospectra.train %>%
  select(all_of(spectral.column.names)) %>%
  as.matrix() %>%
  resemble::ortho_projection(Xr = .,
                             Xu = NULL,
                             Yr = NULL,
                             method = "pca",
                             pc_selection = list(method = "cumvar",
                                                 value = 0.99),
                             center = TRUE, scale = TRUE)

names(train.pca)

# How many components where retained?
plot(train.pca, col = "#D42B08CC")

# We can get the scores and exp. variance back and make a plot
train.pca.scores <- train.pca$scores %>%
  as_tibble() %>%
  bind_cols({neospectra.train %>%
      select(id.sample_local_c,
             location.country_iso.3166_txt,
             oc_usda.c729_w.pct)}, .)

train.pca$variance

train.pca.expvar <- round(train.pca$variance$x_var["explained_var",]*100, 2)

p.pca <- ggplot(train.pca.scores) +
  geom_point(aes(x = pc_1, y = pc_2, color = log1p(oc_usda.c729_w.pct)),
             alpha = 0.5, size = 0.5) +
  scale_colour_gradient(low = "gold",
                        high = "darkred",) +
  labs(title = "Neospectra PCA compression",
       x = paste0("PC1 (", train.pca.expvar[1], "%)"),
       y = paste0("PC2 (", train.pca.expvar[2], "%)")) +
  theme_light() +
  theme(legend.position = "bottom"); p.pca

# Predicting the test samples and plotting
test.pca.scores <- neospectra.test %>%
  select(all_of(spectral.column.names)) %>%
  as.matrix() %>%
  predict(train.pca, newdata = .)
  
test.pca.scores <- neospectra.test %>%
    select(id.sample_local_c,
           location.country_iso.3166_txt,
           oc_usda.c729_w.pct) %>%
  bind_cols(test.pca.scores)

# Adding test samples to PCA plot
p.pca.test <- p.pca +
  geom_point(data = test.pca.scores,
             aes(x = pc_1, y = pc_2),
             size = 0.75) +
  labs(subtitle = "Black dots represent testing points"); p.pca.test


## ----save, eval=TRUE--------------------------------------------------------------
# Saving the PCA plot
ggsave(file.path(my.wd, "plot_pca_train_test.png"), p.pca.test,
       dpi = 300, width = 4, height = 3, units = "in", scale = 1.5)

# Saving the train set
write_csv(train.pca.scores, file.path(my.wd, "train.csv"))

# Saving the test set
write_csv(test.pca.scores, file.path(my.wd, "test.csv"))


## ----backtransform, message=FALSE, error=FALSE------------------------------------
# Getting the necessary parameters for backtransformation
train.loadings <- train.pca$X_loadings
train.center <- train.pca$center
train.scale <- train.pca$scale

# Getting scores
train.scores <- predict(train.pca,
                       newdata = {
                         neospectra.train %>%
                           select(all_of(spectral.column.names)) %>%
                           as.matrix()})

test.scores <- predict(train.pca,
                       newdata = {
                         neospectra.test %>%
                           select(all_of(spectral.column.names)) %>%
                           as.matrix()})

# Backtransforming
train.bt <- train.scores %*% train.loadings
test.bt <- test.scores %*% train.loadings

# Rescaling
train.bt <- sweep(x = train.bt,
                  MARGIN = 2, FUN = "*", STATS = train.scale)

test.bt <- sweep(x = test.bt,
                 MARGIN = 2, FUN = "*", STATS = train.scale)

# Recentering
train.bt <- sweep(x = train.bt,
                  MARGIN = 2, FUN = "+", STATS = train.center)

test.bt <- sweep(x = test.bt,
                 MARGIN = 2, FUN = "+", STATS = train.center)

# Binding id columns
train.bt <- as_tibble(train.bt) %>%
  bind_cols({neospectra.train %>%
      select(-all_of(spectral.column.names))}, .)

test.bt <- as_tibble(test.bt) %>%
  bind_cols({neospectra.test %>%
      select(-all_of(spectral.column.names))}, .)

# Quick visualization
selected.id <- "65140"

plot.data <- bind_rows(
  {neospectra.test %>%
      filter(id.sample_local_c == selected.id) %>%
      mutate(source = "original", .before = 1)},
  {test.bt %>%
      filter(id.sample_local_c == selected.id) %>%
      mutate(source = "compressed", .before = 1)})

plot.data %>%
  pivot_longer(any_of(spectral.column.names),
               names_to = "wavelength",
               values_to = "reflectance") %>%
  mutate(wavelength = as.numeric(wavelength),
         reflectance = as.numeric(reflectance)) %>%
  ggplot(data = .) +
  geom_line(aes(x = wavelength, y = reflectance,
                group = source, color = source),
            alpha = 0.5, linewidth = 1) +
  theme_light() + theme(legend.position = "bottom")


## ----q_stats, message=FALSE, error=FALSE------------------------------------------
neospectra.train.spectra <- neospectra.train %>%
  arrange(id.sample_local_c) %>%
  select(all_of(spectral.column.names)) %>%
  as.matrix()

train.bt.spectra <- train.bt %>%
  arrange(id.sample_local_c) %>%
  select(all_of(spectral.column.names)) %>%
  as.matrix()

neospectra.test.spectra <- neospectra.test %>%
  arrange(id.sample_local_c) %>%
  select(all_of(spectral.column.names)) %>%
  as.matrix()

test.bt.spectra <- test.bt %>%
  arrange(id.sample_local_c) %>%
  select(all_of(spectral.column.names)) %>%
  as.matrix()

# Q stats - sum of squared differences - only for test
test.q.stats <- apply((neospectra.test.spectra-test.bt.spectra)^2,
                       MARGIN = 1, sum)

# Binding q.stats to compressed data
test.pca.scores <- test.pca.scores %>%
  mutate(q_stats = test.q.stats, .before = pc_1)


## ----q_critical, message=FALSE, error=FALSE---------------------------------------
E <- cov(neospectra.train.spectra-train.bt.spectra)
teta1 <- sum(diag(E))^1
teta2 <- sum(diag(E))^2
teta3 <- sum(diag(E))^3
h0 <- 1-((2*teta1*teta3)/(3*teta2^2))
Ca <- 2.57 # 1% significance level
Qa <- teta1*(1-(teta2*h0*((1-h0)/teta1^2))+((sqrt(Ca*(2*teta2*h0^2)))/teta1))^(1/h0)
Qa


## ----q_flag, message=FALSE, error=FALSE-------------------------------------------
test.pca.scores <- test.pca.scores %>%
  mutate(represented = ifelse(q_stats <= Qa, TRUE, FALSE), .after = q_stats)

p.pca +
  geom_point(data = test.pca.scores,
             aes(x = pc_1, y = pc_2, fill = represented),
             shape = 21, size = 1.5) +
  labs(fill = "Represented?")

