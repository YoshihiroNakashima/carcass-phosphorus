################################################################################
# Carcass-derived phosphorus supply: global analysis pipeline
#
# Manuscript: "Conserved in quantity, lost in structure: megafauna loss reshapes
#              the spatial heterogeneity of carcass-derived phosphorus supply"
# Author: Yoshihiro Nakashima (Nihon University)
#
# OVERVIEW
#   This single file reproduces the full analysis, from raw spatial/demographic
#   inputs to every value reported in the manuscript. It is organised into
#   sequential blocks that mirror the original modular scripts (01-07). Each
#   block reads the outputs of earlier blocks, so the file can be run top to
#   bottom (sourcing the whole file) or block by block.
#
#   The estimand is the supply of phosphorus (P) to the site of death. Supply is
#   the SKELETAL P retained at the carcass (soft-tissue P is assimilated by
#   consumers and exported, so it is not local supply). From per-species P-supply
#   fluxes we map (i) the global supply field, (ii) its body-size centroid, and
#   (iii) its spatial heterogeneity (variance-to-mean ratio, VMR), and contrast
#   continents that lost their megafauna (South America, Australia) with one that
#   retained it (Africa), within matched WWF biomes.
#
# DATA INPUTS (see DATA AVAILABILITY in the repository README for access)
#   - IUCN Red List range polygons (registration required; NOT redistributable).
#       Block 1 converts them to a light grid-level occupancy table that IS
#       provided, so Blocks 2+ run without the raw polygons.
#   - MODIS MOD17A3HGF NPP GeoTIFFs (2015-2020).
#   - malddaba demographic database (life-history records).
#   - Etard et al. consolidated trait data (Dryad doi:10.5061/dryad.05qfttfdq).
#   - Greenspoon et al. (2023) per-species global abundances.
#   - WWF Terrestrial Ecoregions of the World (Olson et al. 2001).
#
# KEY INTERMEDIATE OUTPUTS (provided in the repository for reproducibility)
#   cell_species_long.csv      grid-level species occupancy (derived from IUCN)
#   grid_template.csv          0.5-degree grid definition
#   body_mass_master.csv       per-species adult body mass
#   calibration_table.csv      life-table metrics from malddaba
#   predictions_all.csv        extrapolated mortality parameters (all species)
#   turnover_final.csv         per-species annual turnover
#   species_flux_props_NP.csv  per-species P-supply properties (input to 05-07)
#
# REQUIREMENTS
#   R (>= 4.1). Packages: terra, sf, dplyr, tidyr, viridis.
#
# REPRODUCIBILITY NOTES
#   - Random seeds are set where Monte Carlo or model fitting is involved.
#   - The Monte Carlo iteration count (N_ITER) is set to 500 for the reported
#     confidence intervals; reduce only for quick tests.
################################################################################


################################################################################
# CONFIG: all file paths and global constants are defined here.
#         Edit the paths to match your local data layout, then run the blocks.
################################################################################

## ---- Paths (relative; adjust to your environment) --------------------------
DIR_DATA   <- "./data/"     # MODIS NPP GeoTIFFs and malddaba .rdata live here
DIR_OUT    <- "./out/"      # intermediate grid-level outputs are written here
DIR_FIG    <- "./figure/"   # figure outputs

PATH_IUCN_POLYGON <- "MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp"
PATH_NPP_SINGLE   <- "MODIS_NPP_05deg_2015_2020_mean.tif"  # optional single-layer NPP
PATTERN_NPP_TIFF  <- "MOD17A3H_Y_NPP.*\\.(tiff|TIFF)$"     # multi-year NPP tiles
PATH_MALDDABA     <- "data/file.rdata"
DIR_TRAIT         <- "./doi_10_5061_dryad_05qfttfdq__v20260407/"
PATH_TRAIT        <- paste0(DIR_TRAIT, "trait_data_reported.csv")
PATH_GREENSPOON   <- "mammal_biomass-master/results/wild_land_mammal_biomass_inluding_populations.csv"
PATH_ECOREGION    <- "wwf_terr_ecos.shp"  # WWF Terrestrial Ecoregions (Olson et al. 2001)

dir.create(DIR_OUT, showWarnings = FALSE)
dir.create(DIR_FIG, showWarnings = FALSE)

## ---- Grid -------------------------------------------------------------------
RES <- 0.5   # degrees; the analysis grid resolution throughout

## ---- Phosphorus-per-carcass coefficients (literature-based) -----------------
## Dry skeletal mass scales with body mass as M_skel(dry) = A_SKEL * W^B_SKEL
## (Prange et al. 1979). P is an inorganic structural constituent, so skeletal P
## is computed on a dry-skeletal-mass basis (insensitive to water/marrow-fat).
B_SKEL      <- 1.09     # skeletal-mass allometric exponent (Prange et al. 1979)
A_SKEL      <- 0.061    # dry-skeletal-mass coefficient   (Prange et al. 1979)
P_BONE      <- 0.10     # P fraction of dry skeletal mass (McDowell 2003)
P_SOFT_FRAC <- 0.0015   # P fraction of fresh soft-tissue mass (reference only)

## ---- Delivery ---------------------------------------------------------------
## Skeletal P is not removed by predation (the skeleton resists consumption and
## stays at the death site); only local loss to fragmentation, transport and
## weathering is represented, as a uniform, body-size-independent coefficient.
BONE_RETAIN <- 0.8

## ---- Monte Carlo uncertainty (Block 5 and Block 7) --------------------------
## Components are split into shared (common to all species/cells; dominate the
## global CI; cancel in the between-continent contrasts) and independent (per
## species; average down across species).
N_ITER               <- 500    # iterations for reported CIs
SIGMA_BIAS           <- 0.5    # shared: systematic density bias (log scale)
SIGMA_PCONC          <- 0.30   # shared: skeletal P concentration / allometry
SIGMA_TURNOVER_ALLO  <- 0.20   # shared: turnover allometry coefficients
SIGMA_TURNOVER_RESID <- 0.567  # independent: turnover predictive residual SD

## ---- Body-mass filtering ----------------------------------------------------
BM_MIN_G <- 1.0          # exclude < 1 g
BM_MAX_G <- 15.0e6       # exclude > 15 t (terrestrial-mammal sanity bound)


################################################################################
# Shared helper functions
################################################################################

## Build the global 0.5-degree grid with a unique integer id per cell.
make_grid <- function() {
  g <- terra::rast(xmin = -180, xmax = 180, ymin = -90, ymax = 90,
                   resolution = RES, crs = "EPSG:4326")
  terra::values(g) <- seq_len(terra::ncell(g))
  names(g) <- "cell_id"
  g
}

## Read the multi-year MODIS NPP tiles, mask fill/negative values, and return
## the temporal mean as a per-cell table aligned to the analysis grid.
npp_per_cell <- function(grid) {
  npp_files <- list.files(DIR_DATA, pattern = PATTERN_NPP_TIFF,
                          full.names = TRUE, ignore.case = TRUE)
  npp <- terra::rast(npp_files)
  npp[npp >= 99999] <- NA      # fill value (ocean, ice/snow, unvegetated)
  npp[npp < 0]      <- NA
  npp_mean <- terra::mean(npp, na.rm = TRUE)
  data.frame(cell_id = seq_len(terra::ncell(grid)),
             npp = terra::values(npp_mean)[, 1])
}

## Robust column-name matcher (ignores case and non-alphanumeric characters).
find_col <- function(df, candidates) {
  norm <- function(s) gsub("[^a-z0-9]", "", tolower(s))
  nn <- norm(names(df))
  for (c in candidates) {
    h <- which(nn == norm(c))
    if (length(h) > 0) return(names(df)[h[1]])
  }
  NA_character_
}

## Quadrat (cell) area in km^2 at a given latitude band on a sphere.
R_EARTH_KM <- 6371.0
cell_area_km2 <- function(lat) {
  R_EARTH_KM^2 * (pi * RES / 180) *
    (sin((lat + RES / 2) * pi / 180) - sin((lat - RES / 2) * pi / 180))
}


# ----- Block 1: Rasterize IUCN ranges to a grid-level occupancy table
################################################################################
# Convert the IUCN terrestrial-mammal range polygons into a sparse cell-by-
# species occupancy table on the 0.5-degree grid. A cell is "occupied" if a
# range polygon touches it at all (touches = TRUE), which avoids dropping
# species on this coarse grid. Only extant (presence = 1) and native or
# reintroduced (origin in {1, 2}) ranges are used.
#
# This block is the ONLY step that needs the raw IUCN polygons. Its output
# (cell_species_long.csv) is provided in the repository, so all later blocks
# reproduce the results without requiring the polygons.
#
# Input : IUCN range polygons (+ optional single-layer NPP for the npp column)
# Output: cell_species_long.csv, grid_template.csv, range_cell_counts.csv
################################################################################

run_block1_rasterize_ranges <- function() {
  library(terra); library(sf)

  grid <- make_grid()
  xy <- xyFromCell(grid, seq_len(ncell(grid)))
  grid_template <- data.frame(cell_id = seq_len(ncell(grid)),
                              lon = xy[, 1], lat = xy[, 2])

  ## Optional: attach NPP (only needed if you build the occupancy table with an
  ## NPP column here; the main pipeline reads multi-year NPP later regardless).
  if (file.exists(PATH_NPP_SINGLE)) {
    npp <- rast(PATH_NPP_SINGLE)
    npp_aligned <- resample(npp, grid, method = "average")
    grid_template$npp <- values(npp_aligned)[, 1]
  }

  message("Reading IUCN shapefile (large; this can take several minutes)...")
  polys <- st_read(PATH_IUCN_POLYGON, quiet = TRUE)

  ## Column names vary across IUCN versions; detect the species-name column.
  name_col <- if ("binomial" %in% names(polys)) "binomial" else
              if ("sci_name" %in% names(polys)) "sci_name" else names(polys)[1]

  if ("presence" %in% names(polys)) polys <- polys[polys$presence == 1, ]
  if ("origin"   %in% names(polys)) polys <- polys[polys$origin %in% c(1, 2), ]

  polys$sp_name <- polys[[name_col]]
  species_list <- unique(polys$sp_name)
  message(sprintf("Species to rasterize: %d", length(species_list)))

  ## Rasterize each species and record occupied cell ids (memory-efficient loop).
  records <- vector("list", length(species_list))
  for (i in seq_along(species_list)) {
    sp_poly <- polys[polys$sp_name == species_list[i], ]
    r <- rasterize(vect(sp_poly), grid, field = 1, touches = TRUE)
    occupied <- which(!is.na(values(r)))
    if (length(occupied) > 0)
      records[[i]] <- data.frame(cell_id = occupied, binomial = species_list[i])
  }
  cell_species <- do.call(rbind, records)
  cell_species <- merge(cell_species, grid_template, by = "cell_id", all.x = TRUE)

  write.csv(cell_species, paste0(DIR_OUT, "cell_species_long.csv"), row.names = FALSE)
  write.csv(grid_template, paste0(DIR_OUT, "grid_template.csv"), row.names = FALSE)

  range_size <- as.data.frame(table(cell_species$binomial))
  names(range_size) <- c("binomial", "n_cells")
  write.csv(range_size, paste0(DIR_OUT, "range_cell_counts.csv"), row.names = FALSE)
  message("Block 1 done: cell_species_long.csv, grid_template.csv, range_cell_counts.csv")
}


# ----- Block 2a: Body-mass master
################################################################################
# Build the per-species adult body-mass table from the Etard et al. consolidated
# trait data (Dryad doi:10.5061/dryad.05qfttfdq), whose adult_mass_g column
# preserves decimal points correctly. Marine mammals are excluded; only values
# within [1 g, 15 t] are kept.
#
# A small manual supplement adds large species that have IUCN ranges and
# Greenspoon abundances but lack a mass in the trait CSV (e.g. the forest
# elephant Loxodonta cyclotis, sex-averaged adult ~2,700 kg, consistent with the
# savanna elephant L. africana at ~4,565 kg in the same source). Only species
# absent from the master are added; existing values are never overwritten.
#
# Input : Etard trait_data_reported.csv
# Output: body_mass_master.csv (binom, bm_g, bm_source)
################################################################################

run_block2a_body_mass <- function() {
  suppressPackageStartupMessages({ library(dplyr) })

  df <- read.csv(PATH_TRAIT, stringsAsFactors = FALSE, check.names = FALSE)
  sp_col     <- find_col(df, c("iucn2020_binomial", "binomial", "scientificName", "species"))
  bm_col     <- find_col(df, c("adult_mass_g", "body_mass_g", "mass_g"))
  marine_col <- find_col(df, c("marine"))
  if (is.na(sp_col) || is.na(bm_col))
    stop("Could not detect species-name or body-mass column in the trait CSV.")

  master <- data.frame(
    binom  = gsub("_", " ", trimws(as.character(df[[sp_col]]))),
    bm_g   = suppressWarnings(as.numeric(df[[bm_col]])),
    marine = if (!is.na(marine_col)) df[[marine_col]] else 0,
    stringsAsFactors = FALSE)

  master <- master %>%
    filter(!is.na(binom), binom != "", !is.na(bm_g),
           bm_g >= BM_MIN_G, bm_g <= BM_MAX_G,
           is.na(marine) | marine != 1) %>%
    distinct(binom, .keep_all = TRUE) %>%
    mutate(bm_source = "Etard_trait_data_reported") %>%
    select(binom, bm_g, bm_source)

  ## Manual literature supplement (species absent from the trait CSV only).
  mass_supplement <- data.frame(
    binom     = c("Loxodonta cyclotis"),
    bm_g      = c(2.70e6),  # forest elephant, sex-averaged adult ~2,700 kg
    bm_source = c("manual_literature(forest_elephant_mean)"),
    stringsAsFactors = FALSE) %>%
    filter(bm_g >= BM_MIN_G, bm_g <= BM_MAX_G, !(binom %in% master$binom))
  if (nrow(mass_supplement) > 0) master <- bind_rows(master, mass_supplement)

  write.csv(master, "body_mass_master.csv", row.names = FALSE)
  message(sprintf("Block 2a done: body_mass_master.csv (%d species)", nrow(master)))
}


# ----- Block 2b: Life tables from malddaba
################################################################################
# From the raw malddaba records, build a per-species calibration table used to
# train the mortality extrapolations. For each species we derive, from survival
# records, the survivorship curve l_x and then:
#   e0             life expectancy at birth (integral of l_x); basis of turnover
#   mean_age_death mean age of dying individuals
#   q0             juvenile (first-year) mortality
# and, from reproduction (Mx) records, peak annual daughters (mbar).
#
# malddaba records differ in which fields are populated and in array lengths;
# each branch checks that the relevant field is non-empty, equal in length to
# the age vector, and not all-NA, and otherwise skips the record. q0 is
# interpolated only when age = 1 lies within the observed range.
#
# Body mass is joined from body_mass_master.csv (Block 2a) and taxonomic order
# from the same trait CSV used for mass (kept on a single consistent source).
#
# Input : malddaba .rdata, body_mass_master.csv, trait_data_reported.csv (order)
# Output: calibration_table.csv
################################################################################

run_block2b_life_tables <- function() {
  suppressPackageStartupMessages({ library(dplyr) })

  ## Coerce a (possibly NULL / ragged) field to a numeric vector.
  safe_num <- function(x) {
    if (is.null(x) || length(x) == 0) return(NULL)
    suppressWarnings(as.numeric(unlist(lapply(x, function(v) if (is.null(v)) NA else v))))
  }

  ## Decode a PostGIS WKB point (hex string) to c(lon, lat).
  decode_wkb_point <- function(hexstr) {
    if (is.null(hexstr) || is.na(hexstr) || nchar(hexstr) < 42) return(c(NA, NA))
    bytes <- strtoi(substring(hexstr, seq(1, nchar(hexstr) - 1, 2),
                                      seq(2, nchar(hexstr), 2)), 16L)
    raw_b <- as.raw(bytes); little <- raw_b[1] == as.raw(1)
    rd_dbl <- function(off) { b <- raw_b[(off+1):(off+8)]; if (!little) b <- rev(b)
                              readBin(b, "double", size = 8, endian = "little") }
    rd_int <- function(off) { b <- raw_b[(off+1):(off+4)]; if (!little) b <- rev(b)
                              readBin(b, "integer", size = 4, endian = "little") }
    type_int <- rd_int(1); offset <- 5
    if (bitwAnd(type_int, 0x20000000L) != 0) offset <- offset + 4
    c(rd_dbl(offset), rd_dbl(offset + 8))
  }

  ## Survival record -> survivorship curve l_x (NULL if record is ineligible).
  build_lx <- function(rec) {
    if (rec$trait$trait_type != "survival") return(NULL)
    lt <- rec$life_table; age <- safe_num(lt$start_age)
    if (is.null(age) || all(is.na(age))) return(NULL)
    trait <- rec$trait$trait
    if (trait == "Survival_rate") {
      mp <- safe_num(lt$mean_parameter)
      if (is.null(mp) || length(mp) != length(age) || all(is.na(mp))) return(NULL)
      o <- order(age); age <- age[o]; mp <- mp[o]
      sx <- pmin(pmax(mp, 0), 1); sx[is.na(sx)] <- 1
      lx <- c(1, cumprod(sx[-length(sx)]))
    } else if (trait == "N_alive") {
      na <- safe_num(lt$n_alive)
      if (is.null(na) || length(na) != length(age)) return(NULL)
      o <- order(age); age <- age[o]; na <- na[o]
      valid <- which(na > 0 & !is.na(na)); if (length(valid) == 0) return(NULL)
      na <- na[valid[1]:length(na)]; age <- age[valid[1]:length(age)]; lx <- na / na[1]
    } else if (trait == "N_dead") {
      nd <- safe_num(lt$n_dead)
      if (is.null(nd) || length(nd) != length(age) || sum(nd, na.rm = TRUE) == 0) return(NULL)
      o <- order(age); age <- age[o]; nd <- nd[o]; nd[is.na(nd)] <- 0
      dx <- nd / sum(nd); lx <- c(1, 1 - cumsum(dx)[-length(dx)])
    } else return(NULL)
    keep <- !is.na(age) & !is.na(lx); age <- age[keep]; lx <- lx[keep]
    if (length(age) < 2) return(NULL)
    data.frame(age = age, lx = lx)
  }

  ## l_x -> (e0, mean_age_death, q0).
  life_table_metrics <- function(lxdf) {
    age <- lxdf$age; lx <- lxdf$lx
    dage <- diff(c(age, age[length(age)] + median(diff(age))))
    e0 <- sum(lx * dage)
    dx <- -diff(c(lx, 0))
    mad <- if (sum(dx) > 0) sum(age * dx) / sum(dx) else NA
    q0  <- if (min(age) <= 1 && max(age) >= 1) 1 - approx(age, lx, xout = 1)$y else NA
    c(e0 = e0, mean_age_death = mad, q0 = q0)
  }

  load(PATH_MALDDABA)   # provides 'rdata' (a list of records)
  recs <- rdata
  surv_list <- list(); mx_list <- list()
  for (rec in recs) {
    sp <- paste(rec$species$genus, rec$species$species)
    lxdf <- tryCatch(build_lx(rec), error = function(e) NULL)
    if (!is.null(lxdf)) {
      m <- life_table_metrics(lxdf)
      gps <- tryCatch(decode_wkb_point(rec$location$gps_coordinates),
                      error = function(e) c(NA, NA))
      surv_list[[length(surv_list) + 1]] <- data.frame(
        binom = sp, e0 = m["e0"], mean_age_death = m["mean_age_death"],
        q0 = m["q0"], lon = gps[1], lat = gps[2])
    }
    if (!is.null(rec$trait$trait) && rec$trait$trait == "Mx") {
      mxv <- safe_num(rec$life_table$mean_parameter)
      if (!is.null(mxv) && any(!is.na(mxv)))
        mx_list[[length(mx_list) + 1]] <- data.frame(binom = sp, mx_peak = max(mxv, na.rm = TRUE))
    }
  }

  surv_tbl <- bind_rows(surv_list) %>% group_by(binom) %>%
    summarise(e0 = median(e0, na.rm = TRUE),
              mean_age_death = median(mean_age_death, na.rm = TRUE),
              q0 = median(q0, na.rm = TRUE),
              lon = median(lon, na.rm = TRUE), lat = median(lat, na.rm = TRUE),
              .groups = "drop")
  mx_tbl <- bind_rows(mx_list) %>% group_by(binom) %>%
    summarise(mbar = max(mx_peak, na.rm = TRUE), .groups = "drop")

  bm_master <- read.csv("body_mass_master.csv") %>% rename(bm = bm_g)
  tr <- read.csv(PATH_TRAIT, stringsAsFactors = FALSE, check.names = FALSE)
  sp_col_tr    <- find_col(tr, c("iucn2020_binomial", "binomial", "scientificName", "species"))
  order_col_tr <- find_col(tr, c("order"))
  order_tbl <- data.frame(
    binom = gsub("_", " ", trimws(as.character(tr[[sp_col_tr]]))),
    order = trimws(as.character(tr[[order_col_tr]])), stringsAsFactors = FALSE)
  order_tbl <- order_tbl[!is.na(order_tbl$order) & order_tbl$order != "" &
                           !is.na(order_tbl$binom) & order_tbl$binom != "", ]
  order_tbl <- order_tbl[!duplicated(order_tbl$binom), ]

  cal <- surv_tbl %>%
    left_join(mx_tbl, by = "binom") %>%
    left_join(bm_master %>% select(binom, bm), by = "binom") %>%
    left_join(order_tbl, by = "binom") %>%
    filter(!is.na(bm)) %>%
    mutate(logbm = log10(bm))

  write.csv(cal, "calibration_table.csv", row.names = FALSE)
  message(sprintf("Block 2b done: calibration_table.csv (%d species)", nrow(cal)))
}


# ----- Block 2c: Extrapolate mortality parameters to all species
################################################################################
# Train extrapolation models on the calibration table and predict for every
# species in the body-mass master:
#   mean_age_death : log10(mean_age) ~ log10(body mass) + order
#   q0             : logit(q0)       ~ log10(body mass) + order
# Including order materially improves the fit. Species whose order is absent
# from the training set fall back to a body-mass-only model.
#
# Input : calibration_table.csv, body_mass_master.csv, trait_data_reported.csv
# Output: predictions_all.csv (binom, order, bm, logbm, pred_meanage, pred_q0)
################################################################################

run_block2c_extrapolate <- function() {
  suppressPackageStartupMessages({ library(dplyr) })
  set.seed(20260616)

  cal <- read.csv("calibration_table.csv")

  ## mean_age_death model (mass + order), with a mass-only fallback.
  d_shape <- cal %>% filter(is.finite(mean_age_death), is.finite(logbm), !is.na(order)) %>%
    mutate(y = log10(mean_age_death))
  shape_model   <- lm(y ~ logbm + order, data = d_shape)
  shape_bm_only <- lm(y ~ logbm,         data = d_shape)
  shape_orders  <- unique(d_shape$order)

  ## q0 model on the logit scale (mass + order), with a mass-only fallback.
  d_q0 <- cal %>% filter(is.finite(q0), is.finite(logbm), !is.na(order)) %>%
    mutate(q0c = pmin(pmax(q0, 1e-3), 1 - 1e-3), y = log(q0c / (1 - q0c)))
  q0_model   <- lm(y ~ logbm + order, data = d_q0)
  q0_bm_only <- lm(y ~ logbm,         data = d_q0)
  q0_orders  <- unique(d_q0$order)

  allsp <- read.csv("body_mass_master.csv") %>% rename(bm = bm_g) %>%
    filter(bm > 0, !is.na(bm)) %>% mutate(logbm = log10(bm))

  tr <- read.csv(PATH_TRAIT, stringsAsFactors = FALSE, check.names = FALSE)
  sp_col_tr    <- find_col(tr, c("iucn2020_binomial", "binomial", "scientificName", "species"))
  order_col_tr <- find_col(tr, c("order"))
  order_tbl <- data.frame(
    binom = gsub("_", " ", trimws(as.character(tr[[sp_col_tr]]))),
    order = trimws(as.character(tr[[order_col_tr]])), stringsAsFactors = FALSE)
  order_tbl <- order_tbl[!is.na(order_tbl$order) & order_tbl$order != "" &
                           !is.na(order_tbl$binom) & order_tbl$binom != "", ]
  order_tbl <- order_tbl[!duplicated(order_tbl$binom), ]
  allsp <- allsp %>% left_join(order_tbl, by = "binom")

  ## mean_age_death predictions.
  allsp$pred_logmeanage <- NA_real_
  idx1 <- allsp$order %in% shape_orders
  allsp$pred_logmeanage[idx1]  <- predict(shape_model,   newdata = allsp[idx1, ])
  allsp$pred_logmeanage[!idx1] <- predict(shape_bm_only, newdata = allsp[!idx1, ])
  allsp$pred_meanage <- 10^allsp$pred_logmeanage

  ## q0 predictions.
  allsp$pred_q0_logit <- NA_real_
  idx2 <- allsp$order %in% q0_orders
  allsp$pred_q0_logit[idx2]  <- predict(q0_model,   newdata = allsp[idx2, ])
  allsp$pred_q0_logit[!idx2] <- predict(q0_bm_only, newdata = allsp[!idx2, ])
  allsp$pred_q0 <- 1 / (1 + exp(-allsp$pred_q0_logit))

  out <- allsp %>% select(binom, order, bm, logbm, pred_meanage, pred_q0)
  write.csv(out, "predictions_all.csv", row.names = FALSE)
  message(sprintf("Block 2c done: predictions_all.csv (%d species)", nrow(out)))
}


# ----- Block 2d: Turnover (annual population turnover rate)
################################################################################
# Determine per-species turnover (carcass production = abundance x turnover).
# Under a stationary-population assumption:
#   turnover_death = 1 / e0                     (from survivorship)
#   turnover_max   = annual daughters (Mx)      (reproductive ceiling)
#   turnover_final = min(turnover_death, turnover_max)
# e0 is extrapolated from body mass + order (mass-only fallback for unseen
# orders); the reproductive ceiling uses the body-mass allometry of fecundity.
# The ceiling acts as a safeguard and binds only a few small-bodied species.
#
# Input : calibration_table.csv, predictions_all.csv
# Output: turnover_final.csv
################################################################################

run_block2d_turnover <- function() {
  suppressPackageStartupMessages({ library(dplyr) })

  cal  <- read.csv("calibration_table.csv")
  pred <- read.csv("predictions_all.csv")

  ## e0 model (mass + order) with a mass-only fallback.
  d_e0 <- cal %>% filter(is.finite(e0), is.finite(logbm), !is.na(order)) %>%
    mutate(log_e0 = log(pmax(e0, 0.1)))
  e0_model   <- lm(log_e0 ~ logbm + order, data = d_e0)
  e0_bm_only <- lm(log_e0 ~ logbm,         data = d_e0)
  e0_orders  <- unique(d_e0$order)

  ## Reproductive ceiling: log(annual daughters) ~ log(body mass in kg).
  d_mx <- cal %>% filter(is.finite(mbar), mbar > 0, is.finite(logbm)) %>%
    mutate(log_daughters = log(mbar), log_bm_e = log(bm / 1000))
  repro_model <- lm(log_daughters ~ log_bm_e, data = d_mx)
  repro_coef  <- coef(repro_model)

  pred <- pred %>% mutate(bm_kg = bm / 1000, log_bm_e = log(bm_kg), logbm = log10(bm))

  ## Supply a known order for manually added species (e.g. forest elephant) so
  ## the order model can apply; if that order is still unseen in training, the
  ## mass-only fallback is used (no unsupported assumption is introduced).
  order_supplement <- c("Loxodonta cyclotis" = "Proboscidea")
  miss <- is.na(pred$order) & (pred$binom %in% names(order_supplement))
  if (any(miss)) pred$order[miss] <- order_supplement[pred$binom[miss]]

  pred$known <- pred$order %in% e0_orders
  pred$log_e0_pred <- NA_real_
  pred$log_e0_pred[pred$known]  <- predict(e0_model,   newdata = pred[pred$known, ])
  pred$log_e0_pred[!pred$known] <- predict(e0_bm_only, newdata = pred[!pred$known, ])
  pred$e0_pred <- exp(pred$log_e0_pred)

  pred <- pred %>% mutate(
    turnover_death = 1 / pmax(e0_pred, 0.15),
    turnover_max   = exp(repro_coef[1] + repro_coef[2] * log_bm_e),
    turnover_final = pmin(turnover_death, turnover_max))

  out <- pred %>% select(binom, order, bm_kg, turnover_death, turnover_max, turnover_final)
  write.csv(out, "turnover_final.csv", row.names = FALSE)
  message(sprintf("Block 2d done: turnover_final.csv (%d species)", nrow(out)))
}


# ----- Block 3: Density maps (NPP-weighted within-range allocation)
################################################################################
# Build the global mammal density field by distributing each species' fixed
# Greenspoon total abundance across the cells of its IUCN range in proportion to
# NPP:
#   cell abundance = total abundance x (cell NPP / sum of NPP within range)
# More productive cells receive higher density; deserts become low-density
# automatically. Species with no valid NPP cells in range fall back to a uniform
# split. The allocation is purely relative, so the absolute units of NPP do not
# matter, and each species' total is conserved at the Greenspoon value.
#
# Input : Greenspoon abundances, cell_species_long.csv, grid_template.csv, NPP
# Output: cell_density_both.csv, density_npp.tif
################################################################################

run_block3_density_maps <- function() {
  library(terra); library(dplyr)

  grid <- make_grid()
  npp_by_cell <- npp_per_cell(grid)

  cs <- read.csv(paste0(DIR_OUT, "cell_species_long.csv"))
  gp <- read.csv(PATH_GREENSPOON)
  pop <- gp %>% transmute(binomial,
                          pop    = estimated_population,
                          pop_lo = estimated_population_2pt5,
                          pop_hi = estimated_population_97pt5)

  cs <- cs %>% filter(binomial %in% pop$binomial) %>%
    left_join(npp_by_cell, by = "cell_id") %>%
    mutate(npp_filled = ifelse(is.na(npp), 0, npp), has_npp = !is.na(npp)) %>%
    group_by(binomial) %>%
    mutate(npp_sum = sum(npp_filled), n_cells = n(), valid_count = sum(has_npp),
           weight = if_else(valid_count > 0 & npp_sum > 0,
                            if_else(has_npp, npp_filled / npp_sum, 0), 1 / n_cells)) %>%
    ungroup() %>%
    left_join(pop, by = "binomial") %>%
    mutate(ind = pop * weight, ind_lo = pop_lo * weight, ind_hi = pop_hi * weight)

  cell_density <- cs %>% group_by(cell_id) %>%
    summarise(ind = sum(ind), ind_lo = sum(ind_lo), ind_hi = sum(ind_hi),
              n_species = n_distinct(binomial), .groups = "drop") %>%
    left_join(read.csv(paste0(DIR_OUT, "grid_template.csv")), by = "cell_id")

  write.csv(cell_density, paste0(DIR_OUT, "cell_density_both.csv"), row.names = FALSE)
  r <- rast(grid); values(r) <- NA
  r[cell_density$cell_id] <- cell_density$ind
  writeRaster(r, paste0(DIR_OUT, "density_npp.tif"), overwrite = TRUE)
  message("Block 3 done: cell_density_both.csv, density_npp.tif")
}


# ----- Block 4: P-supply flux pipeline (core)
################################################################################
# Compute the per-species phosphorus-supply properties and the cell-level supply
# field. Supply is the SKELETAL P retained at the death site after local loss:
#
#   P_supply_per_carcass = (skeletal P per carcass) x BONE_RETAIN
#
# Skeletal P per carcass is obtained by integrating the age-specific body mass
# W(a) (von Bertalanffy growth) over the age-at-death distribution d(a), so that
# species in which more individuals survive to maturity carry relatively more P.
# Soft-tissue P is computed only as a reference quantity (whole-body P, carcass
# N:P); it is assimilated by consumers and is NOT part of local supply. Nitrogen
# (N) quantities are likewise carried only as a reference foil.
#
# The cell-level supply flux is, per species, (abundance x turnover) x
# P_supply_per_carcass, summed over species; the NPP-weighted within-range
# allocation matches Block 3.
#
# Input : predictions_all.csv, turnover_final.csv, cell_species_long.csv, NPP,
#         Greenspoon abundances
# Output: species_flux_props_NP.csv, cell_flux_NP.csv, cell_flux_NP_perkm2.csv
################################################################################

## ---- Growth and elemental-content models (used by Block 4) -----------------

## von Bertalanffy body mass at age, with simple body-mass-scaled parameters.
vb_weight <- function(age, W_inf, k, b) W_inf * pmax(1 - b * exp(-k * age), 0)^3
growth_params <- function(W_inf_kg) {
  W0 <- 0.06 * W_inf_kg^0.92
  k  <- exp(-0.6 - 0.25 * log(W_inf_kg))
  b  <- 1 - (W0 / W_inf_kg)^(1/3)
  list(W0 = W0, k = k, b = b)
}

## Skeletal (supply) P and reference quantities as functions of body mass (kg).
skeletal_mass <- function(W_kg) A_SKEL * W_kg^B_SKEL
conc_P_bone   <- function(W_kg) skeletal_mass(W_kg) * P_BONE      # supply component
conc_P_soft   <- function(W_kg) W_kg * P_SOFT_FRAC               # reference only
conc_N        <- function(W_kg, dry = 0.30) W_kg * dry * 0.10    # reference foil

## Per-carcass properties: integrate growth x elemental content over the
## age-at-death distribution (first age class = juvenile mortality q0, then
## exponential decay on the mean-age-of-death scale).
species_props <- function(W_inf_g, mean_age, q0) {
  W_inf_kg <- W_inf_g / 1000
  gp <- growth_params(W_inf_kg)
  ages <- seq(0.5, max(mean_age * 3, 2), by = 1.0)
  dx <- numeric(length(ages)); dx[1] <- q0
  if (length(ages) > 1) {
    am <- if (mean_age > 1) mean_age else 1
    w <- exp(-(ages[-1] - 1) / am); w <- w / sum(w) * (1 - q0); dx[-1] <- w
  }
  dx <- dx / sum(dx)
  W_at_age <- vb_weight(ages, W_inf_kg, gp$k, gp$b)
  c(N_per_carcass     = sum(dx * conc_N(W_at_age)),       # reference foil
    Pbone_per_carcass = sum(dx * conc_P_bone(W_at_age)),  # supply (skeletal)
    Psoft_per_carcass = sum(dx * conc_P_soft(W_at_age)),  # reference
    mass_per_carcass  = sum(dx * W_at_age),
    mean_NP_ratio     = sum(dx * conc_N(W_at_age)) /
                        pmax(sum(dx * (conc_P_bone(W_at_age) + conc_P_soft(W_at_age))), 1e-12))
}

run_block4_flux_pipeline <- function() {
  library(terra); library(dplyr)

  grid <- make_grid()
  npp_by_cell <- npp_per_cell(grid)

  pred     <- read.csv("predictions_all.csv")
  turnover <- read.csv("turnover_final.csv")
  stopifnot("bm" %in% names(pred))
  stopifnot(all(c("pred_meanage", "pred_q0") %in% names(pred)))

  props_mat <- t(vapply(seq_len(nrow(pred)), function(i)
    species_props(pred$bm[i], pred$pred_meanage[i], pred$pred_q0[i]),
    FUN.VALUE = c(N_per_carcass = 0, Pbone_per_carcass = 0, Psoft_per_carcass = 0,
                  mass_per_carcass = 0, mean_NP_ratio = 0)))

  props <- data.frame(
    binomial = pred$binom, bm_kg = pred$bm / 1000, pred_meanage = pred$pred_meanage,
    N_per_carcass     = props_mat[, "N_per_carcass"],
    Pbone_per_carcass = props_mat[, "Pbone_per_carcass"],
    Psoft_per_carcass = props_mat[, "Psoft_per_carcass"],
    mass_per_carcass  = props_mat[, "mass_per_carcass"],
    mean_NP_ratio     = props_mat[, "mean_NP_ratio"]) %>%
    mutate(P_per_carcass = Pbone_per_carcass + Psoft_per_carcass) %>%
    left_join(turnover %>% select(binom, turnover_final), by = c("binomial" = "binom")) %>%
    mutate(input_frac_Pbone = BONE_RETAIN,
           P_supply_per_carcass = Pbone_per_carcass * input_frac_Pbone)

  ## Turnover fallback for any species missing from turnover_final.csv: use
  ## 1 / mean_age_death (approximate stationary annual mortality); any remaining
  ## NA is filled conservatively with the median turnover.
  na_turn <- is.na(props$turnover_final) & is.finite(props$pred_meanage) & props$pred_meanage > 0
  props <- props %>% mutate(turnover = ifelse(na_turn, 1 / pred_meanage, turnover_final))
  still_na <- is.na(props$turnover)
  if (any(still_na)) props$turnover[still_na] <- median(props$turnover, na.rm = TRUE)

  write.csv(props %>% select(binomial, bm_kg, N_per_carcass, P_per_carcass,
                             Pbone_per_carcass, Psoft_per_carcass, P_supply_per_carcass,
                             mass_per_carcass, mean_NP_ratio, turnover, input_frac_Pbone),
            "species_flux_props_NP.csv", row.names = FALSE)

  ## Diagnostic: body-size dependence of supply. The realized per-carcass supply
  ## slope is close to 1 (age-structure weighting offsets the mildly superlinear
  ## skeletal allometry); the heterogeneity signal comes from the body-size
  ## distribution and the second-moment form of the VMR index, not from
  ## superlinear stoichiometry.
  fitPskel <- lm(log10(Pbone_per_carcass) ~ log10(bm_kg), data = props[props$Pbone_per_carcass > 0, ])
  message(sprintf("Block 4: per-carcass skeletal-P log-log slope = %.3f", coef(fitPskel)[2]))

  gp <- read.csv(PATH_GREENSPOON)
  pop <- setNames(gp$estimated_population, gp$binomial)

  cs <- read.csv(paste0(DIR_OUT, "cell_species_long.csv")) %>%
    filter(binomial %in% names(pop), binomial %in% props$binomial) %>%
    left_join(npp_by_cell, by = "cell_id") %>%
    mutate(npp_filled = ifelse(is.na(npp), 0, npp), has_npp = !is.na(npp)) %>%
    group_by(binomial) %>%
    mutate(npp_sum = sum(npp_filled), n_cells = n(), valid_count = sum(has_npp),
           weight = if_else(valid_count > 0 & npp_sum > 0,
                            if_else(has_npp, npp_filled / npp_sum, 0), 1 / n_cells)) %>%
    ungroup() %>%
    mutate(ind = pop[binomial] * weight) %>%
    left_join(props, by = "binomial") %>%
    mutate(deaths = ind * turnover,
           flux_N      = deaths * N_per_carcass,                          # reference
           flux_P_body = deaths * P_per_carcass,                          # reference (whole-body P)
           flux_P      = deaths * Pbone_per_carcass * input_frac_Pbone)   # supply (skeletal x retention)

  cell <- cs %>% group_by(cell_id) %>%
    summarise(flux_N = sum(flux_N), flux_P = sum(flux_P), flux_P_body = sum(flux_P_body),
              n_species = n_distinct(binomial),
              massP_weighted_mean = sum(bm_kg * flux_P) / sum(flux_P), .groups = "drop") %>%
    left_join(read.csv(paste0(DIR_OUT, "grid_template.csv")), by = "cell_id") %>%
    mutate(area_km2 = cell_area_km2(lat),
           flux_P_per_km2      = flux_P / area_km2,
           flux_P_body_per_km2 = flux_P_body / area_km2,
           skel_supply_ratio_P = flux_P / pmax(flux_P_body, 1e-12))

  message(sprintf("Block 4: global P supply (skeletal x %.2f) = %.4f Mt P/yr",
                  BONE_RETAIN, sum(cell$flux_P) / 1e9))

  write.csv(cell %>% select(cell_id, lon, lat, flux_N, flux_P, flux_P_body,
                            massP_weighted_mean, skel_supply_ratio_P),
            "cell_flux_NP.csv", row.names = FALSE)
  write.csv(cell %>% select(cell_id, lon, lat, area_km2, flux_P_per_km2, flux_P_body_per_km2),
            "cell_flux_NP_perkm2.csv", row.names = FALSE)
  message("Block 4 done: species_flux_props_NP.csv, cell_flux_NP.csv, cell_flux_NP_perkm2.csv")
}


# ----- Block 5: Uncertainty propagation (Monte Carlo)
################################################################################
# Propagate confidence intervals on the P-supply flux by Monte Carlo, separating
# independent and shared components:
#   independent (per species; average down across species):
#     - per-species density error, back-calculated from the Greenspoon CIs
#     - turnover predictive residual
#   shared (common to all species/cells; dominate the global CI; cancel in the
#   between-continent contrasts):
#     - systematic density bias
#     - skeletal P concentration / allometry
#     - turnover allometry coefficients
# Each iteration draws the shared components once and the independent components
# per species. Because P is not subject to a predation filter, production equals
# supply (the retention coefficient is already in the central estimate).
#
# Input : species_flux_props_NP.csv, cell_species_long.csv, NPP, Greenspoon CIs
# Output: cell_flux_ci_P.csv (per-cell median, bounds, relative uncertainty)
################################################################################

run_block5_uncertainty <- function() {
  library(terra); library(dplyr)
  set.seed(42)

  grid <- make_grid()
  npp_by_cell <- npp_per_cell(grid)

  props <- read.csv("species_flux_props_NP.csv")
  gp <- read.csv(PATH_GREENSPOON) %>%
    mutate(log_sd = (log(pmax(estimated_population_97pt5, 1)) -
                     log(pmax(estimated_population_2pt5, 1))) / (2 * 1.96))
  pop    <- setNames(gp$estimated_population, gp$binomial)
  pop_sd <- setNames(gp$log_sd, gp$binomial)

  cs <- read.csv(paste0(DIR_OUT, "cell_species_long.csv")) %>%
    filter(binomial %in% names(pop), binomial %in% props$binomial) %>%
    left_join(npp_by_cell, by = "cell_id") %>%
    mutate(npp_filled = ifelse(is.na(npp), 0, npp), has_npp = !is.na(npp)) %>%
    group_by(binomial) %>%
    mutate(npp_sum = sum(npp_filled), n_cells = n(), valid_count = sum(has_npp),
           weight = if_else(valid_count > 0 & npp_sum > 0,
                            if_else(has_npp, npp_filled / npp_sum, 0), 1 / n_cells)) %>%
    ungroup() %>%
    left_join(props, by = "binomial") %>%
    mutate(base_flux = weight * pop[binomial] * turnover * P_supply_per_carcass,
           sp_logsd  = ifelse(is.na(pop_sd[binomial]), 0.1, pop_sd[binomial]))

  cells   <- sort(unique(cs$cell_id))
  species <- unique(cs$binomial)
  cs$ci <- match(cs$cell_id, cells)
  cs$si <- match(cs$binomial, species)
  ncells <- length(cells); nsp <- length(species)
  sp_logsd <- cs %>% distinct(si, sp_logsd) %>% arrange(si) %>% pull(sp_logsd)
  rec_ci <- cs$ci; rec_si <- cs$si; rec_base <- cs$base_flux

  flux_samples <- matrix(0, nrow = N_ITER, ncol = ncells)
  for (it in 1:N_ITER) {
    shared_mult <- exp(rnorm(1, 0, SIGMA_BIAS) + rnorm(1, 0, SIGMA_PCONC) +
                       rnorm(1, 0, SIGMA_TURNOVER_ALLO))
    sp_indep <- exp(rnorm(nsp, 0, sp_logsd) + rnorm(nsp, 0, SIGMA_TURNOVER_RESID))
    agg <- tapply(rec_base * shared_mult * sp_indep[rec_si], rec_ci, sum)
    flux_samples[it, as.integer(names(agg))] <- agg
  }

  result <- data.frame(
    cell_id     = cells,
    flux_median = apply(flux_samples, 2, median),
    flux_lo     = apply(flux_samples, 2, quantile, 0.025),
    flux_hi     = apply(flux_samples, 2, quantile, 0.975)) %>%
    mutate(ci_ratio = flux_hi / pmax(flux_lo, 1e-9)) %>%
    left_join(read.csv(paste0(DIR_OUT, "grid_template.csv")), by = "cell_id")
  write.csv(result, "cell_flux_ci_P.csv", row.names = FALSE)

  global_samples <- rowSums(flux_samples) / 1e9
  message(sprintf("Block 5: global flux median = %.4f Mt P/yr, 95%% CI [%.4f, %.4f]",
                  median(global_samples),
                  quantile(global_samples, 0.025), quantile(global_samples, 0.975)))
  message("Block 5 done: cell_flux_ci_P.csv")
}


# ----- Block 6: Figures
################################################################################
# Produce the manuscript figures and their underlying data tables:
#   fig_carcass_P            global carcass-P supply field        (main Fig. 1)
#   fig_size_spectra_P       body-size spectra at four sites      (main Fig. 2)
#   fig_VMR                  spatial heterogeneity (VMR)          (main Fig. 3)
#   fig_standing_P           standing P in living mammals         (Extended Data)
#   fig_P_ci                 relative uncertainty of P flux       (Extended Data)
#   fig_Pweighted_bodymass   P-supply-weighted mean body mass     (supporting)
# All use the same NPP-weighted within-range allocation as Block 4.
#
# Input : species_flux_props_NP.csv, cell_species_long.csv, NPP, Greenspoon,
#         grid_template.csv, cell_flux_ci_P.csv (optional, for the CI map)
# Output: figure/*.png and figure/*.csv
################################################################################

run_block6_figures <- function() {
  suppressPackageStartupMessages({ library(terra); library(dplyr); library(viridis) })

  grid <- make_grid()
  grid_template <- read.csv(paste0(DIR_OUT, "grid_template.csv"))

  ## Map legacy palette names onto viridis options.
  pal_cols <- function(n, name, rev = FALSE) {
    opt <- switch(name, YlGnBu = "mako", Magma = "magma",
                  Viridis = "viridis", YlOrRd = "rocket", "viridis")
    viridis::viridis(n, option = opt, direction = if (rev) -1 else 1)
  }
  to_raster <- function(df, col) { r <- rast(grid); values(r) <- NA; r[df$cell_id] <- df[[col]]; r }

  ## Whole-body P (standing stock) uses the same stoichiometry as Block 4.
  conc_P_total <- function(W_kg) A_SKEL * W_kg^B_SKEL * P_BONE + W_kg * P_SOFT_FRAC

  ## Rebuild cell-by-species with the NPP-weighted allocation (matches Block 4).
  props <- read.csv("species_flux_props_NP.csv")
  gp <- read.csv(PATH_GREENSPOON); pop <- setNames(gp$estimated_population, gp$binomial)
  npp_by_cell <- npp_per_cell(grid)
  cs <- read.csv(paste0(DIR_OUT, "cell_species_long.csv")) %>%
    filter(binomial %in% names(pop), binomial %in% props$binomial) %>%
    left_join(npp_by_cell, by = "cell_id") %>%
    mutate(npp_filled = ifelse(is.na(npp), 0, npp), has_npp = !is.na(npp)) %>%
    group_by(binomial) %>%
    mutate(npp_sum = sum(npp_filled), n_cells = n(), valid_count = sum(has_npp),
           weight = if_else(valid_count > 0 & npp_sum > 0,
                            if_else(has_npp, npp_filled / npp_sum, 0), 1 / n_cells)) %>%
    ungroup() %>%
    mutate(ind = pop[binomial] * weight) %>%
    left_join(props, by = "binomial") %>%
    mutate(standing_P_kg = ind * conc_P_total(bm_kg),
           carcass_P_kg  = ind * turnover * P_supply_per_carcass)

  ## --- Standing P (Extended Data) ---
  cell_sp <- cs %>% group_by(cell_id) %>%
    summarise(standing_P_kg = sum(standing_P_kg, na.rm = TRUE), .groups = "drop") %>%
    left_join(grid_template, by = "cell_id") %>%
    mutate(area_km2 = cell_area_km2(lat), standing_P_per_km2 = standing_P_kg / area_km2)
  write.csv(cell_sp %>% select(cell_id, lon, lat, standing_P_kg, standing_P_per_km2),
            paste0(DIR_FIG, "fig_standing_P.csv"), row.names = FALSE)
  r <- to_raster(cell_sp, "standing_P_per_km2")
  v <- cell_sp$standing_P_per_km2[cell_sp$standing_P_per_km2 > 0]
  brk <- unique(as.numeric(c(0, quantile(v, c(.5,.7,.8,.88,.93,.97,.99,.997,1), na.rm = TRUE))))
  png(paste0(DIR_FIG, "fig_standing_P.png"), width = 1500, height = 680, res = 120)
  par(mar = c(3.2,3.2,2.6,7))
  plot(r, breaks = brk, col = pal_cols(length(brk)-1, "YlGnBu", rev = TRUE),
       main = "Standing phosphorus in living mammals (kg P/km2)")
  dev.off()

  ## --- Carcass-P supply field (main Fig. 1) ---
  cell_cp <- cs %>% group_by(cell_id) %>%
    summarise(carcass_P_kg = sum(carcass_P_kg, na.rm = TRUE), .groups = "drop") %>%
    left_join(grid_template, by = "cell_id") %>%
    mutate(area_km2 = cell_area_km2(lat), carcass_P_per_km2 = carcass_P_kg / area_km2)
  write.csv(cell_cp %>% select(cell_id, lon, lat, carcass_P_kg, carcass_P_per_km2),
            paste0(DIR_FIG, "fig_carcass_P.csv"), row.names = FALSE)
  r <- to_raster(cell_cp, "carcass_P_per_km2")
  v <- cell_cp$carcass_P_per_km2[cell_cp$carcass_P_per_km2 > 0]
  brk <- unique(as.numeric(c(0, quantile(v, c(.5,.7,.8,.88,.93,.97,.99,.997,1), na.rm = TRUE))))
  png(paste0(DIR_FIG, "fig_carcass_P.png"), width = 1500, height = 680, res = 120)
  par(mar = c(3.2,3.2,2.6,7))
  plot(r, breaks = brk, col = pal_cols(length(brk)-1, "YlGnBu", rev = TRUE),
       main = "Carcass-supplied phosphorus (kg P/km2/yr)")
  dev.off()

  ## --- VMR: spatial heterogeneity (main Fig. 3) ---
  ## Each 0.5-degree cell is a compound-Poisson landscape; carcasses fall as
  ## points, each depositing P_i (= P_supply_per_carcass). With carcass_P_kg =
  ## F_i * P_i, the variance-to-mean ratio is sum(F_i P_i^2) / sum(F_i P_i),
  ## independent of quadrat area.
  cell_vmr <- cs %>% filter(carcass_P_kg > 0) %>% group_by(cell_id) %>%
    summarise(sumFP  = sum(carcass_P_kg),                          # sum F_i P_i
              sumFP2 = sum(carcass_P_kg * P_supply_per_carcass),   # sum F_i P_i^2
              n_sp   = n(), .groups = "drop") %>%
    mutate(VMR = sumFP2 / pmax(sumFP, 1e-12)) %>%
    left_join(grid_template, by = "cell_id")
  write.csv(cell_vmr %>% select(cell_id, lon, lat, sumFP, sumFP2, n_sp, VMR),
            paste0(DIR_FIG, "fig_VMR.csv"), row.names = FALSE)
  r <- to_raster(cell_vmr, "VMR")
  v <- cell_vmr$VMR[is.finite(cell_vmr$VMR) & cell_vmr$VMR > 0]
  brk <- unique(as.numeric(c(0, quantile(v, c(.3,.5,.7,.8,.9,.95,.98,.995,1), na.rm = TRUE))))
  png(paste0(DIR_FIG, "fig_VMR.png"), width = 1500, height = 680, res = 120)
  par(mar = c(3.2,3.2,2.6,7))
  plot(r, breaks = brk, col = pal_cols(length(brk)-1, "Magma", rev = TRUE),
       main = "Spatial heterogeneity of carcass P supply: VMR = E[P^2]/E[P] (kg P)")
  dev.off()

  ## --- P-supply-weighted mean body mass (supporting) ---
  cell_wm <- cs %>% group_by(cell_id) %>%
    summarise(P_weighted_bodymass = sum(bm_kg * carcass_P_kg, na.rm = TRUE) /
                                     pmax(sum(carcass_P_kg, na.rm = TRUE), 1e-12),
              totP = sum(carcass_P_kg, na.rm = TRUE), .groups = "drop") %>%
    filter(totP > 0) %>% left_join(grid_template, by = "cell_id")
  write.csv(cell_wm %>% select(cell_id, lon, lat, P_weighted_bodymass),
            paste0(DIR_FIG, "fig_Pweighted_bodymass.csv"), row.names = FALSE)
  r <- to_raster(cell_wm, "P_weighted_bodymass")
  brk <- c(0,1,5,10,25,50,100,250,500,1000,2000,
           max(2001, max(cell_wm$P_weighted_bodymass, na.rm = TRUE)))
  png(paste0(DIR_FIG, "fig_Pweighted_bodymass.png"), width = 1500, height = 680, res = 120)
  par(mar = c(3.2,3.2,2.6,7))
  plot(r, breaks = brk, col = pal_cols(length(brk)-1, "Viridis", rev = FALSE),
       main = "P-supply-weighted mean body mass (kg)")
  dev.off()

  ## --- Body-size spectra at four representative sites (main Fig. 2) ---
  sites <- list(
    "Africa rainforest (Lope, Gabon)" = c(-0.2, 11.6),
    "Savanna (Serengeti)"             = c(-2.3, 34.8),
    "S. America rainforest (Amazon)"  = c(-4.0, -62.0),
    "Australia"                       = c(-25.0, 148.0))
  size_bins <- 10^seq(-2, 4, length.out = 25)
  size_ctr  <- sqrt(head(size_bins, -1) * tail(size_bins, -1))
  spectra_long <- list()
  png(paste0(DIR_FIG, "fig_size_spectra_P.png"), width = 1700, height = 520, res = 120)
  par(mfrow = c(1, 4), mar = c(4, 4, 3.5, 1))
  for (nm in names(sites)) {
    ll <- sites[[nm]]
    sub <- cs %>% filter(abs(lat - ll[1]) < 0.75, abs(lon - ll[2]) < 0.75)
    sub$bin <- cut(sub$bm_kg, breaks = size_bins, labels = size_ctr)
    g <- sub %>% group_by(bin) %>%
      summarise(P = sum(carcass_P_kg, na.rm = TRUE), .groups = "drop") %>%
      mutate(size = as.numeric(as.character(bin))) %>% filter(!is.na(size))
    spectra_long[[nm]] <- g %>% mutate(site = nm) %>% select(site, size, P)
    ymax <- max(g$P / 1e3, 0.01, na.rm = TRUE)
    plot(g$size, g$P / 1e3, type = "h", log = "x", lwd = 6, col = "#2166ac",
         xlim = c(0.01, 10000), ylim = c(0, ymax * 1.05),
         xlab = "Body mass (kg)", ylab = "Carcass P (10^3 kg P/yr per size class)",
         main = sprintf("%s\nmax %.0f kg", nm, ifelse(nrow(sub) > 0, max(sub$bm_kg), NA)))
    abline(v = 100, col = "red", lty = 3)  # megafauna threshold (100 kg)
  }
  dev.off()
  do.call(rbind, spectra_long) %>% write.csv(paste0(DIR_FIG, "fig_size_spectra_P.csv"), row.names = FALSE)

  ## --- Relative uncertainty map (Extended Data), if Block 5 has been run ---
  if (file.exists("cell_flux_ci_P.csv")) {
    ci <- read.csv("cell_flux_ci_P.csv")
    png(paste0(DIR_FIG, "fig_P_ci.png"), width = 1500, height = 680, res = 120)
    par(mar = c(3.2,3.2,2.6,7))
    plot(to_raster(ci, "ci_ratio"), col = pal_cols(50, "YlOrRd", rev = TRUE),
         main = "Relative uncertainty of carcass P flux (95% CI ratio)")
    dev.off()
  }
  message("Block 6 done: figures written to ", DIR_FIG)
}


# ----- Block 7: Reported values
################################################################################
# Recompute every value reported in the manuscript and print it with a [[TAG]]
# label (also saved to report_values.txt). Continent/biome quantities use the
# WWF realm/biome of each cell centre. The body-size centroid and VMR contrasts,
# and the global flux with its CI decomposition, are produced here.
#
# Sensitivity to the allocation rule is evaluated for the NPP-weighted rule and a
# uniform within-range rule.
#
# Input : species_flux_props_NP.csv, cell_species_long.csv, NPP, Greenspoon,
#         grid_template.csv, WWF Terrestrial Ecoregions
# Output: report_values.txt, centroid_summary.csv, vmr_summary.csv,
#         sensitivity_summary.csv
################################################################################

run_block7_report_values <- function() {
  suppressPackageStartupMessages({ library(terra); library(dplyr); library(tidyr) })

  grid <- make_grid()
  npp_by_cell <- npp_per_cell(grid)

  props <- read.csv("species_flux_props_NP.csv")
  gp <- read.csv(PATH_GREENSPOON)
  pop <- setNames(gp$estimated_population, gp$binomial)
  gp$log_sd <- (log(pmax(gp$estimated_population_97pt5, 1)) -
                log(pmax(gp$estimated_population_2pt5, 1))) / (2 * 1.96)
  pop_logsd <- setNames(gp$log_sd, gp$binomial)
  grid_template <- read.csv(paste0(DIR_OUT, "grid_template.csv"))

  LOG <- file("report_values.txt", open = "wt")
  say <- function(...) { msg <- sprintf(...); cat(msg, "\n"); cat(msg, "\n", file = LOG) }
  hdr <- function(t)  say("\n========== %s ==========", t)

  ## Rebuild cell-by-species under a chosen within-range allocation rule.
  build_weighted <- function(weight_col = c("npp", "uniform")) {
    weight_col <- match.arg(weight_col)
    read.csv(paste0(DIR_OUT, "cell_species_long.csv")) %>%
      filter(binomial %in% names(pop), binomial %in% props$binomial) %>%
      left_join(npp_by_cell, by = "cell_id") %>%
      mutate(npp_filled = ifelse(is.na(npp), 0, npp), has_npp = !is.na(npp)) %>%
      group_by(binomial) %>%
      mutate(npp_sum = sum(npp_filled), n_cells = n(), valid_count = sum(has_npp),
             weight = if (weight_col == "uniform") 1 / n_cells else
               if_else(valid_count > 0 & npp_sum > 0,
                       if_else(has_npp, npp_filled / npp_sum, 0), 1 / n_cells)) %>%
      ungroup() %>%
      mutate(ind = pop[binomial] * weight) %>%
      left_join(props, by = "binomial") %>%
      mutate(deaths = ind * turnover, carcass_P = deaths * P_supply_per_carcass, Fi = deaths)
  }
  cs <- build_weighted("npp")

  ## Attach WWF realm/biome to each cell centre (NULL if the shapefile is absent).
  norm_realm <- function(x) dplyr::case_when(
    x %in% c("AT","Afrotropic","Afrotropical","Afrotropics") ~ "Africa",
    x %in% c("NT","Neotropic","Neotropical","Neotropics")    ~ "SAmerica",
    x %in% c("AA","Australasia","Australasian","Oceania")    ~ "Australia",
    TRUE ~ NA_character_)
  norm_biome <- function(x) dplyr::case_when(
    x %in% c("1","Tropical & Subtropical Moist Broadleaf Forests","Tropical and Subtropical Moist Broadleaf Forests") ~ "TropicalMoistForest",
    x %in% c("2","Tropical & Subtropical Dry Broadleaf Forests","Tropical and Subtropical Dry Broadleaf Forests") ~ "TropicalDryForest",
    x %in% c("4","Temperate Broadleaf & Mixed Forests","Temperate Broadleaf and Mixed Forests") ~ "TemperateBroadleaf",
    x %in% c("7","Tropical & Subtropical Grasslands, Savannas & Shrublands","Tropical and Subtropical Grasslands, Savannas and Shrublands") ~ "TropicalGrasslandSavanna",
    x %in% c("8","Temperate Grasslands, Savannas & Shrublands","Temperate Grasslands, Savannas and Shrublands") ~ "TemperateGrassland",
    x %in% c("12","Mediterranean Forests, Woodlands & Scrub","Mediterranean Forests, Woodlands and Scrub") ~ "Mediterranean",
    x %in% c("13","Deserts & Xeric Shrublands","Deserts and Xeric Shrublands") ~ "Desert",
    TRUE ~ NA_character_)

  get_realm_biome <- function() {
    if (!file.exists(PATH_ECOREGION)) {
      say("WARNING: WWF ecoregion shapefile not found (%s).", PATH_ECOREGION)
      say("         realm/biome-dependent values are skipped.")
      return(NULL)
    }
    eco <- vect(PATH_ECOREGION)
    realm_col <- find_col(eco, c("REALM","WWF_REALM","WWF_REALM2","REALM_NAME","RealmName","REALM_1"))
    biome_col <- find_col(eco, c("BIOME","BIOME_NAME","WWF_MHTNAM","BIOME_NUM","BIOME_NUMBER"))
    if (is.na(realm_col) || is.na(biome_col)) stop("WWF realm/biome attribute columns not found.")
    pts <- grid_template %>% filter(!is.na(lon), !is.na(lat))
    cell_pts <- vect(data.frame(cell_id = pts$cell_id, lon = pts$lon, lat = pts$lat),
                     geom = c("lon", "lat"), crs = "EPSG:4326")
    ex <- terra::extract(eco[, c(realm_col, biome_col)], cell_pts)
    data.frame(cell_id = pts$cell_id,
               realm = as.character(ex[[realm_col]]),
               biome = as.character(ex[[biome_col]]))
  }
  rb <- get_realm_biome()

  ## --- [[GLOBAL_FLUX]] global flux + Monte Carlo CI ---
  hdr("[[GLOBAL_FLUX]]")
  global_point <- sum(cs$carcass_P, na.rm = TRUE) / 1e9
  say("[[GLOBAL_FLUX_MEDIAN]] %.4f Mt P/yr", global_point)

  set.seed(42)
  cs <- cs %>% mutate(sp_logsd = ifelse(is.na(pop_logsd[binomial]), 0.1, pop_logsd[binomial]))
  species <- unique(cs$binomial); nsp <- length(species)
  cs$si <- match(cs$binomial, species)
  sp_logsd_vec <- cs %>% distinct(si, sp_logsd) %>% arrange(si) %>% pull(sp_logsd)
  base <- cs$carcass_P
  mc_global <- function(shared_on = TRUE, indep_on = TRUE) {
    g <- numeric(N_ITER)
    for (it in 1:N_ITER) {
      sm <- if (shared_on) exp(rnorm(1,0,SIGMA_BIAS)+rnorm(1,0,SIGMA_PCONC)+rnorm(1,0,SIGMA_TURNOVER_ALLO)) else 1
      si <- if (indep_on) exp(rnorm(nsp,0,sp_logsd_vec)+rnorm(nsp,0,SIGMA_TURNOVER_RESID)) else rep(1, nsp)
      g[it] <- sum(base * sm * si[cs$si], na.rm = TRUE) / 1e9
    }
    g
  }
  g_all <- mc_global(TRUE, TRUE)
  ci_lo <- quantile(g_all, 0.025); ci_hi <- quantile(g_all, 0.975)
  say("[[GLOBAL_FLUX_CI]] 95%% CI [%.4f, %.4f] Mt P/yr (width %.1fx)", ci_lo, ci_hi, ci_hi / ci_lo)

  hdr("[[CI_DECOMP]]")
  fold <- function(x) as.numeric(quantile(x, 0.975) / quantile(x, 0.025))
  say("[[CI_DECOMP_SHARED]] shared only: %.1fx", fold(mc_global(TRUE, FALSE)))
  say("[[CI_DECOMP_INDEP]]  independent only: %.1fx", fold(mc_global(FALSE, TRUE)))

  ## --- realm/biome-dependent values ---
  if (!is.null(rb)) {
    csb <- cs %>% left_join(rb, by = "cell_id") %>%
      mutate(continent = norm_realm(realm), biome_lab = norm_biome(biome)) %>%
      filter(!is.na(continent))
    percell <- csb %>% group_by(cell_id, continent, biome_lab, lat) %>%
      summarise(sumFP = sum(carcass_P, na.rm = TRUE),
                sumFP2 = sum(carcass_P * P_supply_per_carcass, na.rm = TRUE),
                cell_P = sum(carcass_P, na.rm = TRUE), .groups = "drop") %>%
      mutate(VMR = sumFP2 / pmax(sumFP, 1e-12),
             area_km2 = cell_area_km2(lat), P_per_km2 = cell_P / area_km2)

    hdr("[[AREA_NORM]] kg P/km2/yr")
    an <- percell %>% group_by(continent) %>%
      summarise(area_norm = sum(cell_P, na.rm = TRUE) / sum(area_km2, na.rm = TRUE), .groups = "drop")
    for (k in seq_len(nrow(an))) say("[[AREA_NORM_%s]] %.3f", an$continent[k], an$area_norm[k])
    if (all(c("Africa","SAmerica") %in% an$continent))
      say("[[AREA_NORM_RATIO_AFR_SAM]] %.2f",
          an$area_norm[an$continent=="Africa"] / an$area_norm[an$continent=="SAmerica"])

    hdr("[[CENTROID]] supply-weighted mean body mass (kg)")
    centroid_tbl <- csb %>% group_by(continent, biome_lab) %>%
      summarise(centroid = sum(bm_kg*carcass_P, na.rm=TRUE)/pmax(sum(carcass_P, na.rm=TRUE), 1e-12),
                n_cells = n_distinct(cell_id), .groups = "drop") %>% filter(!is.na(biome_lab))
    centroid_all <- csb %>% group_by(continent) %>%
      summarise(centroid = sum(bm_kg*carcass_P, na.rm=TRUE)/pmax(sum(carcass_P, na.rm=TRUE), 1e-12),
                .groups = "drop") %>% mutate(biome_lab = "ALL")
    centroid_full <- bind_rows(centroid_tbl, centroid_all) %>% arrange(biome_lab, continent)
    write.csv(centroid_full, "centroid_summary.csv", row.names = FALSE)
    for (bm in unique(centroid_full$biome_lab)) {
      sub <- centroid_full %>% filter(biome_lab == bm)
      for (k in seq_len(nrow(sub)))
        say("[[CENTROID_%s_%s]] %.1f kg (n_cells=%s)", bm, sub$continent[k],
            sub$centroid[k], ifelse("n_cells" %in% names(sub), sub$n_cells[k], NA))
      if (all(c("Africa","SAmerica") %in% sub$continent))
        say("[[CENTROID_RATIO_%s_AFR_SAM]] %.1fx", bm,
            sub$centroid[sub$continent=="Africa"] / sub$centroid[sub$continent=="SAmerica"])
    }

    hdr("[[VMR]] kg P (per-cell median, max)")
    vmr_tbl <- percell %>% group_by(continent, biome_lab) %>%
      summarise(VMR_median = median(VMR, na.rm=TRUE), VMR_max = max(VMR, na.rm=TRUE),
                n_cells = n(), .groups = "drop") %>% filter(!is.na(biome_lab))
    vmr_all <- percell %>% group_by(continent) %>%
      summarise(VMR_median = median(VMR, na.rm=TRUE), VMR_max = max(VMR, na.rm=TRUE),
                n_cells = n(), .groups = "drop") %>% mutate(biome_lab = "ALL")
    vmr_full <- bind_rows(vmr_tbl, vmr_all) %>% arrange(biome_lab, continent)
    write.csv(vmr_full, "vmr_summary.csv", row.names = FALSE)
    for (bm in unique(vmr_full$biome_lab)) {
      sub <- vmr_full %>% filter(biome_lab == bm)
      for (k in seq_len(nrow(sub)))
        say("[[VMR_%s_%s]] median=%.4f max=%.3f (n_cells=%d)", bm, sub$continent[k],
            sub$VMR_median[k], sub$VMR_max[k], sub$n_cells[k])
      if (all(c("Africa","SAmerica") %in% sub$continent))
        say("[[VMR_RATIO_%s_AFR_SAM]] %.1fx", bm,
            sub$VMR_median[sub$continent=="Africa"] / sub$VMR_median[sub$continent=="SAmerica"])
    }
  } else {
    say("(realm/biome values skipped: shapefile not present)")
  }

  ## --- [[SLOPE_P]] body-size dependence of per-carcass P ---
  hdr("[[SLOPE_P]] P-bodymass log-log slope")
  fitP     <- lm(log10(P_supply_per_carcass) ~ log10(bm_kg), data = props[props$P_supply_per_carcass > 0, ])
  fitPbody <- lm(log10(P_per_carcass)        ~ log10(bm_kg), data = props[props$P_per_carcass > 0, ])
  say("[[SLOPE_P_SUPPLY]] %.3f  (supply = skeletal P)", coef(fitP)[2])
  say("[[SLOPE_P_BODY]] %.3f  (reference: whole-body P)", coef(fitPbody)[2])

  ## --- [[SENSITIVITY]] allocation-rule sensitivity (NPP vs uniform) ---
  hdr("[[SENSITIVITY]] centroid / VMR ratio (Africa / S. America)")
  if (!is.null(rb)) {
    summarise_contrast <- function(weight_col) {
      csx <- build_weighted(weight_col) %>% left_join(rb, by = "cell_id") %>%
        mutate(continent = norm_realm(realm)) %>% filter(continent %in% c("Africa","SAmerica"))
      cen <- csx %>% group_by(continent) %>%
        summarise(centroid = sum(bm_kg*carcass_P, na.rm=TRUE)/sum(carcass_P, na.rm=TRUE), .groups = "drop")
      pc <- csx %>% group_by(cell_id, continent) %>%
        summarise(VMR = sum(carcass_P*P_supply_per_carcass, na.rm=TRUE)/pmax(sum(carcass_P, na.rm=TRUE), 1e-12),
                  .groups = "drop") %>%
        group_by(continent) %>% summarise(VMR_med = median(VMR, na.rm = TRUE), .groups = "drop")
      data.frame(rule = weight_col,
                 centroid_ratio = cen$centroid[cen$continent=="Africa"] / cen$centroid[cen$continent=="SAmerica"],
                 vmr_ratio      = pc$VMR_med[pc$continent=="Africa"] / pc$VMR_med[pc$continent=="SAmerica"])
    }
    sens <- bind_rows(summarise_contrast("npp"), summarise_contrast("uniform"))
    write.csv(sens, "sensitivity_summary.csv", row.names = FALSE)
    for (k in seq_len(nrow(sens)))
      say("[[SENSITIVITY_%s]] centroid_ratio=%.1f vmr_ratio=%.1f",
          toupper(sens$rule[k]), sens$centroid_ratio[k], sens$vmr_ratio[k])
  } else {
    say("(sensitivity skipped: shapefile not present)")
  }

  close(LOG)
  message("Block 7 done: report_values.txt, centroid_summary.csv, vmr_summary.csv, sensitivity_summary.csv")
}


################################################################################
# RUN ALL
#   Uncomment to run the full pipeline top to bottom. Blocks 1-3 need the raw
#   inputs; if you start from the provided intermediates, run Blocks 4-7 only.
################################################################################

run_all <- function(from_intermediates = TRUE) {
  if (!from_intermediates) {
    run_block1_rasterize_ranges()
    run_block2a_body_mass()
    run_block2b_life_tables()
    run_block2c_extrapolate()
    run_block2d_turnover()
    run_block3_density_maps()
  }
  run_block4_flux_pipeline()
  run_block5_uncertainty()
  run_block6_figures()
  run_block7_report_values()
}

## Example:
# run_all(from_intermediates = TRUE)   # uses provided intermediate CSVs
# run_all(from_intermediates = FALSE)  # full rebuild from raw inputs
