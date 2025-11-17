# ==============================================================================
# Data Preparation Script for RSV Package
# ==============================================================================
# Main steps:
# 1. Load and merge location data (SHRUG)
# 2. Load experimental study data
# 3. Add SECC consumption/income data
# 4. Add nighttime lights data (VIIRS)
# 5. Add satellite imagery features (MOSAIKS)
# 6. Create spillover indicators
# 7. Generate experimental and observational samples
# ==============================================================================

rm(list = ls())

# ------------------------------------------------------------------------------
# Load Libraries
# ------------------------------------------------------------------------------
# Install redivis if needed
tryCatch({
  find.package("redivis")
  print("redivis is installed.")
}, error = function(e) {
  print("Installing redivis...")
  devtools::install_github("redivis/redivis-r", ref = "main")
})
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(haven)
  library(tidyr)
  library(sf)
  library(geosphere)
  library(igraph)
})
options(readr.show_col_types = FALSE)


# ==============================================================================
# Helper Function: Spillover Affected Units
# ==============================================================================
#' Identify spillover-affected units using bipartite matching
#'
#' @param D Treatment indicator (0, 1, or NA)
#' @param centroid_lat Latitude of unit centroids
#' @param centroid_lon Longitude of unit centroids
#' @param b_km Distance threshold in kilometers
#' @return Logical vector indicating spillover-affected units
spillover_affected_bkm <- function(D, centroid_lat, centroid_lon, b_km) {
  # 1) Preprocess D (NA -> 0) and coerce to {0,1}
  D[is.na(D)] <- 0
  D <- as.integer(D)
  if (!all(D %in% c(0, 1))) stop("D must take values in {0,1,NA}.")
  
  d <- data.frame(
    D    = D,
    centroid_lon  = centroid_lon,
    centroid_lat  = centroid_lat
  )
  n <- nrow(d)
  
  if (n <= 1) {
    spillover_affected <- rep(FALSE, n)
    return(spillover_affected)
  }
  
  # 2) Build geodesic neighbor structure (within b_km)
  pts <- st_as_sf(d, coords = c("centroid_lon", "centroid_lat"), crs = 4326, remove = FALSE)
  
  # For each i, indices j within b (includes i)
  nb_list <- st_is_within_distance(pts, pts, dist = b_km * 1000, sparse = TRUE)
  
  # 3) Conflict edges: connect i<->j if dist<=b AND D_i != D_j
  edges_i <- integer(0)
  edges_j <- integer(0)
  for (i in seq_len(n)) {
    js <- nb_list[[i]]
    if (length(js)) {
      js <- js[js > i]  # skip self; ensure i<j
      if (length(js)) {
        diffD <- (D[i] != D[js])
        if (any(diffD)) {
          pick <- js[diffD]
          edges_i <- c(edges_i, rep.int(i, length(pick)))
          edges_j <- c(edges_j, pick)
        }
      }
    }
  }
  
  # If no conflicts, keep everyone
  if (length(edges_i) == 0) {
    spillover_affected <- rep(FALSE, n)
    return(spillover_affected)
  }
  
  # 4) Bipartition by D
  U_idx <- which(D == 0)   # left part
  V_idx <- which(D == 1)   # right part
  types <- rep(NA, n)
  types[U_idx] <- TRUE
  types[V_idx] <- FALSE
  
  # 5) Graph with all vertices present; add conflict edges
  g <- make_empty_graph(n = n, directed = FALSE)
  g <- add_edges(g, as.vector(rbind(edges_i, edges_j)))
  V(g)$type <- types
  
  # 6) Maximum bipartite matching (robust across igraph versions)
  fn <- NULL
  if ("max_bipartite_match" %in% getNamespaceExports("igraph")) {
    fn <- max_bipartite_match
  } else if ("maximum_bipartite_matching" %in% getNamespaceExports("igraph")) {
    fn <- maximum_bipartite_matching
  } else if ("maximum.bipartite.matching" %in% getNamespaceExports("igraph")) {
    fn <- maximum.bipartite.matching
  } else stop("No bipartite matching function found in your igraph installation.")
  
  mm <- fn(g, types = V(g)$type)
  
  # 7) Unify to a full 'mate' vector of length n (0 = unmatched)
  mate <- rep(0, n)
  m <- mm$matching
  if (is.null(m)) stop("Matching result missing.")
  
  vn <- V(g)$name
  fill_sym <- function(idx, vals) {
    vals2 <- as.integer(vals); vals2[is.na(vals2)] <- 0
    mate[idx] <<- ifelse(vals2 > 0, vals2, 0)
    ok <- which(mate[idx] > 0)
    if (length(ok)) {
      u <- idx[ok]; v <- mate[idx][ok]
      mate[v] <<- u
    }
  }
  
  if (length(m) == n) {
    # igraph 2.x: one entry per vertex; NA = unmatched
    fill_sym(seq_len(n), m)
  } else if (length(m) == length(U_idx)) {
    # older igraph: vector for U side
    fill_sym(U_idx, m)
  } else if (length(m) == length(V_idx)) {
    # older igraph: vector for V side
    fill_sym(V_idx, m)
  } else if (!is.null(names(m))) {
    # map by vertex names
    idx <- if (is.null(vn)) as.integer(names(m)) else match(names(m), vn)
    fill_sym(idx, m)
  } else {
    stop("Unrecognized matching format from igraph; please update igraph.")
  }
  
  # Optional validation if available
  if ("is_matching" %in% getNamespaceExports("igraph")) {
    ok <- is_matching(g, ifelse(mate == 0, NA_integer_, mate), types = V(g)$type)
    if (!isTRUE(ok)) warning("Constructed 'mate' does not validate as a matching.")
  }
  
  # 8) Alternating-BFS to get a minimum vertex cover (Kőnig)
  # Free U-vertices:
  match_u_to_v <- rep(0, n); match_v_to_u <- rep(0, n)
  match_u_to_v[U_idx] <- mate[U_idx]
  match_v_to_u[V_idx] <- mate[V_idx]
  
  adj <- adjacent_vertices(g, V(g), mode = "all")
  
  visitedU <- rep(FALSE, n)
  visitedV <- rep(FALSE, n)
  
  queue <- U_idx[match_u_to_v[U_idx] == 0]  # unmatched on U side
  if (length(queue)) visitedU[queue] <- TRUE
  
  while (length(queue) > 0) {
    u <- queue[1]; queue <- queue[-1]
    
    # Traverse unmatched edges u -> v
    nbrs <- as.integer(adj[[u]])
    if (length(nbrs)) {
      for (v in nbrs) {
        if (!visitedV[v] && match_u_to_v[u] != v) {
          visitedV[v] <- TRUE
          # If v is matched to u2, traverse back via matched edge
          u2 <- match_v_to_u[v]
          if (u2 != 0 && !visitedU[u2]) {
            visitedU[u2] <- TRUE
            queue <- c(queue, u2)
          }
        }
      }
    }
  }
  
  # Minimum vertex cover: C = (U \ visitedU) ∪ (V ∩ visitedV)
  min_cover <- c(U_idx[!visitedU[U_idx]], V_idx[visitedV[V_idx]])
  
  # 9) Maximum independent set = complement of min vertex cover
  mis_idx <- setdiff(seq_len(n), min_cover)
  spillover_affected <- !(seq_len(n) %in% mis_idx)
  return(spillover_affected)
}

# ==============================================================================
# 1. Load and Clean Muralidharan et al. Study Data
# ==============================================================================

sample_smartcard_levels <- c(
  "1" = "Experimental: Treated (2010)",
  "2" = "Experimental: Untreated (2011)", 
  "3" = "Experimental: Untreated (2012)",
  "D" = "Observational (N/A)"
)

study_data <- read_dta(
  "data-raw/muralidharanetal_replication_package/balance-for-ap-mandal-comparison.dta",
  col_select = c(uniqueM, district_name, subdistrict_name, wave) 
) %>%
  rename(clusters = ) %>%
  mutate(
    district_name    = tolower(district_name),
    subdistrict_name = gsub("[()]", "", gsub(".", " ", tolower(subdistrict_name), fixed = TRUE)),
    `Sample (Smartcard)` = recode_factor(wave, !!!sample_smartcard_levels),
    D                = case_when(
      `Sample (Smartcard)` == "Experimental: Treated (2010)"              ~ 1, # Experimental: Treated
      `Sample (Smartcard)` %in% c("Experimental: Untreated (2012)", "Experimental: Untreated (2011)") ~ 0, # Experimental: Untreated
      TRUE                             ~ NA # Holdout (Observational)
    ),
    clusters = paste(subdistrict_name, district_name, sep = " ")
  ) %>% 
  select(-wave, -uniqueM) 

# ==============================================================================
# 2. Get Shrid IDs
# ==============================================================================

# Name mappings for consistency
district_name_map <- c("sri potti sriramulu nellore" = "nellore", 
                       "ysr kadapa" = "kadapa")
subdistrict_name_map <- c("sirpur town" = "sirpur t")

shrug_loc <- read.table(
  "https://dataverse.harvard.edu/api/access/datafile/10742739",
  sep = "\t", header = TRUE
) %>%
  select(shrid2, state_name, district_name, subdistrict_name) %>%
  filter(state_name == "andhra pradesh") %>%
  mutate(
    district_name    = recode(district_name, !!!district_name_map),
    subdistrict_name = recode(subdistrict_name, !!!subdistrict_name_map)
  ) %>%
  select(-state_name)


# ==============================================================================
# 3. Load and Merge SECC Data
# ==============================================================================

# Consumption data
secc_cons <- read.table(
  "https://dataverse.harvard.edu/api/access/datafile/10742743",
  sep = "\t", header = TRUE
) %>%
  select(shrid2, Ycons_raw = secc_cons_pc_rural)

# Income data
secc_income <- read_dta(
  "https://dataverse.harvard.edu/api/access/datafile/10742876",
  col_select = c(shrid2, tot_p, tot_f, inc_5k_plus_share, inc_10k_plus_share)
) %>%
  rename(
    Ylowinc_raw = inc_5k_plus_share, 
    Ymidinc_raw = inc_10k_plus_share
  ) %>%
  mutate(tot_f = if_else(tot_p == 0, 0, tot_f))

# Merge consumption and income
secc <- full_join(secc_cons, secc_income, by = "shrid2")


# ==============================================================================
# 4. Merge All Data Sources
# ==============================================================================

base_data <- shrug_loc%>%
  inner_join(study_data, by = c("district_name", "subdistrict_name")) %>%
  left_join(secc, by = "shrid2") %>%
  filter(tot_p >= 100) %>%
  mutate(
    Ylowinc = as.integer(Ylowinc_raw == 0),
    Ymidinc = as.integer(Ymidinc_raw == 0),
    Ycons   = as.integer(Ycons_raw <= quantile(Ycons_raw, 0.25, na.rm = TRUE))
  )

# Check consumption threshold
cat(sprintf(
  "25th percentile of consumption: %.2f\n", 
  quantile(base_data$Ycons_raw, 0.25, na.rm = TRUE)
))

# Clean up and create cluster identifiers
base_data <- base_data %>% 
  select(-ends_with("_raw"))


# ==============================================================================
# 5. Load VIIRS Nighttime Lights Data
# ==============================================================================

viirs_annual <- read_dta("https://dataverse.harvard.edu/api/access/datafile/10742856") %>%
  filter(
    year %in% 2012:2021,
    category == "median-masked"
  ) %>%
  rename_with(~ sub("viirs_annual", "luminosity", .x)) %>%
  select(-category) %>%
  pivot_wider(
    id_cols = shrid2,
    names_from = year,
    values_from = c(
      luminosity_min,
      luminosity_max,
      luminosity_mean,
      luminosity_sum,
      luminosity_num_cells
    ),
    names_sep = "."
  ) %>% 
  relocate(
    paste0("luminosity_min.", 2012:2021), 
    paste0("luminosity_max.", 2012:2021), 
    paste0("luminosity_mean.", 2012:2021), 
    paste0("luminosity_sum.", 2012:2021), 
    paste0("luminosity_num_cells.", 2012:2021),  
    .after = shrid2
  )

# ==============================================================================
# 6. Load Spatial Data, Compute Centroids, and get MOSAIKS Satellite Features
# ==============================================================================

shrid_ids <- unique(base_data$shrid2)

shrids <- st_read(
  "data-raw/shrug-shrid-poly-gpkg/shrid2_open.gpkg", 
  quiet = TRUE
) %>%
  filter(shrid2 %in% shrid_ids)

# Compute centroids
centroids <- st_centroid(st_make_valid(shrids))
coords <- st_coordinates(centroids)

centroid_coords <- shrids %>%
  mutate(
    centroid_lon = coords[, 1], 
    centroid_lat = coords[, 2]
  ) %>%
  as.data.frame() %>%
  select(shrid2, centroid_lat, centroid_lon) %>%
  arrange(shrid2)

# Query MOSAIKS Satellite Features
dataset <- redivis::redivis$organization("SDSS_data_repository")$dataset("mosaiks:8bqm")

# Build bounding box for query
buffer <- 0.05
lon_min <- round(min(centroid_coords$centroid_lon), 2) - buffer
lon_max <- round(max(centroid_coords$centroid_lon), 2) + buffer
lat_min <- round(min(centroid_coords$centroid_lat), 2) - buffer
lat_max <- round(max(centroid_coords$centroid_lat), 2) + buffer

condition_str <- sprintf(
  "(lon >= %s AND lon <= %s) AND (lat >= %s AND lat <= %s)",
  lon_min, lon_max, lat_min, lat_max
)

# Query MOSAIKS coordinates
query <- dataset$query(sprintf("
  SELECT lon, lat
  FROM mosaiks_2019_planet:ergr
  WHERE %s
", condition_str))

mosaik_coords <- query$to_tibble()

# Match each centroid to nearest MOSAIKS coordinate
chunk_size <- 100
chunks <- list()
num_rows <- nrow(centroid_coords)
num_chunks <- ceiling(num_rows / chunk_size)

for (i in seq_len(num_chunks)) {
  cat(sprintf("Matching chunk %d of %d\n", i, num_chunks))
  start_row <- (i - 1) * chunk_size + 1
  end_row <- min(i * chunk_size, num_rows)
  chunk <- centroid_coords[start_row:end_row, ]
  
  # Compute distance matrix
  dist_matrix <- distm(
    chunk %>% select(centroid_lon, centroid_lat) %>% as.matrix(),
    mosaik_coords %>% as.matrix(),
    fun = distHaversine
  )
  
  # Find nearest MOSAIKS coordinate
  min_idx <- apply(dist_matrix, 1, which.min)
  min_dist_km <- dist_matrix[cbind(1:nrow(dist_matrix), min_idx)] / 1000
  
  chunk <- chunk %>%
    mutate(
      mosaik_idx = min_idx,
      distance_km = min_dist_km,
      lat = mosaik_coords$lat[mosaik_idx],
      lon = mosaik_coords$lon[mosaik_idx]
    )
  
  chunks[[i]] <- chunk
}

coordinates <- bind_rows(chunks) %>% 
  filter(distance_km <= 1) %>%
  select(-mosaik_idx, -distance_km)
  
# Download feature vectors in chunks
chunk_size <- 1000
num_rows <- nrow(coordinates)
num_chunks <- ceiling(num_rows / chunk_size)
chunks <- list()

for (i in seq_len(num_chunks)) {
  cat(sprintf("Querying features for chunk %d of %d\n", i, num_chunks))
  start_row <- (i - 1) * chunk_size + 1
  end_row <- min(i * chunk_size, num_rows)
  chunk <- coordinates[start_row:end_row, ]
  
  # Build WHERE clause for this chunk
  where_clause <- apply(chunk, 1, function(row) {
    sprintf("(lon = %s AND lat = %s)", row["lon"], row["lat"])
  }) %>%
    paste(collapse = " OR ")
  
  # Query features
  query <- dataset$query(sprintf("
    SELECT * EXCEPT(shapeGroup, adm1_shapeID_geoBoundaries, adm2_shapeID_geoBoundaries)
    FROM mosaiks_2019_planet:ergr
    WHERE %s
  ", where_clause))
  
  features <- query$to_tibble()
  colnames(features)[grep("X_", colnames(features))] <- paste0("satellite_", 1:4000)
  
  chunks[[i]] <- left_join(chunk, features, by = c("lat", "lon"))
}

satellite_features <- bind_rows(chunks)

# ==============================================================================
# 7. Create Final Dataset with Spillover Indicators
# ==============================================================================

smartcard_data <- base_data %>%
  inner_join(select(satellite_features, shrid2, centroid_lat, centroid_lon), by = "shrid2") %>%
  mutate(spillover_20km = spillover_affected_bkm(D, centroid_lat, centroid_lon, b_km = 20)) %>%
  inner_join(viirs_annual, by = "shrid2") %>%
  inner_join(satellite_features, by = "shrid2") %>%
  select(shrid2, spillover_20km, tot_p, tot_f, `Sample (Smartcard)`, D, Ycons, Ylowinc, Ymidinc, clusters, 
         starts_with("luminosity_"), starts_with("satellite_"))

# Strip Stata attributes
smartcard_data <- zap_labels(smartcard_data)
smartcard_data <- zap_formats(smartcard_data) 
for (i in seq_along(smartcard_data)) {
  attr(smartcard_data[[i]], "label") <- NULL
  attr(smartcard_data[[i]], "format.stata") <- NULL
}

# ==============================================================================
# 8. Save Datasets
# ==============================================================================
remote_vars_p1 <- smartcard_data %>% select(shrid2, starts_with("luminosity"), paste0("satellite_", 1:2000))
remote_vars_p2 <- smartcard_data %>% select(shrid2, paste0("satellite_", 2001:4000))
save(remote_vars_p1, file = "data/remote_vars_p1.rda", compress = "xz", compression_level = 9)
save(remote_vars_p2, file = "data/remote_vars_p2.rda", compress = "xz", compression_level = 9)

smartcard_data <- smartcard_data %>% select(-starts_with("luminosity"), -starts_with("satellite_"))
usethis::use_data(smartcard_data, overwrite = TRUE, compress = "xz")

cat("\nData preparation complete!\n")