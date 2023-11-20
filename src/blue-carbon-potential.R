# Modelling Blue carbon potential in South Australia

# install-packages
required_packages <- c("tidyverse", "terra", "smoothr", "pivottabler", "sf")
for (package in required_packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    require(package, character.only = TRUE)
  }
}

# set-temp-path
# Location to save temporary files
temp_dir <- "bc-potential-intermediary-data"
terra::terraOptions(tempdir = temp_dir)

# Calculating flooding extent and sea level rise To identify the maximum
# potential extent of potential future coastal wetlands, we extracted areas that
# fall between the 2020: MHWS tide level + 0 m SLR and 2100 : MHWS tide level +
# 1m SLR + 1% AEP storm surge level and wave set-up (MHWS 2020 minus ARI 2100 ->
# 'Full extent')

# Additional scenarios were also modelled:
# 2020MHWS: MHWS tide level + 0m SLR
# 2020ARI:  MHWS tide level + 0m SLR+1% AEP storm surge level and wave set-up
# 2050MHWS: MHWS tide level + 0.3m SLR
# 2050ARI:  MHWS tide level + 0.3m SLR+1% AEP storm surge level and wave set-up
# 2100MHWS: MHWS tide level + 1m SLR
# 2100ARI:  MHWS tide level + 1m SLR+1% AEP storm surge level and wave set-up
#  2020MHWS minus  2020ARI
#  2050MHWS minus  2050ARI
#  2100MHWS minus  2100ARI

# define-flood-extent-function
WriteFloodExtentVector <- function(
    flood_data_path,
    area_code_col, min_tide_col, max_tide_col,
    area_code,
    dem_path,
    out_dir, out_suffix) {
  flood_data <- read.csv(flood_data_path)
  flood_data_sub <- flood_data %>%
    dplyr::filter(.data[[area_code_col]] == area_code) %>%
    dplyr::summarise(
      min_tide = min(.data[[min_tide_col]], na.rm = TRUE),
      max_tide = max(.data[[max_tide_col]], na.rm = TRUE)
    ) %>%
    dplyr::select(min_tide, max_tide)

  flood_data_sub$reclass_value <- 1
  dem <- terra::rast(dem_path)
  dem_reclassed <- terra::classify(dem, flood_data_sub, others = NA)
  dem_reclassed <- terra::classify(dem_reclassed, cbind(0, NA))
  dem_reclassed_polygons <- terra::as.polygons(dem_reclassed)
  out_path <- file.path(
    out_dir, paste0("flood_extent_", area_code, "_", out_suffix, ".shp")
  )
  terra::writeVector(dem_reclassed_polygons, out_path, overwrite = TRUE)
}

# run-flood-extents
flood_data_path <- "flood_data.csv" # file w Min_Tides, Max_Tides, Area_Code
dem_dir <- "digital_elevation_models" # directory containing DEMs by area code
out_dir <- "flood_extent_output" # directory for saving flood extent shapefiles

# Load data
flood_data <- read.csv(flood_data_path)
area_codes <- unique(flood_data$Area_Code)
dem_paths <- list.files(dem_path, full.names = TRUE)

# Order DEM list by area code (assuming area_code is in the file name of DEM)
dem_paths[sapply(area_codes, function(x) {
  grep(x, dem_paths)
})]

# Set up df with parameter combinations
flood_cells <- data.frame(
  flood_data_path = flood_data_path,
  area_code_col = "Area_Code",
  area_code = area_codes,
  dem_path = dem_paths,
  out_dir = out_dir
)

min_tide_col <- list("X2020MHWS", "X2020MHWS", "X2050MHWS", "X2100MHWS")
max_tide_col <- list("X2100ARI", "X2020ARI", "X2050ARI", "X2100ARI")
out_suffix <- list("full_extent", "2020", "2050", "2100")
scenarios <- cbind(min_tide_col, max_tide_col, out_suffix)

# Create DF of all combo's of scenarios and flood cells
parameters <- tidyr::expand_grid(flood_cells, scenarios)

col_order <- c(
  "flood_data_path",
  "area_code_col",
  "min_tide_col",
  "max_tide_col",
  "area_code",
  "dem_path",
  "out_dir",
  "out_suffix"
)
parameters <- parameters[col_order]

apply(parameters, 1, WriteFloodExtentVector)

# Clean up
do.call(file.remove, list(list.files(temp_dir, full.names = TRUE)))

# Running veg classification and carbon potential To estimate and compare the
# potential carbon capture and storage potential project sites across the South
# Australian coastline, we compared the extent of existing coastal wetland
# vegetation (pre-intervention/present day) to the predicted coastal wetland
# vegetation (based on land elevation) for post-intervention (blockage removed –
# and sea level rise) scenarios.

# Pre-intervention scenario layers are constructed in ArcGIS with the following
# steps
# - Relative carbon values applied to each habitat type in the feature layer for
#   coastal wetland vegetation
# - Feature layer converted to raster format, where raster value equals applied 
#   relative carbon values
# - 'Extract by Mask' tool used to extract the above raster grid into each flood
#   cell (with unique area code)
# - Output: 'Pre-intervention' raster layer created for each 
#   flood cell/area code

# Post-intervention scenarios were constructed in R Studio (below), by
# reclassifying DEMs for each flood cell, with reclassification values relating
# to specific elevations bands for key wetland habitat types (veg
# classification), and constructed for years 2020, 2050 and 2100 by accounting
# for sea level rise.

# Carbon potential (below) was calculated by subtracting the 'pre-intervention
# scenario' from 'post-intervention scenario' raster grid to estimate the
# difference in carbon capture and storage over time

# define-carbon-potential-function
CalculateCarbonPotential <- function(
    flood_data_path,
    pre_path,
    area_code_col,
    area_code,
    year_col,
    low_col,
    high_col,
    reclass_col,
    dem_path,
    out_dir) {
  flood_data <- read.csv(flood_data_path)
  pre_data <- read.csv(pre_path)

  for (year in c("2020", "2050", "2050")) {
    flood_data_sub <- flood_data %>%
      dplyr::filter((.data[[area_code_col]] == area_code) &
        (.data[[year_col]] == year)) %>%
      dplyr::select(low_col, high_col, reclass_col)

    dem <- terra::rast(dem_path)
    dem_reclassed <-
      terra::classify(dem, flood_data_sub, others = NA)
    dem_reclassed <- terra::classify(dem_reclassed, cbind(NA, 0))
    pre_resampled <- terra::resample(pre_data, dem_reclassed)

    potential <- dem_reclassified - pre_resampled
    potential_poly <- terra::as.polygons(potential)

    tif_path <- file.path(out_dir, "raster")
    shp_path <- file.path(out_dir, "vector")
    if (!file.exists(tif_path)) {
      dir.create(tif_path, recursive = TRUE)
    }
    if (!file.exists(shp_path)) {
      dir.create(shp_path, recursive = TRUE)
    }
    file_name <- paste0("c_potential_", area_code, "_by_", year)
    tif_full <- file.path(tif_path, paste0(filename, ".tif"))
    shp_full <- file.path(shp_path, paste0(filename, ".shp"))

    terra::writeVector(potential, tif_full, overwrite = TRUE)
    terra::writeVector(potential_poly, shp_full, overwrite = TRUE)
  }
}

# run-carbon-potential
flood_data_path <- "flood_data.csv" # file containing LOW, HIGH, RECL, Area_Code
dem_dir <- "digital_elevation_models" # directory containing DEMs by area code
pre_path <- "veg_data"
out_dir <- "carbon_potential" # directory for saving output

# Load data
flood_data <- read.csv(flood_data_path)
area_codes <- unique(flood_data$Area_Code)
dem_paths <- list.files(dem_path, full.names = TRUE)
pre_paths <- list.files(pre_path, full.names = TRUE)
# Order DEM & pre list by area code (assuming area_code is in file name of DEM)
dem_paths[sapply(area_codes, function(x) {
  grep(x, dem_paths)
})]
pre_paths[sapply(area_codes, function(x) {
  grep(x, pre_paths)
})]

# Set up df with parameter combinations
parameters <- data.frame(
  flood_data_path = flood_data_path,
  pre_path = pre_path,
  area_code_col = "Area_Code",
  area_code = area_codes,
  year_col <- "Year",
  low_col <- "LOW",
  high_col <- "HIGH",
  reclass_col <- "RECL",
  dem_path = dem_paths,
  out_dir = out_dir,
)

apply(parameters, 1, CalculateCarbonPotential)

# Clean up
do.call(file.remove, list(list.files(temp_dir, full.names = TRUE)))
