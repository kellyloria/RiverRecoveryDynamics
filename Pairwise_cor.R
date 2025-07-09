# Step 1: Extract site metadata and time series
site_coords <- GPP_wide_day[, c("Site_ID", "Lat", "Lon")]
gpp_matrix <- as.matrix(GPP_wide_day[, 4:ncol(GPP_wide_day)])  # Time series (rows = sites, cols = days)

# Step 2: Compute pairwise correlations between rows (sites)
# Use pairwise.complete.obs to handle NAs
cor_matrix <- cor(t(gpp_matrix), use = "pairwise.complete.obs")  # transpose so rows = time points

# Step 3: Create a dataframe of unique site pairs and their correlation
site_ids <- site_coords$Site_ID
pairwise_corrs <- as.data.frame(as.table(cor_matrix))  # Convert matrix to long format
names(pairwise_corrs) <- c("Site1", "Site2", "GPP_corr")

# Step 4: Filter to keep only unique, non-self pairs
pairwise_corrs <- pairwise_corrs %>%
  dplyr::filter(Site1 != Site2) %>%
  dplyr::mutate(pair = paste(pmin(Site1, Site2), pmax(Site1, Site2), sep = "_")) %>%
  dplyr::distinct(pair, .keep_all = TRUE)

# Step 5 (optional): Add spatial distances between site pairs
library(geosphere)

# Create a site lookup table for coordinates
site_lookup <- site_coords %>%
  dplyr::select(Site_ID, Lat, Lon)

# Merge lat/lon for both sites in the pair
pairwise_corrs <- pairwise_corrs %>%
  left_join(site_lookup, by = c("Site1" = "Site_ID")) %>%
  rename(Lat1 = Lat, Lon1 = Lon) %>%
  left_join(site_lookup, by = c("Site2" = "Site_ID")) %>%
  rename(Lat2 = Lat, Lon2 = Lon)

# Compute distances (in kilometers)
pairwise_corrs$dist_km <- distHaversine(
  matrix(c(pairwise_corrs$Lon1, pairwise_corrs$Lat1), ncol = 2),
  matrix(c(pairwise_corrs$Lon2, pairwise_corrs$Lat2), ncol = 2)
) / 1000  # convert meters to km

# Final result
head(pairwise_corrs)
