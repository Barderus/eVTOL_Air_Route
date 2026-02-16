import geopandas as gpd
import pandas as pd

# Load shapefile and population file
gdf = gpd.read_file("../population/tl_2025_17_bg/tl_2025_17_bg.shp")
density_df = pd.read_excel("../population/population_data.xlsx")

# Copy to avoid modifying original
pop_df_clean = density_df.copy()

# Remove header-like row
pop_df_clean = pop_df_clean[pop_df_clean["GEO_ID"] != "Geography"].copy()

# Extract GEOID from GEO_ID
pop_df_clean["GEOID"] = pop_df_clean["GEO_ID"].str.split("US").str[-1]

# Convert population to numeric
pop_df_clean["population"] = pd.to_numeric(
    pop_df_clean["B01003_001E"], errors="coerce"
)

# Keep only what we need
pop_df_clean = pop_df_clean[["GEOID", "population"]]
print(pop_df_clean.shape)

# Merge both dataframe
gdf_merged = gdf.merge(pop_df_clean, on="GEOID", how="left")
print(gdf_merged["population"].isna().sum())

# Calculation of population density
gdf_merged["density_p_km2"] = (
    gdf_merged["population"] / (gdf_merged["ALAND"] / 1_000_000)
)
print(gdf_merged["density_p_km2"].describe())

# Let's get density risk now by using the quantiles
gdf_merged_crs = gdf_merged.to_crs("EPSG:3857")
q_low, q_med = gdf_merged_crs["density_p_km2"].quantile([0.70, 0.90])

gdf_merged_crs["density_risk"] = "low"
gdf_merged_crs.loc[gdf_merged_crs["density_p_km2"] >= q_low, "density_risk"] = "medium" # between 25 and 50 quantile
gdf_merged_crs.loc[gdf_merged_crs["density_p_km2"] >= q_med, "density_risk"] = "high" # > 50 quantile
print(q_low, q_med)
print(gdf_merged_crs["density_risk"].value_counts())

# Save as GeoJSON for Leaflet
out_geojson = "il_blockgroups_population_density.geojson"
gdf_merged_crs.to_file(out_geojson, driver="GeoJSON")
print("Saved:", out_geojson)

