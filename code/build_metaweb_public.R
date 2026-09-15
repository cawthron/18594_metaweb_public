

library(here)
library(cheddar)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(igraph)
library(cowplot)
library(furrr)
library(marmap)
library(ggsn)
library(ggrepel)
library(tidyverse)
library(progressr)
library(rglobi)

plan(multisession, workers = 4)
dir.create("results", showWarnings = FALSE)
source(here("code/functions/metaweb_functions.R"))
source(here("code/functions/map_functions.R"))
source(here("code/functions/resource_plots.R"))
source(here("code/functions/localfoodwebs.R"))
source(here("code/functions/nullmodels.R"))

# =============================================================================
# Create Map of Lakes
# =============================================================================


sd.plankton <- read.csv("data/SupplementaryTable1.csv", h=T)

sd.plankton$Latitude <- as.numeric(sd.plankton$Latitude)
sd.plankton$Longitude <- as.numeric(sd.plankton$Longitude)

sd.plankton$lat <- sd.plankton$Latitude
sd.plankton$long <- sd.plankton$Longitude


getNOAA.bathy(lon1=165.5,lon2=179,lat1=-48,lat2=-34, resolution=1, keep = T) -> NZ

getNOAA.bathy(lon1=168,lon2=172.2,lat1=-45,lat2=-42.3, resolution=1, keep = T) -> West.Coast


West.Coast.p2 <- plot_west_coast_map(
  sd.plankton = sd.plankton,
  nz_bathy = NZ,
  wc_bathy = West.Coast
)

# =============================================================================
# Build Edge List
# =============================================================================

master_taxa <- read.csv(here("data/master_taxa.csv"))

#Building the edge list only using the highest confidence list consumers.

edge_list_high  <- build_edge_list(master_taxa = master_taxa,
                                        all_consumers =  FALSE,
                                        test = FALSE)

edge_list_high <- edge_list_high$edge_list

edge_list_high |>
  group_by(consumer) |>
  summarise("no_resources"= n()) |>
  pull(no_resources) |>
  hist()

#Building the edge list only using the all consumers.

edge_list_complete  <- build_edge_list(master_taxa = master_taxa,
                                        all_consumers =  TRUE,
                                        test = FALSE)

edge_list_complete <- edge_list_complete$edge_list

edge_list_complete |>
  group_by(consumer) |>
  summarise("no_resources"= n()) |>
  pull(no_resources) |>
  hist()


result <- plot_resource_histograms(
  edge_list_high,
  edge_list_complete,
  save_plot = TRUE,
  save_table = TRUE
)


# =============================================================================
# Create Adjacency Matrix
# =============================================================================

#Build a adjacency matrix for the high confidence list 

metaweb_high <- build_metaweb(edge_list = edge_list_high,
                               master_taxa = master_taxa,
                               plot=TRUE)

adjacency_matrix_high <- metaweb_high$adjacency_matrix
community_hp <- metaweb_high$community


write.csv(adjacency_matrix_high, file=here("results/metaweb_binary_high.csv"))

#Build a adjacency matrix for the all consumer list 

metaweb_complete <- build_metaweb(edge_list = edge_list_complete,
                               master_taxa = master_taxa,
                               plot=TRUE)

adjacency_matrix_complete <- metaweb_complete$adjacency_matrix
community_hp <- metaweb_complete$community


write.csv(adjacency_matrix_complete, file=here("results/metaweb_binary_complete.csv"))


# =============================================================================
# LOAD Adjacency Data
# =============================================================================

adjacency_high     <- read.csv(here("results/metaweb_binary_high.csv"),
                               row.names = 1, check.names = FALSE)
adjacency_complete <- read.csv(here("results/metaweb_binary_complete.csv"),
                               row.names = 1, check.names = FALSE)

# =============================================================================
# LOAD eDNA Lake Data
# =============================================================================

pa_table <- read.csv("data/pa_table_final.csv", row.names = 1)
# row.names = 1 assumes first column is species names

# Transpose so rows = lakes, cols = species
pa_table_t <- t(pa_table) |> as.data.frame()

# Check dimensions
dim(pa_table_t)        # should be 15 lakes × n species
rownames(pa_table_t)   # should be your lake names
colnames(pa_table_t)


colnames(pa_table_t) <- colnames(pa_table_t) |>
  str_replace_all(" ", "_") |>
  str_trim()

lake_names <- rownames(pa_table_t)
cat("Lakes:", length(lake_names), "\n")
cat("Species in PA table:", ncol(pa_table_t), "\n")

# =============================================================================
# BUILD IGRAPH OBJECTS
# =============================================================================

build_igraph <- function(adj_matrix) {
  graph_from_adjacency_matrix(
    as.matrix(adj_matrix),
    mode     = "directed",
    diag     = FALSE,
    weighted = NULL
  )
}

metaweb_high     <- build_igraph(adjacency_high)
metaweb_complete <- build_igraph(adjacency_complete)

V(metaweb_high)$name     <- V(metaweb_high)$name     |> str_replace_all(" ", "_") |> str_trim()
V(metaweb_complete)$name <- V(metaweb_complete)$name |> str_replace_all(" ", "_") |> str_trim()

cat("Metaweb (high) nodes:",     vcount(metaweb_high),
    "| edges:", ecount(metaweb_high), "\n")
cat("Metaweb (complete) nodes:", vcount(metaweb_complete),
    "| edges:", ecount(metaweb_complete), "\n")

# Sanity check — these should differ
stopifnot(
  "Metaweb edge counts are identical — check that both CSVs were saved correctly" =
    ecount(metaweb_high) != ecount(metaweb_complete)
)



metaweb_high_cheddar <- igraph_to_cheddar(metaweb_high, "HCL")

TopLevelNodes(metaweb_high_cheddar)

BasalNodes(metaweb_high_cheddar)


metaweb_complete_cheddar <- igraph_to_cheddar(metaweb_complete, "AC")

TopLevelNodes(metaweb_complete_cheddar)

BasalNodes(metaweb_complete_cheddar)

# =============================================================================
# SPECIES MATCHING DIAGNOSTICS
# =============================================================================

match_high     <- diagnose_matching(pa_table_t, metaweb_high,     "High Potential")
match_complete <- diagnose_matching(pa_table_t, metaweb_complete, "Complete")


# =============================================================================
# COMPUTE METRICS FOR METAWEBS
# =============================================================================

metrics_metaweb_high <- compute_network_metrics(
  metaweb_high,
  community_name = "Metaweb_High_Potential"
)

metrics_metaweb_complete <- compute_network_metrics(
  metaweb_complete,
  community_name = "Metaweb_Complete"
)

# Combine into one summary table
metrics_metawebs <- bind_rows(
  data.frame(metaweb = "High Potential", metrics_metaweb_high),
  data.frame(metaweb = "Complete",       metrics_metaweb_complete)
)

# Optional: save
write.csv(
  metrics_metawebs,
  here("results/metaweb_network_metrics.csv"),
  row.names = FALSE
)

# =============================================================================
# VERIFY CHEDDAR EDGE DIRECTION ON A TEST LAKE
# =============================================================================
# Run this block once and check that known top predators appear in TopLevelNodes()
# If basal/top are swapped, flip "consumer"/"resource" in igraph_to_cheddar()

cat("\n--- Cheddar edge direction check ---\n")
test_lake   <- lake_names[1]
test_g      <- build_local_webs_single <- {
  mw_spp  <- V(metaweb_complete)$name
  present <- colnames(pa_table_t)[pa_table_t[test_lake, ] == 1]
  spp     <- intersect(present, mw_spp)
  induced_subgraph(metaweb_complete, vids = spp)
}
test_ch <- igraph_to_cheddar(test_g, test_lake)
cat("Basal nodes (should be primary producers):\n")
print(head(BasalNodes(test_ch), 10))
cat("Top nodes (should be apex predators):\n")
print(head(TopLevelNodes(test_ch), 10))
cat("Trophic levels (basal ~1, top = highest):\n")
print(sort(PreyAveragedTrophicLevel(test_ch)) |> tail(10))




# =============================================================================
# BUILD LOCAL WEBS (INDUCED SUBGRAPHS)
# =============================================================================

pa_table_t <- pa_table_t %>%
  rownames_to_column("sample_id") %>%
  filter(sample_id %in% c("Lake_Alexandrina",
                              "Lake_Douglas",
                              "Lake_Hartley_Tarn", 
                              "Lake_Hochstetter", 
                              "Lake_Haupiri",
                              "Lake_Ianthe", 
                              "Lake_Kangaroo", 
                              "Lake_Lady",
                              "Lake_Mahinapua", 
                              "Lake_Mapourika",
                              "Lake_Mcgregor", 
                              "Lake_Moeraki",
                              "Lake_Paringa",
                              "Lake_Poerua", 
                              "Stoney_Tarn"))  %>%
  column_to_rownames("sample_id")


local_webs_high     <- build_local_webs(pa_table_t, metaweb_high,     "High Potential")
local_webs_complete <- build_local_webs(pa_table_t, metaweb_complete, "Complete")

# =============================================================================
# Local Observed Metrics
# =============================================================================

#For High confidence list (HCL)
observed_high <- map_dfr(
  names(local_webs_high),
  function(lake) {
    compute_network_metrics(local_webs_high[[lake]], lake) |>
      mutate(lake = lake)
  }
) |> mutate(metaweb = "HCL")


observed_high_summary <- observed_high %>%
  summarise(across(where(is.numeric), list(
    mean = ~mean(.x, na.rm = TRUE),
    min  = ~min(.x, na.rm = TRUE),
    max  = ~max(.x, na.rm = TRUE),
    sd   = ~sd(.x, na.rm = TRUE)
  ))) |>
  pivot_longer(everything(),
               names_to  = c("metric", ".value"),
               names_sep = "_(?=[^_]+$)")

write.csv(observed_high_summary, "results/observed_high_summary.csv", row.names = F)

#For All consumers
observed_complete <- map_dfr(
  names(local_webs_complete),
  function(lake) {
    compute_network_metrics(local_webs_complete[[lake]], lake) |>
      mutate(lake = lake)
  }
) |> mutate(metaweb = "AC")



observed_complete_summary <- observed_complete %>%
  summarise(across(where(is.numeric), list(
    mean = ~mean(.x, na.rm = TRUE),
    min  = ~min(.x, na.rm = TRUE),
    max  = ~max(.x, na.rm = TRUE),
    sd   = ~sd(.x, na.rm = TRUE)
  ))) |>
  pivot_longer(everything(),
               names_to  = c("metric", ".value"),
               names_sep = "_(?=[^_]+$)")

write.csv(observed_complete_summary, "results/observed_complete_summary.csv", row.names = F)

observed_metrics <- bind_rows(observed_high, observed_complete)

write.csv(observed_metrics,
          here("results/observed_network_metrics.csv"),
          row.names = FALSE)

cat("\nObserved metrics summary:\n")
print(observed_metrics)

# =============================================================================
# NULL MODEL — RANDOM SUBWEBS FROM METAWEB
# =============================================================================

# Non-redundant metrics for SES
metrics_cols <- c(
  "connectance",
 # "clustering",
  "modularity",
  "prop_basal",
  "prop_top",
 # "prop_int",
#  "prop_omnivory",
  "trophic_level_mean",
  "trophic_level_max",
  "max_chain_length",
  "mean_chain_length"
)

null_high     <- run_null_model(metaweb_high,     local_webs_high,     "High Potential", n_null = 999)
null_complete <- run_null_model(metaweb_complete, local_webs_complete, "Complete",       n_null = 999)

saveRDS(null_high,     here("results/null_distributions_high.rds"))
saveRDS(null_complete, here("results/null_distributions_complete.rds"))

# =============================================================================
# STANDARDISED EFFECT SIZES (SES)
# FDR correction applied WITHIN each metric across lakes (not globally)
# =============================================================================

ses_high <- compute_ses(observed_high, null_high, "high", metrics_cols)

ses_high_summary <- data.frame(
  metric       = c("connectance", "modularity", "prop_basal", 
                   "prop_top", "trophic_level_mean", "trophic_level_max",
                   "mean_chain_length", "max_chain_length"),
  N_lakes      = nrow(ses_high)
) |>
  rowwise() |>
  mutate(
    col          = paste0("SES_", metric),
    Positive     = sum(ses_high[[col]] > 0,     na.rm = TRUE),
    Positive_Sig = sum(ses_high[[col]] > 1.96,  na.rm = TRUE),
    Negative     = sum(ses_high[[col]] < 0,     na.rm = TRUE),
    Negative_Sig = sum(ses_high[[col]] < -1.96, na.rm = TRUE)
  ) |>
  ungroup() |>
  select(metric, Positive, Positive_Sig, Negative, Negative_Sig, N_lakes)

write.csv(ses_high_summary, "results/ses_high.summary.csv", row.names = F)

ses_complete <- compute_ses(observed_complete, null_complete, "complete",     metrics_cols)

ses_complete_summary <- data.frame(
  metric       = c("connectance", "modularity", "prop_basal", 
                   "prop_top", "trophic_level_mean", "trophic_level_max",
                   "mean_chain_length", "max_chain_length"),
  N_lakes      = nrow(ses_high)
) |>
  rowwise() |>
  mutate(
    col          = paste0("SES_", metric),
    Positive     = sum(ses_complete[[col]] > 0,     na.rm = TRUE),
    Positive_Sig = sum(ses_complete[[col]] > 1.96,  na.rm = TRUE),
    Negative     = sum(ses_complete[[col]] < 0,     na.rm = TRUE),
    Negative_Sig = sum(ses_complete[[col]] < -1.96, na.rm = TRUE)
  ) |>
  ungroup() |>
  select(metric, Positive, Positive_Sig, Negative, Negative_Sig, N_lakes)

write.csv(ses_complete_summary, "results/ses_complete.summary.csv", row.names = F)

ses_all      <- bind_rows(ses_high, ses_complete)

ses_all <- ses_all %>% filter(lake == "Lake_Alexandrina" |
                              lake == "Lake_Douglas" |
                              lake == "Lake_Hartley_Tarn" | 
                              lake == "Lake_Haupiri" | 
                              lake == "Lake_Hochstetter" | 
                              lake == "Lake_Ianthe" | 
                              lake == "Lake_Kangaroo" | 
                              lake == "Lake_Lady" | 
                              lake == "Lake_Mahinapua" | 
                              lake == "Lake_Mapourika" |
                              lake == "Lake_Mcgregor" | 
                              lake == "Lake_Moeraki" |
                              lake == "Lake_Paringa" | 
                              lake == "Lake_Poerua" | 
                              lake == "Stoney_Tarn")

# FDR correction WITHIN each metric across lakes (correct approach)
ses_long_clean <- ses_all |>
  select(lake, metaweb, S,
         all_of(paste0("SES_",     metrics_cols)),
         all_of(paste0("p_upper_", metrics_cols))) |>
  pivot_longer(
    cols          = c(starts_with("SES_"), starts_with("p_upper_")),
    names_to      = c(".value", "metric"),
    names_pattern = "(SES|p_upper)_(.*)"
  ) |>
  group_by(metaweb, metric) |>                      # correct within metric
  mutate(padj = p.adjust(p_upper, method = "BH")) |>
  ungroup()

write.csv(ses_all,        here("results/SES_results_full.csv"),    row.names = FALSE)
write.csv(ses_long_clean, here("results/SES_results_long.csv"),    row.names = FALSE)

# Significant results
cat("\n=== Significant results (padj < 0.05) ===\n")
ses_long_clean |>
  filter(padj < 0.05) |>
  select(lake, metaweb, metric, SES, p_upper, padj) |>
  arrange(metric, padj) |>
  print(n = Inf)


# =============================================================================
# VISUALISATION
# =============================================================================

# --- 13a. SES heatmap ---

ses_long_clean_renamed <- ses_long_clean

ses_long_clean_renamed$lake <- gsub("Lake_", "", ses_long_clean_renamed$lake)

ses_long_clean_renamed$lake <- gsub("_Tarn", " Tarn", ses_long_clean_renamed$lake)

ses_long_clean_renamed$metric <- gsub("connectance", "SES Connectance", ses_long_clean_renamed$metric)
ses_long_clean_renamed$metric <- gsub("modularity", "SES Modularity", ses_long_clean_renamed$metric)
ses_long_clean_renamed$metric <- gsub("prop_basal", "SES % Basal", ses_long_clean_renamed$metric)
ses_long_clean_renamed$metric <- gsub("prop_top", "SES % Top", ses_long_clean_renamed$metric)
ses_long_clean_renamed$metric <- gsub("trophic_level_mean", "SES Mean Trophic Level", ses_long_clean_renamed$metric)
ses_long_clean_renamed$metric <- gsub("trophic_level_max", "SES Max Trophic Level", ses_long_clean_renamed$metric)
ses_long_clean_renamed$metric <- gsub("mean_chain_length", "SES Mean Chain Length", ses_long_clean_renamed$metric)
ses_long_clean_renamed$metric <- gsub("max_chain_length", "SES Max Chain Length", ses_long_clean_renamed$metric)


ses_long_clean_renamed$metric <- factor(ses_long_clean_renamed$metric, levels = c("SES Connectance", "SES Modularity", "SES % Basal", 
                   "SES % Top", "SES Mean Trophic Level", 
                   "SES Mean Chain Length", "SES Max Chain Length"))

p_ses_heatmap <- ses_long_clean_renamed |>
  filter(metric %in% c("SES Connectance", "SES Modularity", "SES % Basal", 
                   "SES % Top", "SES Mean Trophic Level",
                   "SES Mean Chain Length", "SES Max Chain Length")) |>
  mutate(
    sig     = padj < 0.05,
    metaweb = factor(metaweb, levels = c("high", "complete"),
                     labels = c("HCL", "AC"))
  ) |>
  ggplot(aes(x = metric, y = lake, fill = SES)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(data = ~ filter(., sig),
            aes(label = "*"), colour = "black", size = 3.5, vjust = 0.75) +
  scale_fill_gradient2(
    low      = "#175c7d",
    mid      = "white",
    high     = "#c22c43",
    midpoint = 0,
    limits   = c(-4, 4),
    oob      = scales::squish,
    name     = "SES"
  ) +
  facet_wrap(~metaweb) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x      = element_text(angle = 40, hjust = 1),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 16)
  ) +
  labs(x = NULL, y = NULL)

ggsave(here("results/SES_heatmap.jpeg"),
       p_ses_heatmap, width = 12, height = 7, dpi = 300, bg = "white")

# --- 13b. Uncertainty envelope: high vs complete metaweb ---

p_ses_envelope <- ses_long_clean_renamed |>
  filter(metric %in% c("SES Connectance", "SES Modularity", "SES % Basal", 
                   "SES % Top", "SES Mean Trophic Level",
                   "SES Mean Chain Length", "SES Max Chain Length")) |>
  mutate(metaweb = factor(metaweb, levels = c("high", "complete"),
                          labels = c("HCL", "AC"))) |>
  ggplot(aes(x = SES, y = lake, colour = metaweb, shape = metaweb)) +
  geom_vline(xintercept = c(-1.96, 1.96),
             linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
  geom_line(aes(group = lake), colour = "grey70", linewidth = 0.5) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_colour_manual(values = c("HCL" = "#c22c43",
                                 "AC"       = "#175c7d"),
                      name = "Metaweb") +
  scale_shape_manual(values = c("HCL" = 16, "AC" = 17),
                     name = "Metaweb") +
  facet_wrap(~metric, scales = "free_x", nrow = 1) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 16),
    legend.position  = "bottom"
  ) +
  labs(x = "SES", y = NULL)

ggsave(here("results/SES_envelope.jpeg"),
       p_ses_envelope, width = 20, height = 5, dpi = 300, bg = "white")

# --- 13c. Null distributions with observed value (per lake, key metrics) ---

for (m in c("connectance",
            "modularity",
            "prop_basal",
            "prop_top",
            "trophic_level_mean",
            "mean_chain_length",
            "max_chain_length")) {
  p <- plot_null_dist(null_high, observed_high,
                      metric = m, metaweb_label = "HCL")
  ggsave(here(paste0("results/null_dist_", m, "_high.jpeg")),
         p, width = 14, height = 10, dpi = 300, bg = "white")
}

for (m in c("connectance",
            "modularity",
            "prop_basal",
            "prop_top",
            "trophic_level_mean",
            "mean_chain_length",
            "max_chain_length")) {
  p <- plot_null_dist(null_complete, observed_complete,
                      metric = m, metaweb_label = "AC")
  ggsave(here(paste0("results/null_dist_", m, "_complete.jpeg")),
         p, width = 14, height = 10, dpi = 300, bg = "white")
}



################Constrained models#############



lake_elevation_lookup <- tibble::tribble(
  ~pa_column,             ~sd_lake,
  "Lake_45434",           "L45434",
  "Lake_46774",           "L46774",
  "Lake_54655",           "L54655",
  "Lake_Alexandrina",     "Alexandrina",
  "Alpine_Lake",          "Alpine",
  "Lake_Douglas",         "Douglas",
  "Lake_Dumbbell",        "Dumb-bell",
  "Lake_Duncan_Stream_1", "Duncan Stream 1",
  "Lake_Duncan_Stream_2", "Duncan Stream 2",
  "Lake_Duncan_Stream_3", "Duncan Stream 3",
  "Lake_Eggeling",        "Eggeling",
  "Lake_Gault",           "Gault",
  "Lake_Greaney",         "Greaney",
  "Lake_Hartley_Tarn",    "Hartley Tarn",
  "Lake_Haupiri",         "Haupiri",
  "Lake_Hochstetter",     "Hochstetter",
  "Lake_Ianthe",          "Ianthe",
  "Lake_Irish_Stream",    "Irish Stream",
  "Lake_Kangaroo",        "Kangaroo",
  "Lake_Kaurapataka",     "Kaurapataka",
  "Lake_Lady",            "Lady",
  "Little_Lake_Mckay",    "Little Mckay",
  "Lake_Lochnagar",       "Lochnagar",
  "Lake_Mahinapua",       "Mahinapua",
  "Lake_Mapourika",       "Mapourika",
  "Lake_Mcgregor",        "Mcgregor",
  "Lake_Mckay",           "Mckay",
  "Lake_Moeraki",         "Moeraki",
  "Lake_Morgan",          "Morgan",
  "Lake_Paringa",         "Paringa",
  "Lake_Poerua",          "Poerua",
  "Skifington_Swamp",     "Skifington Swamp",
  "Stoney_Tarn",          "Stoney Tarn"
) |>
  left_join(sd.plankton |> select(Lake, `Elevation.Type`),
            by = c("sd_lake" = "Lake")) |>
  rename(Elevation_Type = `Elevation.Type`)

stopifnot(
  "Some lakes did not match an Elevation Type in SupplementaryTable1" =
    all(!is.na(lake_elevation_lookup$Elevation_Type)),
  "Lookup table does not cover every lake column in pa_table" =
    setequal(lake_elevation_lookup$pa_column, colnames(pa_table))
)

alpine_lakes  <- lake_elevation_lookup$pa_column[lake_elevation_lookup$Elevation_Type == "Alpine"]
lowland_lakes <- lake_elevation_lookup$pa_column[lake_elevation_lookup$Elevation_Type == "Lowland"]

# Local-web subset: the same 15 lakes used throughout build_metaweb_public.R
local_web_lakes <- c(
  "Lake_Alexandrina", "Lake_Douglas", "Lake_Hartley_Tarn",
  "Lake_Hochstetter",  "Lake_Haupiri", "Lake_Ianthe",
  "Lake_Kangaroo",     "Lake_Lady",    "Lake_Mahinapua",
  "Lake_Mapourika",    "Lake_Mcgregor","Lake_Moeraki",
  "Lake_Paringa",      "Lake_Poerua",  "Stoney_Tarn"
)

alpine_local_lakes  <- intersect(alpine_lakes,  local_web_lakes)
lowland_local_lakes <- intersect(lowland_lakes, local_web_lakes)

stopifnot(
  "local_web_lakes not fully covered by alpine_lakes / lowland_lakes" =
    setequal(c(alpine_local_lakes, lowland_local_lakes), local_web_lakes)
)

cat("\nAlpine lakes used for local webs (",  length(alpine_local_lakes),  "):\n", sep = "");  print(alpine_local_lakes)
cat("\nLowland lakes used for local webs (", length(lowland_local_lakes), "):\n", sep = ""); print(lowland_local_lakes)

# =============================================================================
# 2. Subset the PA table by Elevation Type (full habitat lake set, not just
#    the 15 local-web lakes - this defines each habitat's species pool AND
#    the basis for habitat-specific prevalence weights)
# =============================================================================

subset_pa_table <- function(pa_table, lakes) {
  sub <- pa_table[, lakes, drop = FALSE]
  sub[rowSums(sub) > 0, , drop = FALSE]
}

pa_table_alpine  <- subset_pa_table(pa_table, alpine_lakes)
pa_table_lowland <- subset_pa_table(pa_table, lowland_lakes)

cat("\nAlpine PA table:  ",  nrow(pa_table_alpine),  "taxa x", ncol(pa_table_alpine),  "lakes\n")
cat("Lowland PA table: ",    nrow(pa_table_lowland), "taxa x", ncol(pa_table_lowland), "lakes\n")

taxa_alpine  <- rownames(pa_table_alpine)
taxa_lowland <- rownames(pa_table_lowland)

# =============================================================================
# 3. Subset master_taxa to match, filtering consumer-list columns too
# =============================================================================

filter_consumer_list <- function(x, keep) {
  if (is.na(x) || x == "") return(x)
  spp <- str_split(x, ";\\s*")[[1]]
  spp <- spp[spp %in% keep]
  if (length(spp) == 0) return(NA_character_)
  paste(spp, collapse = ";")
}

build_master_taxa_subset <- function(master_taxa, taxa_keep) {
  master_taxa |>
    filter(taxa_name %in% taxa_keep) |>
    rowwise() |>
    mutate(
      highest_potential_consumers = filter_consumer_list(highest_potential_consumers, taxa_keep),
      all_consumers               = filter_consumer_list(all_consumers, taxa_keep)
    ) |>
    ungroup()
}

master_taxa_alpine  <- build_master_taxa_subset(master_taxa, taxa_alpine)
master_taxa_lowland <- build_master_taxa_subset(master_taxa, taxa_lowland)

# =============================================================================
# 4. Rebuild metawebs (HCL + AC), restricted to each habitat's taxa
# =============================================================================

build_restricted_metaweb <- function(master_taxa_sub, all_consumers) {
  el <- build_edge_list(master_taxa = master_taxa_sub,
                        all_consumers = all_consumers,
                        test = FALSE)$edge_list
  mw <- build_metaweb(edge_list = el, master_taxa = master_taxa_sub, plot = FALSE)
  g  <- build_igraph(mw$adjacency_matrix)
  V(g)$name <- V(g)$name |> str_replace_all(" ", "_") |> str_trim()
  g
}

metaweb_alpine_high      <- build_restricted_metaweb(master_taxa_alpine,  all_consumers = FALSE)
metaweb_alpine_complete  <- build_restricted_metaweb(master_taxa_alpine,  all_consumers = TRUE)
metaweb_lowland_high     <- build_restricted_metaweb(master_taxa_lowland, all_consumers = FALSE)
metaweb_lowland_complete <- build_restricted_metaweb(master_taxa_lowland, all_consumers = TRUE)

cat("\nAlpine metaweb (HCL):   ", vcount(metaweb_alpine_high),      "nodes,", ecount(metaweb_alpine_high),      "edges\n")
cat("Alpine metaweb (AC):    ",   vcount(metaweb_alpine_complete),  "nodes,", ecount(metaweb_alpine_complete),  "edges\n")
cat("Lowland metaweb (HCL):  ",   vcount(metaweb_lowland_high),     "nodes,", ecount(metaweb_lowland_high),     "edges\n")
cat("Lowland metaweb (AC):   ",   vcount(metaweb_lowland_complete), "nodes,", ecount(metaweb_lowland_complete), "edges\n")

# =============================================================================
# 5. Rebuild local food webs (restricted to the 15 local_web_lakes)
# =============================================================================

pa_table_t_full <- t(pa_table) |> as.data.frame()
colnames(pa_table_t_full) <- colnames(pa_table_t_full) |>
  str_replace_all(" ", "_") |>
  str_trim()

pa_table_t_alpine  <- pa_table_t_full[alpine_local_lakes,  , drop = FALSE]
pa_table_t_lowland <- pa_table_t_full[lowland_local_lakes, , drop = FALSE]

local_webs_alpine_high      <- build_local_webs(pa_table_t_alpine,  metaweb_alpine_high,      "Alpine - HCL")
local_webs_alpine_complete  <- build_local_webs(pa_table_t_alpine,  metaweb_alpine_complete,  "Alpine - AC")
local_webs_lowland_high     <- build_local_webs(pa_table_t_lowland, metaweb_lowland_high,     "Lowland - HCL")
local_webs_lowland_complete <- build_local_webs(pa_table_t_lowland, metaweb_lowland_complete, "Lowland - AC")

# =============================================================================
# 6. Observed network metrics for the restricted local webs
# =============================================================================

compute_observed <- function(local_webs, metaweb_label) {
  map_dfr(names(local_webs), function(lake) {
    compute_network_metrics(local_webs[[lake]], lake) |> mutate(lake = lake)
  }) |> mutate(metaweb = metaweb_label)
}

observed_alpine_high      <- compute_observed(local_webs_alpine_high,      "high")
observed_alpine_complete  <- compute_observed(local_webs_alpine_complete,  "complete")
observed_lowland_high     <- compute_observed(local_webs_lowland_high,     "high")
observed_lowland_complete <- compute_observed(local_webs_lowland_complete, "complete")

# =============================================================================
# 7. Habitat-specific prevalence weights and shared taxonomic group
# =============================================================================
# Prevalence computed WITHIN each habitat's full lake set (not just the 15
# local-web lakes) - "how common is this species among alpine lakes
# generally", Laplace-smoothed so no species gets a hard-zero probability.

species_prevalence_alpine <- rowSums(pa_table_alpine) / ncol(pa_table_alpine)
species_prevalence_alpine <- species_prevalence_alpine + (0.5 / ncol(pa_table_alpine))
names(species_prevalence_alpine) <- rownames(pa_table_alpine)

species_prevalence_lowland <- rowSums(pa_table_lowland) / ncol(pa_table_lowland)
species_prevalence_lowland <- species_prevalence_lowland + (0.5 / ncol(pa_table_lowland))
names(species_prevalence_lowland) <- rownames(pa_table_lowland)

#species_group <- setNames(master_taxa$coarse_group, master_taxa$taxa_name)
#species_group[is.na(species_group) | species_group == ""] <- "Unknown"

# =============================================================================
# 8. Constrained null model function (habitat pool is implicit in whichever
#    metaweb is passed in; weights/groups add prevalence + taxonomic
#    composition constraints on top)
# =============================================================================

run_null_model_constrained <- function(metaweb, local_webs, label, n_null = 999,
                                       weights = NULL, groups = NULL) {
  
  cat("\nRunning CONSTRAINED null model for:", label, "(n =", n_null, ")\n")
  all_spp <- V(metaweb)$name
  
  w <- if (!is.null(weights)) weights[all_spp] else NULL
  if (!is.null(w)) {
    w[is.na(w) | w <= 0] <- min(w[w > 0], na.rm = TRUE) / 10
  }
  g <- if (!is.null(groups)) groups[all_spp] else NULL
  if (!is.null(g)) g[is.na(g)] <- "Unknown"
  
  with_progress({
    p <- progressor(steps = length(local_webs))
    
    null_results <- map(names(local_webs), function(lake) {
      
      p(message = lake)
      S_obs <- vcount(local_webs[[lake]])
      if (S_obs < 3) return(NULL)
      
      obs_spp <- V(local_webs[[lake]])$name
      
      future_map_dfr(1:n_null, function(i) {
        
        if (!is.null(g)) {
          g_obs               <- g[obs_spp]
          g_obs[is.na(g_obs)] <- "Unknown"
          obs_group_counts    <- table(g_obs)
          
          rand_spp <- unlist(lapply(names(obs_group_counts), function(grp) {
            pool     <- all_spp[g == grp]
            n_needed <- obs_group_counts[[grp]]
            if (length(pool) < n_needed) pool <- all_spp
            if (!is.null(w)) {
              sample(pool, size = n_needed, replace = FALSE, prob = w[pool])
            } else {
              sample(pool, size = n_needed, replace = FALSE)
            }
          }))
          rand_spp <- unique(rand_spp)
          
        } else if (!is.null(w)) {
          rand_spp <- sample(all_spp, size = S_obs, replace = FALSE, prob = w)
        } else {
          rand_spp <- sample(all_spp, size = S_obs, replace = FALSE)
        }
        
        rand_web <- induced_subgraph(metaweb, vids = rand_spp)
        compute_network_metrics(rand_web, community_name = paste0("null_", i))
        
      }, .options = furrr_options(seed = TRUE)) |>
        mutate(iteration = 1:n_null)
      
    }) |> set_names(names(local_webs))
  })
  
  null_results
}

# =============================================================================
# 9. Run the constrained null model per habitat x metaweb
# =============================================================================

null_alpine_high_constrained <- run_null_model_constrained(
  metaweb_alpine_high, local_webs_alpine_high, "Alpine - HCL (prevalence + taxonomic composition)",
  n_null = 999, weights = species_prevalence_alpine#, groups = species_group
)

null_alpine_complete_constrained <- run_null_model_constrained(
  metaweb_alpine_complete, local_webs_alpine_complete, "Alpine - AC (prevalence + taxonomic composition)",
  n_null = 999, weights = species_prevalence_alpine#, groups = species_group
)

null_lowland_high_constrained <- run_null_model_constrained(
  metaweb_lowland_high, local_webs_lowland_high, "Lowland - HCL (prevalence + taxonomic composition)",
  n_null = 999, weights = species_prevalence_lowland#, groups = species_group
)

null_lowland_complete_constrained <- run_null_model_constrained(
  metaweb_lowland_complete, local_webs_lowland_complete, "Lowland - AC (prevalence + taxonomic composition)",
  n_null = 999, weights = species_prevalence_lowland#, groups = species_group
)

saveRDS(null_alpine_high_constrained,      here("results/null_distributions_alpine_high_constrained.rds"))
saveRDS(null_alpine_complete_constrained,  here("results/null_distributions_alpine_complete_constrained.rds"))
saveRDS(null_lowland_high_constrained,     here("results/null_distributions_lowland_high_constrained.rds"))
saveRDS(null_lowland_complete_constrained, here("results/null_distributions_lowland_complete_constrained.rds"))

# =============================================================================
# 10. SES under the habitat + prevalence constrained
#    null model, combined into a single table covering all 15 local_web_lakes
# =============================================================================

ses_alpine_high_constrained      <- compute_ses(observed_alpine_high,      null_alpine_high_constrained,      "high",     metrics_cols)
ses_alpine_complete_constrained  <- compute_ses(observed_alpine_complete,  null_alpine_complete_constrained,  "complete", metrics_cols)
ses_lowland_high_constrained     <- compute_ses(observed_lowland_high,     null_lowland_high_constrained,     "high",     metrics_cols)
ses_lowland_complete_constrained <- compute_ses(observed_lowland_complete, null_lowland_complete_constrained, "complete", metrics_cols)

ses_habitat_constrained <- bind_rows(
  ses_alpine_high_constrained  |> mutate(Elevation_Type = "Alpine"),
  ses_alpine_complete_constrained |> mutate(Elevation_Type = "Alpine"),
  ses_lowland_high_constrained    |> mutate(Elevation_Type = "Lowland"),
  ses_lowland_complete_constrained|> mutate(Elevation_Type = "Lowland")
)

write.csv(ses_habitat_constrained,
          here("results/SES_results_habitat_constrained.csv"), row.names = FALSE)

# =============================================================================
# 11. Comparison vs the original (fully random, non-habitat-restricted) null
#    model, for the SAME 15 lakes - diagnostic table + delta SES plot +
#    null_dist plots. 
# =============================================================================

ses_original_subset <- bind_rows(ses_high, ses_complete) |>
  filter(lake %in% local_web_lakes) |>
  mutate(pool = "original") |>
  left_join(lake_elevation_lookup |> select(pa_column, Elevation_Type),
            by = c("lake" = "pa_column"))

ses_constrained_labelled <- ses_habitat_constrained |>
  mutate(pool = "habitat_constrained")

ses_comparison <- bind_rows(ses_original_subset, ses_constrained_labelled)

write.csv(ses_comparison, here("results/SES_comparison_original_vs_habitat_constrained.csv"), row.names = FALSE)

# --- Diagnostic table: obs / null_mean / null_sd / SES / p_upper side by side

diagnostic_stats <- c("obs", "null_mean", "null_sd", "SES", "p_upper")

ses_diagnostic_wide <- ses_comparison |>
  select(lake, metaweb, pool, S, Elevation_Type,
         all_of(paste0("obs_",       metrics_cols)),
         all_of(paste0("null_mean_", metrics_cols)),
         all_of(paste0("null_sd_",   metrics_cols)),
         all_of(paste0("SES_",       metrics_cols)),
         all_of(paste0("p_upper_",   metrics_cols))) |>
  pivot_longer(
    cols          = c(starts_with("obs_"), starts_with("null_mean_"),
                      starts_with("null_sd_"), starts_with("SES_"),
                      starts_with("p_upper_")),
    names_to      = c(".value", "metric"),
    names_pattern = "(obs|null_mean|null_sd|SES|p_upper)_(.*)"
  ) |>
  pivot_wider(
    id_cols     = c(lake, metaweb, metric, Elevation_Type),
    names_from  = pool,
    values_from = all_of(diagnostic_stats),
    names_glue  = "{.value}_{pool}"
  ) |>
  mutate(delta_SES = SES_habitat_constrained - SES_original)

write.csv(ses_diagnostic_wide,
          here("results/SES_diagnostic_original_vs_habitat_constrained.csv"),
          row.names = FALSE)

cat("\n=== Diagnostic: connectance, original vs habitat-constrained null model ===\n")
ses_diagnostic_wide |>
  filter(metric == "connectance") |>
  arrange(Elevation_Type, metaweb, lake) |>
  print(n = Inf)

# --- Delta SES plot: coloured by metaweb, shaped by Elevation Type ----------

ses_diagnostic_wide$metric <- gsub("connectance", "SES Connectance", ses_diagnostic_wide$metric)
ses_diagnostic_wide$metric <- gsub("modularity", "SES Modularity", ses_diagnostic_wide$metric)
ses_diagnostic_wide$metric <- gsub("prop_basal", "SES % Basal", ses_diagnostic_wide$metric)
ses_diagnostic_wide$metric <- gsub("prop_top", "SES % Top", ses_diagnostic_wide$metric)
ses_diagnostic_wide$metric <- gsub("trophic_level_mean", "SES Mean Trophic Level", ses_diagnostic_wide$metric)
ses_diagnostic_wide$metric <- gsub("mean_chain_length", "SES Mean Chain Length", ses_diagnostic_wide$metric)
ses_diagnostic_wide$metric <- gsub("max_chain_length", "SES Max Chain Length", ses_diagnostic_wide$metric)

ses_diagnostic_wide <- ses_diagnostic_wide %>% filter(metric != "trophic_level_max")

ses_diagnostic_wide$metric <- factor(ses_diagnostic_wide$metric, levels = c("SES Connectance", "SES Modularity", "SES % Basal", 
                                                                            "SES % Top", "SES Mean Trophic Level",
                                                                            "SES Mean Chain Length", "SES Max Chain Length"))


p_delta_ses <- ses_diagnostic_wide |>
  mutate(metaweb = factor(metaweb, levels = c("high", "complete"),
                          labels = c("HCL", "AC"))) |>
  ggplot(aes(x = metric, y = delta_SES, colour = metaweb, shape = Elevation_Type)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.8) +
  scale_colour_manual(values = c("HCL" = "#c22c43", "AC" = "#175c7d"), name = "Metaweb") +
  scale_shape_manual(values = c("Alpine" = 16, "Lowland" = 17), name = "Elevation Type") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
  labs(x = NULL, y = "Delta SES (constrained - original)")

ggsave(here("results/SES_delta_original_vs_habitat_constrained.jpeg"),
       p_delta_ses, width = 12, height = 6, dpi = 300, bg = "white")




compute_metaweb_roles <- function(metaweb) {
  data.frame(species = V(metaweb)$name) |>
    mutate(
      indeg_mw  = degree(metaweb, mode = "out"),
      outdeg_mw = degree(metaweb, mode = "in"),
      metaweb_role = case_when(
        outdeg_mw == 0 & indeg_mw >  0 ~ "basal",
        indeg_mw  == 0 & outdeg_mw > 0 ~ "top",
        indeg_mw  >  0 & outdeg_mw > 0 ~ "intermediate",
        TRUE                            ~ "isolated"
      )
    )
}

compute_role_presence <- function(lakes, pa_table_t, metaweb, habitat_label, mw_label) {
  metaweb_roles <- compute_metaweb_roles(metaweb)
  
  map_dfr(lakes, function(lake) {
    present    <- colnames(pa_table_t)[pa_table_t[lake, ] == 1]
    present    <- intersect(present, V(metaweb)$name)
    lake_roles <- metaweb_roles |> filter(species %in% present)
    
    data.frame(
      lake              = lake,
      habitat           = habitat_label,
      metaweb           = mw_label,
      n_top_present     = sum(lake_roles$metaweb_role == "top"),
      n_top_metaweb     = sum(metaweb_roles$metaweb_role == "top"),
      pct_top_present   = round(100 * sum(lake_roles$metaweb_role == "top") /
                                  sum(metaweb_roles$metaweb_role == "top"), 1),
      n_basal_present   = sum(lake_roles$metaweb_role == "basal"),
      pct_basal_present = round(100 * sum(lake_roles$metaweb_role == "basal") /
                                  sum(metaweb_roles$metaweb_role == "basal"), 1),
      n_int_present     = sum(lake_roles$metaweb_role == "intermediate"),
      pct_int_present   = round(100 * sum(lake_roles$metaweb_role == "intermediate") /
                                  sum(metaweb_roles$metaweb_role == "intermediate"), 1)
    )
  })
}


metaweb_role_presence_constrained <- bind_rows(
  compute_role_presence(alpine_local_lakes,  pa_table_t_alpine,  metaweb_alpine_high,      "Alpine",  "high"),
  compute_role_presence(alpine_local_lakes,  pa_table_t_alpine,  metaweb_alpine_complete,  "Alpine",  "complete"),
  compute_role_presence(lowland_local_lakes, pa_table_t_lowland, metaweb_lowland_high,     "Lowland", "high"),
  compute_role_presence(lowland_local_lakes, pa_table_t_lowland, metaweb_lowland_complete, "Lowland", "complete")
)

ses_long_clean_constrained <- ses_habitat_constrained |>
  select(lake, metaweb, S,
         all_of(paste0("SES_",     metrics_cols)),
         all_of(paste0("p_upper_", metrics_cols))) |>
  pivot_longer(
    cols          = c(starts_with("SES_"), starts_with("p_upper_")),
    names_to      = c(".value", "metric"),
    names_pattern = "(SES|p_upper)_(.*)"
  ) |>
  filter(metric != "trophic_level_max") |>
  group_by(metaweb, metric) |>                      # correct within metric
  mutate(padj = p.adjust(p_upper, method = "BH")) |>
  ungroup()

ses_long_clean_constrained |>
  filter(padj < 0.05) |>
  select(lake, metaweb, metric, SES, p_upper, padj) |>
  arrange(metric, padj) |>
  print(n = Inf)



ses_long_clean_constrained_renamed <- ses_long_clean_constrained

ses_long_clean_constrained_renamed$lake <- gsub("Lake_", "", ses_long_clean_constrained_renamed$lake)

ses_long_clean_constrained_renamed$lake <- gsub("_Tarn", " Tarn", ses_long_clean_constrained_renamed$lake)

ses_long_clean_constrained_renamed$metric <- gsub("connectance", "SES Connectance", ses_long_clean_constrained_renamed$metric)
ses_long_clean_constrained_renamed$metric <- gsub("modularity", "SES Modularity", ses_long_clean_constrained_renamed$metric)
ses_long_clean_constrained_renamed$metric <- gsub("prop_basal", "SES % Basal", ses_long_clean_constrained_renamed$metric)
ses_long_clean_constrained_renamed$metric <- gsub("prop_top", "SES % Top", ses_long_clean_constrained_renamed$metric)
ses_long_clean_constrained_renamed$metric <- gsub("trophic_level_mean", "SES Mean Trophic Level", ses_long_clean_constrained_renamed$metric)
ses_long_clean_constrained_renamed$metric <- gsub("trophic_level_max", "SES Max Trophic Level", ses_long_clean_constrained_renamed$metric)
ses_long_clean_constrained_renamed$metric <- gsub("mean_chain_length", "SES Mean Chain Length", ses_long_clean_constrained_renamed$metric)
ses_long_clean_constrained_renamed$metric <- gsub("max_chain_length", "SES Max Chain Length", ses_long_clean_constrained_renamed$metric)


ses_long_clean_constrained_renamed$metric <- factor(ses_long_clean_constrained_renamed$metric, levels = c("SES Connectance", "SES Modularity", "SES % Basal", 
                                                                                                          "SES % Top", "SES Mean Trophic Level", 
                                                                                                          "SES Mean Chain Length", "SES Max Chain Length"))

p_ses_heatmap <- ses_long_clean_constrained_renamed |>
  filter(metric %in% c("SES Connectance", "SES Modularity", "SES % Basal", 
                       "SES % Top", "SES Mean Trophic Level",
                       "SES Mean Chain Length", "SES Max Chain Length")) |>
  mutate(
    sig     = padj < 0.05,
    metaweb = factor(metaweb, levels = c("high", "complete"),
                     labels = c("HCL", "AC"))
  ) |>
  ggplot(aes(x = metric, y = lake, fill = SES)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(data = ~ filter(., sig),
            aes(label = "*"), colour = "black", size = 3.5, vjust = 0.75) +
  scale_fill_gradient2(
    low      = "#175c7d",
    mid      = "white",
    high     = "#c22c43",
    midpoint = 0,
    limits   = c(-4, 4),
    oob      = scales::squish,
    name     = "SES"
  ) +
  facet_wrap(~metaweb) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x      = element_text(angle = 40, hjust = 1),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 16)
  ) +
  labs(x = NULL, y = NULL)

ggsave(here("results/SES_heatmap_constrained.jpeg"),
       p_ses_heatmap, width = 12, height = 7, dpi = 300, bg = "white")


null_dist_metrics <- c(
  "connectance", "modularity", "prop_basal", "prop_top",
  "trophic_level_mean", "mean_chain_length", "max_chain_length"
)

null_dist_jobs <- list(
  list(null = null_alpine_high_constrained,      obs = observed_alpine_high,      tag = "alpine_high",      label = "Alpine - HCL (constrained)"),
  list(null = null_alpine_complete_constrained,  obs = observed_alpine_complete,  tag = "alpine_complete",  label = "Alpine - AC (constrained)"),
  list(null = null_lowland_high_constrained,     obs = observed_lowland_high,     tag = "lowland_high",     label = "Lowland - HCL (constrained)"),
  list(null = null_lowland_complete_constrained, obs = observed_lowland_complete, tag = "lowland_complete", label = "Lowland - AC (constrained)")
)

for (job in null_dist_jobs) {
  for (m in null_dist_metrics) {
    p <- plot_null_dist(job$null, job$obs, metric = m, metaweb_label = job$label)
    ggsave(here(paste0("results/null_dist_", m, "_", job$tag, "_constrained.jpeg")),
           p, width = 14, height = 10, dpi = 300, bg = "white")
  }
}

# =============================================================================
# 12. Does the SIGNIFICANCE of the trophic-structuring signal survive under
#    the habitat + prevalence + taxonomic composition constrained null
#    model? Reported both pooled and split by Elevation_Type.
#    Addresses reviewer comment (2) directly.
# =============================================================================

sig_summary <- function(ses_df, metrics_cols, pool_label, group_col = NULL) {
  extra_groups <- c("metaweb", group_col)
  ses_df |>
    group_by(across(all_of(extra_groups))) |>
    group_modify(~ map_dfr(metrics_cols, function(m) {
      vals <- .x[[paste0("SES_", m)]]
      data.frame(
        metric      = m,
        pool        = pool_label,
        n_lakes     = sum(!is.na(vals)),
        mean_SES    = mean(vals, na.rm = TRUE),
        pct_sig_pos = mean(vals >  1.96, na.rm = TRUE) * 100,
        pct_sig_neg = mean(vals < -1.96, na.rm = TRUE) * 100,
        pct_sig_any = mean(abs(vals) > 1.96, na.rm = TRUE) * 100
      )
    })) |>
    ungroup()
}

# Pooled (all 15 lakes together), split by metaweb
sig_comparison_pooled <- bind_rows(
  sig_summary(ses_original_subset,      metrics_cols, "original"),
  sig_summary(ses_habitat_constrained,  metrics_cols, "habitat_constrained")
)

write.csv(sig_comparison_pooled,
          here("results/SES_significance_comparison_pooled.csv"),
          row.names = FALSE)

# Split by metaweb AND Elevation_Type
sig_comparison_by_habitat <- bind_rows(
  sig_summary(ses_original_subset,     metrics_cols, "original",            group_col = "Elevation_Type"),
  sig_summary(ses_habitat_constrained, metrics_cols, "habitat_constrained", group_col = "Elevation_Type")
)

write.csv(sig_comparison_by_habitat,
          here("results/SES_significance_comparison_by_habitat.csv"),
          row.names = FALSE)

cat("\n=== % of lakes with significant SES (|SES| > 1.96): pooled (by metaweb), original vs habitat-constrained ===\n")
sig_comparison_pooled |> arrange(metaweb, metric, pool) |> print(n = Inf)

cat("\n=== Same, split by metaweb AND Elevation Type ===\n")
sig_comparison_by_habitat |> arrange(metaweb, metric, Elevation_Type, pool) |> print(n = Inf)

cat("\n=== Connectance specifically ===\n")
sig_comparison_by_habitat |> filter(metric == "connectance") |> arrange(metaweb, Elevation_Type, pool) |> print(n = Inf)


# --- 13d. Richness vs connectance SES (check for richness artefact) ---

ses_long_clean_renamed <- ses_long_clean

ses_long_clean_renamed$lake <- gsub("Lake_", "", ses_long_clean_renamed$lake)

ses_long_clean_renamed$lake <- gsub("_Tarn", " Tarn", ses_long_clean_renamed$lake)

ses_long_clean_renamed$metric <- gsub("connectance", "Connectance", ses_long_clean_renamed$metric)
ses_long_clean_renamed$metric <- gsub("modularity", "Modularity", ses_long_clean_renamed$metric)
ses_long_clean_renamed$metric <- gsub("prop_basal", "% Basal", ses_long_clean_renamed$metric)
ses_long_clean_renamed$metric <- gsub("prop_top", "% Top", ses_long_clean_renamed$metric)
ses_long_clean_renamed$metric <- gsub("trophic_level_mean", "Mean Trophic Level", ses_long_clean_renamed$metric)
ses_long_clean_renamed$metric <- gsub("mean_chain_length", "Mean Chain Length", ses_long_clean_renamed$metric)
ses_long_clean_renamed$metric <- gsub("max_chain_length", "Max Chain Length", ses_long_clean_renamed$metric)

p_richness_sesconnectance <- ses_long_clean_renamed |>
  filter(metric == "Connectance") |>
  mutate(metaweb = factor(metaweb, levels = c("high", "complete"),
                          labels = c("HCL", "AC"))) |>
  ggplot(aes(x = S, y = SES, colour = metaweb, label = lake)) +
 # geom_hline(yintercept = c(-1.96, 1.96),
 #            linetype = "dashed", colour = "grey60") +
  geom_hline(yintercept = 0, colour = "black") +
  geom_smooth(method = "lm", se = F, alpha = 0.15) +
  geom_point(size = 3) +
  #ggrepel::geom_text_repel(size = 2.5, show.legend = FALSE) +
  scale_colour_manual(values = c("HCL" = "#c22c43",
                                 "AC"       = "#175c7d")) +
  theme_classic(base_size = 10) +
  labs(x = "Species richness (S)", y = "SES Connectance") + 
  theme(legend.title = element_blank(),
        legend.position = "",
        axis.title.x = element_blank(),
        axis.text = element_text(size=12),
        axis.title.y = element_text(size=14),
        plot.title = element_text(size=15)) +
  ggtitle("A")

ggsave(here("results/richness_vs_SES_connectance.jpeg"),
       p_richness_sesconnectance, width = 8, height = 5, dpi = 300, bg = "white")

cat("\nDone. Outputs written to results/ and results/\n")


df_model <- ses_long_clean_renamed |>
  filter(metric == "Connectance") |>
  mutate(metaweb = factor(metaweb,
                          levels = c("high", "complete"),
                          labels = c("HCL", "AC")))

lm_list <- df_model |>
  split(df_model$metaweb) |>
  lapply(function(d) lm(SES ~ S, data = d))


lapply(lm_list, summary)


p_richness_sesmodularity <- ses_long_clean_renamed |>
  filter(metric == "Modularity") |>
  mutate(metaweb = factor(metaweb, levels = c("high", "complete"),
                          labels = c("HCL", "AC"))) |>
  ggplot(aes(x = S, y = SES, colour = metaweb, label = lake)) +
 # geom_hline(yintercept = c(-1.96, 1.96),
 #            linetype = "dashed", colour = "grey60") +
  geom_hline(yintercept = 0, colour = "black") +
  geom_smooth(method = "lm", se = F, alpha = 0.15) +
  geom_point(size = 3) +
 # ggrepel::geom_text_repel(size = 2.5, show.legend = FALSE) +
  scale_colour_manual(values = c("HCL" = "#c22c43",
                                 "AC"       = "#175c7d")) +
  theme_classic(base_size = 10) +
  labs(x = "Species richness (S)", y = "SES Modularity") + 
  theme(legend.title = element_blank(),
        legend.position = "",
        axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        plot.title = element_text(size=15)) +
  ggtitle("B")

ggsave(here("results/richness_vs_SES_modularity.jpeg"),
       p_richness_sesmodularity, width = 8, height = 5, dpi = 300, bg = "white")


df_model <- ses_long_clean_renamed |>
  filter(metric == "Modularity") |>
  mutate(metaweb = factor(metaweb,
                          levels = c("high", "complete"),
                          labels = c("HCL", "AC")))

lm_list <- df_model |>
  split(df_model$metaweb) |>
  lapply(function(d) lm(SES ~ S, data = d))


lapply(lm_list, summary)

p_richness_sestrophic_level_mean <- ses_long_clean_renamed |>
  filter(metric == "Mean Trophic Level") |>
  mutate(metaweb = factor(metaweb, levels = c("high", "complete"),
                          labels = c("HCL", "AC"))) |>
  ggplot(aes(x = S, y = SES, colour = metaweb, label = lake)) +
 # geom_hline(yintercept = c(-1.96, 1.96),
 #            linetype = "dashed", colour = "grey60") +
  geom_hline(yintercept = 0, colour = "black") +
  geom_smooth(method = "lm", se = F, alpha = 0.15) +
  geom_point(size = 3) +
 # ggrepel::geom_text_repel(size = 2.5, show.legend = FALSE) +
  scale_colour_manual(values = c("HCL" = "#c22c43",
                                 "AC"       = "#175c7d")) +
  theme_classic(base_size = 10) +
  labs(x = "Species richness (S)", y = "SES Mean Trophic Level") + 
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size=12),
        axis.title.y = element_text(size=14),
        plot.title = element_text(size=15),
        legend.text = element_text(size=14)) +
  ggtitle("C")

ggsave(here("results/richness_vs_SES_trophic_level_mean.jpeg"),
       p_richness_sestrophic_level_mean, width = 8, height = 5, dpi = 300, bg = "white")


df_model <- ses_long_clean_renamed |>
  filter(metric == "Mean Trophic Level") |>
  mutate(metaweb = factor(metaweb,
                          levels = c("high", "complete"),
                          labels = c("HCL", "AC")))

lm_list <- df_model |>
  split(df_model$metaweb) |>
  lapply(function(d) lm(SES ~ S, data = d))


lapply(lm_list, summary)

combined_SES_richness <- plot_grid(p_richness_sesconnectance, p_richness_sesmodularity, p_richness_sestrophic_level_mean, nrow = 1, align = "hv")

ggsave(here("results/combined_SES_richness.jpeg"),
       combined_SES_richness, width = 15, height = 5, dpi = 300, bg = "white")






# --- Null distribution plots (constrained model), split by habitat ----------


# =============================================================================
# 13. Species richness (S) vs SES regression - POOLED ACROSS ALL 15 LAKES,
#    NOT split by habitat, since this is the test for reviewer comment (1)
#    and needs the power. Compares original (unconstrained) vs the habitat +
#    prevalence + taxonomic composition constrained null model.
# =============================================================================

richness_test <- function(ses_df, metrics_cols, pool_label, metaweb_label) {
  map_dfr(metrics_cols, function(m) {
    f   <- as.formula(paste0("SES_", m, " ~ S"))
    fit <- tryCatch(lm(f, data = ses_df), error = function(e) NULL)
    if (is.null(fit) || nrow(fit$model) < 3) return(NULL)
    s <- summary(fit)
    data.frame(
      metric    = m,
      metaweb   = metaweb_label,
      pool      = pool_label,
      n         = nrow(fit$model),
      slope     = coef(fit)[["S"]],
      p_value   = s$coefficients["S", "Pr(>|t|)"],
      r_squared = s$r.squared
    )
  })
}

richness_comparison <- bind_rows(
  richness_test(ses_original_subset     |> filter(metaweb == "high"),     metrics_cols, "original",            "high"),
  richness_test(ses_original_subset     |> filter(metaweb == "complete"), metrics_cols, "original",            "complete"),
  richness_test(ses_habitat_constrained |> filter(metaweb == "high"),     metrics_cols, "habitat_constrained", "high"),
  richness_test(ses_habitat_constrained |> filter(metaweb == "complete"), metrics_cols, "habitat_constrained", "complete")
)

write.csv(richness_comparison,
          here("results/richness_SES_regression_comparison.csv"),
          row.names = FALSE)


p_ses_vs_S <- bind_rows(
  ses_original_subset      |> mutate(pool = "Full"),
  ses_habitat_constrained  |> mutate(pool = "Constrained")
) |>
  select(lake, metaweb, S, pool, all_of(paste0("SES_", metrics_cols))) |>
  pivot_longer(starts_with("SES_"), names_to = "metric", values_to = "SES") |>
  mutate(metric = str_remove(metric, "SES_")) |>
  filter(metric %in% c("connectance", "modularity", "trophic_level_mean")) |>
  mutate(
    metric  = factor(metric,
                     levels = c("connectance", "modularity", "trophic_level_mean"),
                     labels = c("SES Connectance", "SES Modularity", "SES Mean Trophic Level")),
    metaweb = factor(metaweb, levels = c("high", "complete"), labels = c("HCL", "AC")),
    pool    = factor(pool, levels = c("Full", "Constrained"))
  ) |>
  ggplot(aes(x = S, y = SES, colour = metaweb)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = F, aes(fill = metaweb), alpha = 0.15) +
  facet_grid(pool ~ metric, scales = "fixed") +
  scale_colour_manual(values = c("HCL" = "#c22c43", "AC" = "#175c7d"), name = "Metaweb") +
  scale_fill_manual(values = c("HCL" = "#c22c43", "AC" = "#175c7d"), guide = "none") +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 40, hjust = 1),
        strip.text.y = element_text(angle = 0)) +
  labs(x = "Species richness (S)", y = "SES")

ggsave(here("results/diagnostic_SES_vs_S_comparison.jpeg"),
       p_ses_vs_S, width = 14, height = 10, dpi = 300, bg = "white")


####Environmental gradients

response_df <- observed_high |>
  mutate(lake = as.character(lake)) |>
  distinct(lake, .keep_all = TRUE)

cat("\nHCL observed metrics:\n")
print(response_df)

# =============================================================================
# 2. GET ENVIRONMENTAL DATA
# =============================================================================

env_df <- lake_elevation_lookup |>
  filter(pa_column %in% response_df$lake) |>
  dplyr::select(pa_column, sd_lake) |>
  left_join(sd.plankton, by = c("sd_lake" = "Lake")) |>
  rename(lake = pa_column) |>
  arrange(lake) |>
  transmute(
    lake,
    
    # Physical
    depth       = `Max..Depth..m.`,
    elevation   = `Elevation..m.`,
    clarity     = `Secchi.Disk.Depth..m.`,
    
    # Catchment land use
    pasture_pct = `X..High.productivity.Grassland.in.Catchment`,
    native_pct  = `X..Native..non.forest..in.Catchment`,
    
    # Water chemistry
    TN          = `Total.Nitrogen..g.m3.`,
    TP          = `Total.Phosphorus..g.m3.`,
    TOC         = `Total.Organic.Carbon..g.m3.`,
    chl_a       = `Chlorophyll.a...g.m3.`
  )

cat("\nEnvironmental data:\n")
print(env_df)

# =============================================================================
# 3. JOIN NETWORK METRICS + ENVIRONMENTAL DATA
# =============================================================================

analysis_df <- response_df |>
  left_join(env_df, by = "lake")

cat("\nNumber of HCL lakes:", nrow(analysis_df), "\n")

# Check for missing environmental data
cat("\nMissing environmental values:\n")

analysis_df |>
  summarise(
    across(
      c(depth, elevation, clarity, pasture_pct, native_pct,
        TN, TP, TOC, chl_a),
      ~ sum(is.na(.x))
    )
  ) |>
  print()

# Save combined dataset
write.csv(
  analysis_df,
  here("results/HCL_raw_metrics_environmental_data.csv"),
  row.names = FALSE
)

# =============================================================================
# 4. DEFINE RESPONSE VARIABLES
# =============================================================================

# These are RAW observed HCL metrics.
#
# DO NOT replace these with SES_* variables.

response_metrics <- c("S",
                      "connectance",
                      "modularity",
                      # "prop_basal",
                      # "prop_top",
                      "trophic_level_mean"#,
                      # "trophic_level_max",
                      # "mean_chain_length",
                      # "max_chain_length"
)

# Environmental predictors
environmental_vars <- c(
  "depth",
  "elevation",
  "clarity",
  "pasture_pct",
  "native_pct",
  "TN",
  "TP",
  "TOC"
)

# Check that everything exists
stopifnot(
  all(response_metrics %in% names(analysis_df)),
  all(environmental_vars %in% names(analysis_df))
)

# =============================================================================
# 5. STANDARDISE ENVIRONMENTAL VARIABLES
# =============================================================================
#
# Standardisation is useful for comparing effect sizes.
#
# IMPORTANT:
# The response variables remain RAW.
# Only environmental predictors are standardised.

analysis_df <- analysis_df |>
  mutate(
    across(
      all_of(environmental_vars),
      ~ as.numeric(scale(.x)),
      .names = "{.col}_z"
    )
  )

environmental_vars_z <- paste0(environmental_vars, "_z")

# =============================================================================
# 6. CALCULATE SUMS OF SQUARES
# =============================================================================


calculate_ss <- function(data, response, predictor) {
  
  dat <- data |>
    dplyr::select(
      response = all_of(response),
      predictor = all_of(predictor)
    ) |>
    dplyr::filter(
      is.finite(response),
      is.finite(predictor)
    )
  
  if (nrow(dat) < 3) {
    return(
      tibble(
        metric = response,
        environment = predictor,
        n = nrow(dat),
        SS_effect = NA_real_,
        SS_residual = NA_real_,
        SS_total = NA_real_,
        pct_SS = NA_real_,
        R2 = NA_real_,
        F = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  fit <- lm(response ~ predictor, data = dat)
  
  aov_tab <- anova(fit)
  
  # Raw sums of squares
  SS_effect_raw <- aov_tab$`Sum Sq`[1]
  SS_residual_raw <- aov_tab$`Sum Sq`[2]
  SS_total_raw <- SS_effect_raw + SS_residual_raw
  
  # ---------------------------------------------------------------
  # STANDARDISED SUMS OF SQUARES
  #
  # Divide every SS component by total SS.
  # Therefore:
  #
  #   SS_effect + SS_residual = 1
  #   SS_effect = R2
  # ---------------------------------------------------------------
  
  SS_effect_std <- SS_effect_raw / SS_total_raw
  SS_residual_std <- SS_residual_raw / SS_total_raw
  
  tibble(
    metric = response,
    environment = predictor,
    n = nrow(dat),
    
    # Standardised SS
    SS_effect = SS_effect_std,
    SS_residual = SS_residual_std,
    SS_total = 1,
    
    # Percentage explained
    pct_SS = 100 * SS_effect_std,
    
    # Same quantity as standardised effect SS
    R2 = SS_effect_std,
    
    F = aov_tab$`F value`[1],
    p_value = aov_tab$`Pr(>F)`[1]
  )
}
# Remove "_z" from predictor names
ss_results <- purrr::map_dfr(
  response_metrics,
  function(metric) {
    
    purrr::map_dfr(
      environmental_vars_z,
      function(env) {
        
        calculate_ss(
          data = analysis_df,
          response = metric,
          predictor = env
        )
        
      }
    )
  }
)


ss_results <- ss_results |>
  mutate(
    environment = str_remove(environment, "_z$")
  )

# =============================================================================
# 7. FDR CORRECTION
# =============================================================================
#
# Correct across environmental predictors within each network metric.

ss_results <- ss_results |>
  group_by(metric) |>
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  ) |>
  ungroup()

# =============================================================================
# 8. SAVE FULL SS TABLE
# =============================================================================

write.csv(
  ss_results,
  here("results/HCL_environmental_SS_results.csv"),
  row.names = FALSE
)

# Print results
cat("\n============================================================\n")
cat("HCL RAW METRICS — ENVIRONMENTAL SUMS OF SQUARES\n")
cat("============================================================\n\n")

ss_results |>
  arrange(metric, p_adj) |>
  print(n = Inf)

# =============================================================================
# 9. CREATE A CLEAN SUMMARY TABLE
# =============================================================================

ss_summary <- ss_results |>
  mutate(
    metric = dplyr::recode(
      metric,
      connectance         = "Connectance",
      modularity          = "Modularity",
      S = "Richness",
      trophic_level_mean  = "Mean trophic level"
    ),
    environment = dplyr::recode(
      environment,
      depth       = "Depth",
      elevation   = "Elevation",
      clarity     = "Clarity",
      pasture_pct = "HPG",
      native_pct  = "Native",
      TN          = "TN",
      TP          = "TP",
      TOC         = "TOC"
    )
  )

write.csv(
  ss_summary,
  here("results/HCL_environmental_SS_summary.csv"),
  row.names = FALSE
)

# =============================================================================
# 10. SS HEATMAP
# =============================================================================

# Heatmap shows percentage of TOTAL RAW METRIC variation explained
# by each environmental variable.

ss_summary$environment <- factor(ss_summary$environment, levels=c("TN", "TP", "TOC", "Clarity", "Depth", "Elevation", "HPG", "Native"))

p_ss_heatmap <- ss_summary |>
  ggplot(
    aes(
      x = environment,
      y = metric,
      fill = pct_SS
    )
  ) +
  geom_tile(
    colour = "white",
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      label = ifelse(
        p_adj < 0.05,
        paste0(round(pct_SS, 1), "%*"),
        round(pct_SS, 1)
      )
    ),
    size = 3.5
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#365C8DFF",
    limits = c(0, max(ss_summary$pct_SS, na.rm = TRUE)),
    name = "% SS"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 40,
      hjust = 1
    ),
    axis.title = element_blank(),
    legend.position = "right"
  ) +
  ggtitle("B")

ggsave(
  here("results/HCL_environmental_SS_heatmap.jpeg"),
  p_ss_heatmap,
  width = 12,
  height = 7,
  dpi = 300,
  bg = "white"
)




# =============================================================================
# 10. PREPARE SS PLOT DATA
# =============================================================================

ss_plot_data <- ss_results |>
  mutate(
    metric = dplyr::recode(
      metric,
      S                   = "Richness",
      connectance         = "Connectance",
      modularity          = "Modularity",
      trophic_level_mean  = "Mean trophic level"
    ),
    environment = dplyr::recode(
      environment,
      depth       = "Depth",
      elevation   = "Elevation",
      clarity     = "Clarity",
      pasture_pct = "HPG",
      native_pct  = "Native",
      TN          = "TN",
      TP          = "TP",
      TOC         = "TOC"
    )
  )

# Calculate mean SS for each metric
ss_means <- ss_plot_data |>
  group_by(metric) |>
  summarise(
    mean_SS = mean(SS_effect, na.rm = TRUE),
    sd_SS = sd(SS_effect, na.rm = TRUE),
    .groups = "drop"
  )



# =============================================================================
# 11. VIOLIN PLOT WITH MEAN SS
# =============================================================================

p_ss_violin <- ggplot(
  ss_plot_data,
  aes(
    x = metric,
    y = SS_effect
  )
) +
  geom_violin(
    fill = "#365C8DFF",
    colour = "black",
    alpha = 0.65,
    trim = TRUE
  ) +
  
  # Individual environmental variables
  geom_jitter(
    width = 0.08,
    height = 0,
    size = 2,
    alpha = 0.7
  ) +
  
  # Mean line
  geom_crossbar(
    data = ss_means,
    aes(
      x = metric,
      y = mean_SS,
      ymin = mean_SS,
      ymax = mean_SS
    ),
    width = 0.65,
    linewidth = 0.8,
    colour = "black"
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  coord_cartesian(
    ylim = c(0, NA)
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 40,
      hjust = 1
    ),
    axis.title.x = element_blank()
  ) +
  labs(
    y = "Standardised sum of squares (SS)"
  ) +
  ggtitle("A")





library(patchwork)

p_combined <- p_ss_violin +
  p_ss_heatmap +
  plot_layout(
    widths = c(1, 1.5)
  )

ggsave(
  here("results/HCL_environmental_SS_combined.jpeg"),
  p_combined,
  width = 16,
  height = 7,
  dpi = 300,
  bg = "white"
)

p_combined



# ---------------------------------------------------------------------------
# Compare with GLoBI
# ---------------------------------------------------------------------------



master_taxa <- read.csv(here("data/master_taxa.csv"))

master_taxa <- master_taxa %>%
  mutate(
    taxa_clean = taxa_name %>%
      str_remove("_sp\\.?$") %>%
      str_replace_all("_", " ") %>%
      str_trim()
  )

unique_taxa <- unique(master_taxa$taxa_clean)
length(unique_taxa)

get_all_interactions <- function(taxon, page_size = 10000) {
  skip <- 0
  pages <- list()
  
  repeat {
    otherkeys <- list(limit = page_size, skip = skip)
    
    page <- tryCatch(
      get_interactions_by_taxa(sourcetaxon = taxon, otherkeys = otherkeys),
      error = function(e) {
        message("  ! error for '", taxon, "': ", conditionMessage(e))
        NULL
      }
    )
    
    if (is.null(page) || nrow(page) == 0) break
    
    pages[[length(pages) + 1]] <- page
    
    if (nrow(page) < page_size) break  
    skip <- skip + page_size
  }
  
  if (length(pages) == 0) return(NULL)
  dplyr::bind_rows(pages)
}


all_results <- list()

for (taxon in unique_taxa) {
  message("Querying: ", taxon)
  res <- get_all_interactions(taxon)
  if (!is.null(res) && nrow(res) > 0) {
    res$queried_taxon <- taxon   
    all_results[[taxon]] <- res
  }
  Sys.sleep(0.2)  
}

all_results_fixed <- lapply(all_results, function(df) {
  df %>% mutate(across(everything(), as.character))
})

totalInteractions <- dplyr::bind_rows(all_results_fixed)

source_names_in_globi <- unique(totalInteractions$source_taxon_name)


is_binomial <- grepl(" ", unique_taxa)

species_names <- unique_taxa[is_binomial]
genus_names   <- unique_taxa[!is_binomial]

matched_species <- species_names[species_names %in% source_names_in_globi]

source_genus <- sub("^(\\S+).*", "\\1", source_names_in_globi)

matched_genus <- genus_names[sapply(genus_names, function(g) {
  any(source_genus == g)
})]

matched_taxa <- c(matched_species, matched_genus)

genus_match_details <- lapply(genus_names, function(g) {
  source_names_in_globi[source_genus == g]
})
names(genus_match_details) <- genus_names


n_matched <- length(matched_taxa)
n_total   <- length(unique_taxa)

cat(sprintf(
  "%d of %d unique taxa (%.1f%%) were found as source_taxon_name in GloBI interactions\n",
  n_matched, n_total, 100 * n_matched / n_total
))


matched_taxa
unmatched_taxa


master_taxa_matched <- master_taxa %>%
  mutate(globi_match = if_else(taxa_clean %in% source_names_in_globi,
                               "matched", "not matched"))


summary_results <- master_taxa_matched %>%
  group_by(coarse_group, globi_match) %>%
  tally() %>%
  group_by(coarse_group) %>%
  mutate(total = sum(n),
         percentage = round(100 * n / total, 1)) %>%
  ungroup()


summary_wide <- summary_results %>%
  dplyr::select(coarse_group, globi_match, n) %>%
  pivot_wider(names_from = globi_match, values_from = n, values_fill = 0) %>%
  left_join(
    summary_results %>% distinct(coarse_group, total),
    by = "coarse_group"
  ) %>%
  mutate(across(c(matched, `not matched`), 
                ~ round(100 * . / total, 1), 
                .names = "{.col}_pct"))