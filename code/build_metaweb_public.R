
library(dplyr)
library(tidyr)
library(here)
library(stringr)
library(cheddar)
library(cassandRa)
library(ggplot2)
library(purrr)
library(here)
library(tidyverse)
library(igraph)
library(cheddar)
library(cowplot)
library(furrr)
library(progressr)
library(marmap)
library(ggsn)
library(ggrepel)

plan(multisession, workers = 4)

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
test_lake   <- lake_names[4]
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
# METAWEB TROPHIC ROLE PRESENCE
# How much of each trophic role class is represented in each lake?
# =============================================================================

metaweb_roles <- data.frame(
  species   = V(metaweb_high)$name
) |>
  mutate(
    indeg_mw  = degree(metaweb_high, mode = "out"),
    outdeg_mw = degree(metaweb_high, mode = "in"),
    # Assumes consumer -> resource edge direction
    metaweb_role = case_when(
      outdeg_mw == 0 & indeg_mw >  0 ~ "basal",
      indeg_mw  == 0 & outdeg_mw > 0 ~ "top",
      indeg_mw  >  0 & outdeg_mw > 0 ~ "intermediate",
      TRUE                            ~ "isolated"
    )
  )

metaweb_role_presence <- map_dfr(lake_names, function(lake) {
  present    <- colnames(pa_table)[pa_table[lake, ] == 1]
  present    <- intersect(present, V(metaweb_high)$name)
  lake_roles <- metaweb_roles |> filter(species %in% present)

  data.frame(
    lake              = lake,
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

write.csv(metaweb_role_presence,
          here("results/metaweb_role_presence_by_lake.csv"),
          row.names = FALSE)

# =============================================================================
# VISUALISATION
# =============================================================================

# --- 13a. SES heatmap ---

ses_long_clean_renamed <- ses_long_clean

ses_long_clean_renamed$lake <- gsub("Lake_", "", ses_long_clean_renamed$lake)

ses_long_clean_renamed$lake <- gsub("_Tarn", " Tarn", ses_long_clean_renamed$lake)

ses_long_clean_renamed$metric <- gsub("connectance", "Connectance", ses_long_clean_renamed$metric)
ses_long_clean_renamed$metric <- gsub("modularity", "Modularity", ses_long_clean_renamed$metric)
ses_long_clean_renamed$metric <- gsub("prop_basal", "% Basal", ses_long_clean_renamed$metric)
ses_long_clean_renamed$metric <- gsub("prop_top", "% Top", ses_long_clean_renamed$metric)
ses_long_clean_renamed$metric <- gsub("trophic_level_mean", "Mean Trophic Level", ses_long_clean_renamed$metric)
ses_long_clean_renamed$metric <- gsub("trophic_level_max", "Max Trophic Level", ses_long_clean_renamed$metric)
ses_long_clean_renamed$metric <- gsub("mean_chain_length", "Mean Chain Length", ses_long_clean_renamed$metric)
ses_long_clean_renamed$metric <- gsub("max_chain_length", "Max Chain Length", ses_long_clean_renamed$metric)


ses_long_clean_renamed$metric <- factor(ses_long_clean_renamed$metric, levels = c("Connectance", "Modularity", "% Basal", 
                   "% Top", "Mean Trophic Level", 
                   "Mean Chain Length", "Max Chain Length"))

p_ses_heatmap <- ses_long_clean_renamed |>
  filter(metric %in% c("Connectance", "Modularity", "% Basal", 
                   "% Top", "Mean Trophic Level",
                   "Mean Chain Length", "Max Chain Length")) |>
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
  filter(metric %in% c("Connectance", "Modularity", "% Basal", 
                   "% Top", "Mean Trophic Level",
                   "Mean Chain Length", "Max Chain Length")) |>
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
       p_ses_envelope, width = 18, height = 5, dpi = 300, bg = "white")

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




df <- read.csv('data/Boxplot_Foodweb metrics_data.csv', check.names = FALSE, stringsAsFactors = FALSE)
df$GroupLabel  <- as.character(df$GroupLabel)
df$LegendGroup <- as.character(df$LegendGroup)

# --- Fill colours used in the plot ------------------------------------------
fill_values <- c(
  'Connectance' = '#365C8DFF',
  'Modularity' = '#365C8DFF',
  'Mean trophic level' = '#365C8DFF',
  'α-diversity ' = '#365C8DFF'
)

p <- ggplot(df, aes(x = GroupLabel, y = Value, fill = LegendGroup)) +
  geom_boxplot(colour = 'black', outlier.colour = 'black', linewidth = 0.6) +
  scale_fill_manual(values = fill_values, name = 'Groups') +
  labs(
    x = 'Groups',
    y = 'Standardised sums of squares'
  )  + 
  theme_classic(base_size = 10) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size=12),
        axis.title.y = element_text(size=14),
        plot.title = element_text(size=15),
        legend.text = element_text(size=14)) + 
  scale_y_continuous(
          limits = c(0, 0.1837),
          minor_breaks = waiver()
        ) +
  theme(legend.position = "")

ggsave(here("results/Boxplot_Foodweb metrics.jpeg"),
       p, width = 10, height = 4, dpi = 300, bg = "white")
