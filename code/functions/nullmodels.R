run_null_model <- function(metaweb, local_webs, label, n_null = 999) {
  
  cat("\nRunning null model for:", label, "(n =", n_null, ")\n")
  all_spp <- V(metaweb)$name
  
  with_progress({
    p <- progressor(steps = length(local_webs))
    
    null_results <- map(names(local_webs), function(lake) {
      
      p(message = lake)
      S_obs <- vcount(local_webs[[lake]])
      if (S_obs < 3) return(NULL)
      
      future_map_dfr(1:n_null, function(i) {
        rand_spp <- sample(all_spp, size = S_obs, replace = FALSE)
        rand_web <- induced_subgraph(metaweb, vids = rand_spp)
        compute_network_metrics(rand_web, community_name = paste0("null_", i))
      }, .options = furrr_options(seed = TRUE)) |>
        mutate(iteration = 1:n_null)
      
    }) |> set_names(names(local_webs))
  })
  
  null_results
}



compute_ses <- function(observed_metrics, null_metrics, metaweb_label, metrics_cols) {
  
  map_dfr(observed_metrics$lake, function(lake) {
    
    obs  <- observed_metrics[observed_metrics$lake == lake, ]
    null <- null_metrics[[lake]]
    
    if (is.null(null) || nrow(obs) == 0) return(NULL)
    
    row <- data.frame(lake = lake, metaweb = metaweb_label,
                      S = obs$S, L = obs$L)
    
    for (m in metrics_cols) {
      obs_val   <- obs[[m]]
      null_vals <- null[[m]]
      
      # Track how many non-NA null values we actually have
      n_valid   <- sum(!is.na(null_vals))
      null_mean <- mean(null_vals, na.rm = TRUE)
      null_sd   <- sd(null_vals,   na.rm = TRUE)
      
      row[[paste0("obs_",        m)]] <- obs_val
      row[[paste0("null_mean_",  m)]] <- null_mean
      row[[paste0("null_sd_",    m)]] <- null_sd
      row[[paste0("null_n_valid_",m)]]<- n_valid
      row[[paste0("SES_",        m)]] <- if (!is.na(null_sd) && null_sd > 0 && n_valid >= 100)
        (obs_val - null_mean) / null_sd else NA
      # Only compute p-value if enough valid null draws
      row[[paste0("p_upper_",    m)]] <- if (n_valid >= 100)
        mean(null_vals >= obs_val, na.rm = TRUE) else NA
      row[[paste0("p_lower_",    m)]] <- if (n_valid >= 100)
        mean(null_vals <= obs_val, na.rm = TRUE) else NA
    }
    
    row
  })
}




plot_null_dist <- function(null_list, observed_metrics, metric,
                           metaweb_label = "HCL") {
  
  metric_labels <- c(
    connectance        = "Connectance",
    modularity         = "Modularity",
    prop_basal         = "% Basal",
    prop_top           = "% Top",
    trophic_level_mean = "Mean Trophic Level",
    trophic_level_max  = "Max Trophic Level",
    mean_chain_length  = "Mean Chain Length",
    max_chain_length   = "Max Chain Length"
  )
  
  plot_data <- map_dfr(names(null_list), function(lake) {
    if (is.null(null_list[[lake]])) return(NULL)
    null_list[[lake]] |>
      select(all_of(metric)) |>
      mutate(lake = gsub("Lake_|_Tarn", " Tarn", gsub("Lake_", "", lake)))
  })
  
  obs_data <- observed_metrics |>
    select(lake, all_of(metric)) |>
    rename(obs_val = all_of(metric)) |>
    mutate(lake = gsub("_Tarn", " Tarn", gsub("Lake_", "", lake)))
  
  ggplot(plot_data, aes(x = .data[[metric]])) +
    geom_histogram(fill = "#175c7d", colour = "white",
                   bins = 40, alpha = 0.7) +
    geom_vline(data = obs_data,
               aes(xintercept = obs_val),
               colour = "#c22c43", linewidth = 0.8, linetype = "dashed") +
    facet_wrap(~lake, scales = "free") +
    theme_classic(base_size = 14) +
    theme(strip.background = element_blank(),
          strip.text        = element_text(size = 14)) +
    labs(x = metric_labels[[metric]], y = "Count")
}




plot_null_dist <- function(null_list, observed_metrics, metric,
                           metaweb_label = "AC") {
  
  metric_labels <- c(
    connectance        = "Connectance",
    modularity         = "Modularity",
    prop_basal         = "% Basal",
    prop_top           = "% Top",
    trophic_level_mean = "Mean Trophic Level",
    trophic_level_max  = "Max Trophic Level",
    mean_chain_length  = "Mean Chain Length",
    max_chain_length   = "Max Chain Length"
  )
  
  plot_data <- map_dfr(names(null_list), function(lake) {
    if (is.null(null_list[[lake]])) return(NULL)
    null_list[[lake]] |>
      select(all_of(metric)) |>
      mutate(lake = gsub("Lake_|_Tarn", " Tarn", gsub("Lake_", "", lake)))
  })
  
  obs_data <- observed_metrics |>
    select(lake, all_of(metric)) |>
    rename(obs_val = all_of(metric)) |>
    mutate(lake = gsub("_Tarn", " Tarn", gsub("Lake_", "", lake)))
  
  ggplot(plot_data, aes(x = .data[[metric]])) +
    geom_histogram(fill = "#175c7d", colour = "white",
                   bins = 40, alpha = 0.7) +
    geom_vline(data = obs_data,
               aes(xintercept = obs_val),
               colour = "#c22c43", linewidth = 0.8, linetype = "dashed") +
    facet_wrap(~lake, scales = "free") +
    theme_classic(base_size = 12) +
    theme(strip.background = element_blank(),
          strip.text        = element_text(size = 12)) +
    labs(x = metric_labels[[metric]], y = "Count",
         title = paste(metaweb_label, "—", metric_labels[[metric]]))
}