plot_resource_histograms <- function(edge_list_high,
                                     edge_list_complete,
                                     save_plot = FALSE,
                                     save_table = FALSE,
                                     outfile_plot = "results/resources.histo.jpeg",
                                     outfile_table = "results/resources_combined.csv") {

  
  #------------------------
  # HIGH
  #------------------------
  resources_high <- edge_list_high |>
    group_by(consumer) |>
    summarise(no_resources = n(), .groups = "drop")
  
  mean_high <- mean(resources_high$no_resources)
  
  p.high <- ggplot(resources_high, aes(x = no_resources)) +
    geom_histogram(binwidth = 5, boundary = 0,
                   closed = "right",
                   fill = "steelblue",
                   color = "black") +
    theme_classic() +
    labs(x = "Number of Resources",
         y = "Count of Consumers") +
    ggtitle("A")
  
  #------------------------
  # COMPLETE
  #------------------------
  resources_complete <- edge_list_complete |>
    group_by(consumer) |>
    summarise(no_resources = n(), .groups = "drop")
  
  mean_complete <- mean(resources_complete$no_resources)
  
  p.complete <- ggplot(resources_complete, aes(x = no_resources)) +
    geom_histogram(binwidth = 5, boundary = 0,
                   closed = "right",
                   fill = "steelblue",
                   color = "black") +
    theme_classic() +
    labs(x = "Number of Resources",
         y = "Count of Consumers") +
    ggtitle("B")
  
  #------------------------
  # COMBINED PLOT
  #------------------------
  p.combo <- cowplot::plot_grid(p.high, p.complete, align = "hv")
  
  #------------------------
  # COMBINED TABLE
  #------------------------
  resources_combined <- full_join(
    resources_high,
    resources_complete,
    by = "consumer",
    suffix = c("_high", "_complete")
  )
  
  #------------------------
  # OPTIONAL SAVES
  #------------------------
  if (save_plot) {
    ggsave(
      filename = outfile_plot,
      plot = p.combo,
      width = 6,
      height = 2,
      units = "in",
      bg = "transparent"
    )
  }
  
  if (save_table) {
    write.csv(resources_combined,
              outfile_table,
              row.names = FALSE)
  }
  
  #------------------------
  # RETURN EVERYTHING
  #------------------------
  return(list(
    plot = p.combo,
    resources_high = resources_high,
    resources_complete = resources_complete,
    resources_combined = resources_combined,
    mean_high = mean_high,
    mean_complete = mean_complete
  ))
}