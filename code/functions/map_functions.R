
plot_west_coast_map <- function(sd.plankton,
                                nz_bathy,
                                wc_bathy,
                                outfile = "results/West.Coast.map.jpeg",
                                width = 6,
                                height = 6) {
  
  # Base bathy plot (West Coast)
  West.Coast.p <- autoplot.bathy(
    wc_bathy,
    geom = c("contour", "raster"),
    colour = "white",
    size = 0.1,
    show.legend = FALSE
  ) +
    scale_fill_gradient2(low = "#FFFFFF", mid = "#FFFFFF", high = "#a19f9c") +
    theme_void()
  
  # Main map
  West.Coast.p1 <- West.Coast.p +
    geom_point(
      data = sd.plankton,
      aes(x = Longitude, y = Latitude),
      colour = "#029cbd",
      size = 2
    ) +
    geom_text_repel(
      data = sd.plankton,
      aes(x = Longitude, y = Latitude, label = Lake),
      size = 2.5,
      segment.size = 0.1
    ) +
    north(
      x.min = 168, x.max = 172.2,
      y.min = -42.3, y.max = -45,
      scale = 0.1, symbol = 3,
      anchor = c(x = 170.0, y = -42.35)
    ) +
    ggsn::scalebar(
      sd.plankton,
      dist = 25, dist_unit = "km",
      transform = TRUE,
      model = "WGS84",
      height = 0.01,
      st.size = 3,
      border.size = 0.5,
      anchor = c(x = 170.8, y = -42.5)
    ) +
    xlab("Longitude (°)") +
    ylab("Latitude (°)") +
    theme_bw()
  
  # Inset NZ map
  NZ.p <- autoplot.bathy(
    nz_bathy,
    geom = c("contour", "raster"),
    colour = "white",
    size = 0.1,
    show.legend = FALSE
  ) +
    scale_fill_gradient2(low = "#FFFFFF", mid = "#FFFFFF", high = "#a19f9c") +
    theme_void()
  
  West.Coast.p2 <- West.Coast.p1 +
    annotation_custom(
      grob = ggplotGrob(NZ.p),
      xmin = 168.2, xmax = 169.4,
      ymin = -42.3, ymax = -43.5
    ) +
    annotate(
      "rect",
      xmin = 168.2, xmax = 169.4,
      ymin = -42.37, ymax = -43.4,
      alpha = 0,
      color = "black"
    )
  
  ggsave(
    filename = outfile,
    plot = West.Coast.p2,
    width = width,
    height = height,
    units = "in",
    bg = "transparent"
  )
  
  return(West.Coast.p2)
}