build_local_webs <- function(pa_table, metaweb, label) {
  mw_spp     <- V(metaweb)$name
  local_webs <- list()
  
  for (lake in rownames(pa_table)) {
    present            <- colnames(pa_table)[pa_table[lake, ] == 1]
    spp_in_metaweb     <- intersect(present, mw_spp)
    local_webs[[lake]] <- induced_subgraph(metaweb, vids = spp_in_metaweb)
  }
  
  cat("\nLocal webs built for:", label, "\n")
  map_dfr(names(local_webs), function(l) {
    data.frame(lake = l,
               S    = vcount(local_webs[[l]]),
               L    = ecount(local_webs[[l]]))
  }) |> print()
  
  local_webs
}
