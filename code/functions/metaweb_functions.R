# =============================================================================
# Build Edge list
# =============================================================================

build_edge_list <- function(master_taxa, 
                            all_consumers = FALSE,
                            test = FALSE
){
  
  #decide if the edge list must be based on all possible interactions (AI) or just the ones that have been reviewed
  if(all_consumers){
    consumers <- rlang::sym("all_consumers")
  }else{
    consumers <- rlang::sym("highest_potential_consumers")
  }
  
  
  results <- list(
    "taxa_top" = NULL,
    "taxa_bottom" = NULL,
    "cannibals" = NULL,
    "mutual_predation"= NULL,
    "edge_list" = NULL
  )
  
  top_consumers <- master_taxa |>
    filter(is.na(!!consumers) | !!consumers== "" | !!consumers == " ")
  
  #check the taxa we deem as top consumers.
  if (nrow(top_consumers) > 0) {
    results[[1]] <- top_consumers
  }
  
  #check that all taxa have at least one resource unless they are an autotroph. 
  
  master_taxa_resources <- master_taxa |>
    separate_rows(!!consumers, sep = ";\\s*") |>
    distinct(taxa_name, !!consumers, .keep_all = TRUE)|>
    group_by(taxa_name) |>
    mutate("preference"= row_number()) |>
    ungroup()
  
  master_taxa_resources <- master_taxa_resources |>
    group_by(!!consumers) |>
    summarise(resources = paste0(taxa_name,collapse = ";"),
              preferences = paste0(preference, collapse = ";")) |>
    rename("taxa_name"= consumers)
  
  master_taxa <- left_join(master_taxa,
                           master_taxa_resources,
                           by="taxa_name")
  
  #check the ones that are not autotrophs and do not consume any other taxa
  bottom_taxa <- filter(master_taxa, is.na(resources)) |>
    filter(feeding !=  "autotroph" | feeding != "filter-feeder" | feeding != "decomposer" | feeding != "heterotroph" | feeding != "detritivore" | feeding != "phytoplankton")
    #filter(feeding !=  "autotroph")
  
  if(nrow(bottom_taxa) >0){
    warning(
      "There are taxa that do not consume any other taxa in the metaweb and are not autotrophs.They will be dropped.Check bottom_taxa.")
    
    master_taxa <- dplyr::filter(master_taxa,!taxa_name %in% bottom_taxa$taxa_name )
    
    results[[2]] <- bottom_taxa
    
  }
  
  #building edge list
  
  
  #this list is clean and does not have any isolated nodes.
  # First, we transform lack of resources into a character string to avoid problems in string matching
  #  Then, we string match all the taxa from the resource categories. This is why resource categories has to be reviewed in DETAIL.
  
  edge_list <- master_taxa |>
    separate_rows(resources, preferences, sep = ";\\s*") |>
    rename("consumer"= taxa_name,
           "resource"= resources) |>
    select(consumer, resource, preferences) |>
    mutate(preferences = as.numeric(preferences))
  
  #Check if there is any cannibalistic links or intraguild predation
  
  edge_list <- edge_list |>
    mutate("cannibalism"= case_when(consumer == resource ~ 1,
                         TRUE ~0),
           "intraguild_1"= paste0(consumer,"_", resource),
           "intraguild_2"= paste0(resource,"_", consumer),
           "mutual_predation"= case_when(intraguild_1 %in% intraguild_2 ~ 1,
                                         TRUE ~0))
  
  cannibals <- filter(edge_list, cannibalism ==1)
  mutual_predation <- filter(edge_list, mutual_predation== 1)
  
  if(nrow(cannibals) >0){
    warning("There are cannibalistic links present. They will be dropped. Check cannibals.")
    
    edge_list <- dplyr::filter(edge_list, cannibalism ==0 )
    
    results[[3]] <- cannibals
    
  }
  
  if(nrow(mutual_predation) >0){
    warning("There are mutual predation links present.")
    
    #edge_list <- dplyr::filter(edge_list, mutual_predation==0 )
    
    results[[4]] <- mutual_predation
    
  }
  
  edge_list <- select(edge_list,
                      consumer,
                      resource)

  
  if(test==TRUE){
    Frac_obs <- .8
    n_obs = round(nrow(edge_list)*Frac_obs,0)
    all_obs <- seq(1:nrow(edge_list))
    keep_obs <- sample(all_obs, size=n_obs,replace=FALSE)
    edge_list <- edge_list[keep_obs,]
    
  }
  
  
  results[[5]] <- edge_list
  
  return(results)
  
}

# =============================================================================
# Build Metaweb
# =============================================================================

build_metaweb <- function(edge_list,
                          master_taxa,
                          plot= TRUE,
                          remove_isolated_nodes=TRUE,
                          probabilistic = FALSE){
  
  results <- list("isolated_nodes"= NULL, "cannibalism" =NULL, "adjacency_matrix"=NULL, "community" = NULL)
  
  nodes <- master_taxa |>
    rename("node" = taxa_name) |>
    select(node, feeding) #####Added feeding
  
  rownames(nodes) <- nodes$node
  
  
  if(probabilistic){
    trophic.links <- edge_list |>
      select(consumer,
             resource,
             probability) |>
      filter(!is.na(resource),
             probability > 0)
  }else{
    trophic.links <- edge_list |>
      select(consumer,
             resource) |>
      filter(!is.na(resource))
  }
  # we need to drop autotrophs to get it working
  
  properties <- list("title" = "West Coast Lakes Metaweb", "habitat" = "Freshwater")
  
  
  metaweb_community <- cheddar::Community(nodes = nodes,
                                          trophic.links = trophic.links,
                                          properties = properties)

  #check for isolated nodes
  iso <- IsolatedNodes(metaweb_community)
  
  if(length(iso) >0 ){
    results[[1]] <- iso
    if(remove_isolated_nodes){
      metaweb_community <- cheddar::RemoveIsolatedNodes(community=metaweb_community)
    }
  }
  
  #check for cannibalism
  cannibalism <- cheddar::Cannibals(metaweb_community)
  
  if(length(cannibalism)>0){
    results[[2]] <-cannibalism
    metaweb_community <- cheddar::RemoveCannibalisticLinks(community = metaweb_community)
  }
  
  #create matrix
  if(probabilistic){
    adjacency_matrix <-cheddar::PredationMatrix(metaweb_community, weight="probability" )
  }else{
    adjacency_matrix <-cheddar::PredationMatrix(metaweb_community)
  }

  
  results[[3]] <- adjacency_matrix
  results[[4]] <- metaweb_community
  
  cheddar::PlotWebByLevel(metaweb_community)
  
  return(results)
  
}







# =============================================================================
# IGRAPH -> CHEDDAR CONVERSION
# =============================================================================

igraph_to_cheddar <- function(g, community_name = "community") {
  
  nodes <- data.frame(node = V(g)$name)
  
  el <- as_edgelist(g, names = TRUE) |>
    as.data.frame() |>
    # IMPORTANT: adjust this if your edges go resource -> consumer
    # Check: a known top predator should have high outdegree in igraph
    # and should appear in TopLevelNodes() in cheddar after conversion
    setNames(c("resource", "consumer"))
  
  Community(
    nodes         = nodes,
    trophic.links = el,
    properties    = list(title = community_name)
  )
}

# =============================================================================
# SPECIES MATCHING DIAGNOSTICS
# =============================================================================

diagnose_matching <- function(pa_table_t, metaweb, label) {
  
  pa_spp <- colnames(pa_table_t)
  mw_spp <- V(metaweb)$name
  
  cat("\n--- Matching diagnostics:", label, "---\n")
  cat("Species in PA table:    ", length(pa_spp), "\n")
  cat("Species in metaweb:     ", length(mw_spp), "\n")
  cat("Matched (intersection): ", length(intersect(pa_spp, mw_spp)), "\n")
  cat("In PA but NOT metaweb:  ", length(setdiff(pa_spp, mw_spp)), "\n")
  cat("In metaweb but NOT PA:  ", length(setdiff(mw_spp, pa_spp)), "\n")
  
  missing_pa <- setdiff(pa_spp, mw_spp)
  missing_mw <- setdiff(mw_spp, pa_spp)
  
  list(
    summary = map_dfr(rownames(pa_table), function(lake) {
      present    <- pa_spp[pa_table[lake, ] == 1]
      in_metaweb <- intersect(present, mw_spp)
      
      data.frame(
        lake         = lake,
        S_detected   = length(present),
        S_in_metaweb = length(in_metaweb),
        S_missing    = length(present) - length(in_metaweb),
        pct_matched  = round(100 * length(in_metaweb) / max(length(present), 1), 1),
        metaweb      = label
      )
    }),
    missing_in_metaweb = missing_pa,
    missing_in_pa      = missing_mw
  )
}


# =============================================================================
# NETWORK METRICS FUNCTION (igraph + cheddar)
# =============================================================================

# Non-redundant metrics only:
#   DROPPED: link_density, mean_degree, sd_degree, mean_indegree, mean_outdegree
#   (all are linear functions of L/S — identical SES to connectance)
#   ADDED via cheddar: trophic levels, omnivory, chain length, food web roles

compute_network_metrics <- function(g, community_name = "lake") {
  
  S <- vcount(g)
  L <- ecount(g)
  
  empty_row <- data.frame(
    S                  = S,
    L                  = L,
    connectance        = NA,
    link_density = NA,
    clustering         = NA,
    modularity         = NA,
    prop_basal         = NA,
    prop_top           = NA,
    prop_int           = NA,
    prop_omnivory      = NA,
    trophic_level_mean = NA,
    trophic_level_max  = NA,
    trophic_level_sd   = NA,
    max_chain_length   = NA,
    mean_chain_length  = NA
  )
  
  if (S < 3 || L == 0) return(empty_row)
  
  # --- igraph metrics ---
  con   <- L / (S * (S - 1))
  link_density <- L / S 
  clust <- transitivity(g, type = "global")
  
  g_un <- as_undirected(g, mode = "collapse")
  comm  <- tryCatch(cluster_louvain(g_un), error = function(e) NULL)
  mod   <- if (!is.null(comm) && length(unique(membership(comm))) > 1)
    modularity(comm) else NA
  
  max_chain_length <- diameter(g, directed = TRUE) #max chain length
  mean_chain_length <- mean_distance(
    g,
    directed = TRUE,
    unconnected = TRUE
  ) #mean chain length
  
  # --- cheddar metrics ---
  ch <- tryCatch({
    
    ch_comm <- igraph_to_cheddar(g, community_name)
    tl      <- PreyAveragedTrophicLevel(ch_comm)
    
    data.frame(
      prop_basal         = FractionBasalNodes(ch_comm),
      prop_top           = FractionTopLevelNodes(ch_comm),
      prop_int           = FractionIntermediateNodes(ch_comm),
      prop_omnivory      = FractionOmnivorous(ch_comm),
      trophic_level_mean = mean(tl, na.rm = TRUE),
      trophic_level_max  = max(tl,  na.rm = TRUE),
      trophic_level_sd   = sd(tl,   na.rm = TRUE)
    )
    
  }, error = function(e) {
    message("cheddar failed for: ", community_name, " | ", e$message)
    data.frame(
      prop_basal = NA, prop_top = NA, prop_int = NA,
      prop_omnivory = NA, trophic_level_mean = NA,
      trophic_level_max = NA, trophic_level_sd = NA
    )
  })
  
  data.frame(S = S, L = L, connectance = con, link_density = link_density, 
             clustering = clust, modularity = mod, max_chain_length = max_chain_length, mean_chain_length = mean_chain_length) |>
    bind_cols(ch)
}
