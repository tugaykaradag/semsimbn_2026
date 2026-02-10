# ============================================================
# 5 CONNECTED Random DAGs with 3???5 latent variables (subset of A..E)
# - Reproducible via seed
# - Allows 3-node DAGs like A->B, B->C, A->C
# - Ensures the chosen node-set is connected (as undirected)
# - No SEM / no manifest / no aggregation
# - Plots all 5 DAGs
# ============================================================

library(bnlearn)

# ---------- helper: connectivity check (treat DAG as undirected) ----------
is_connected_undirected <- function(dag) {
  am <- amat(dag)
  und <- (am + t(am)) > 0
  p <- nrow(und)
  if (p <= 1) return(TRUE)
  
  seen <- rep(FALSE, p)
  q <- 1
  seen[1] <- TRUE
  
  while (length(q) > 0) {
    v <- q[1]; q <- q[-1]
    neigh <- which(und[v, ])
    new <- neigh[!seen[neigh]]
    if (length(new) > 0) {
      seen[new] <- TRUE
      q <- c(q, new)
    }
  }
  all(seen)
}

# ---------- connected random DAG generator on a given node set ----------
make_connected_random_dag_on_nodes <- function(nodes,
                                               max_parents = 2,
                                               extra_edge_prob = 0.35,
                                               seed = NULL,
                                               max_tries = 5000) {
  if (!is.null(seed)) set.seed(seed)
  
  for (tt in seq_len(max_tries)) {
    ord <- sample(nodes, length(nodes), replace = FALSE)
    
    # Backbone chain guarantees connectivity on this node-set:
    arcs_backbone <- cbind(from = ord[-length(ord)], to = ord[-1])
    
    parent_count <- setNames(rep(0, length(nodes)), nodes)
    for (k in seq_len(nrow(arcs_backbone))) {
      parent_count[arcs_backbone[k, "to"]] <- parent_count[arcs_backbone[k, "to"]] + 1
    }
    
    arcs_extra <- matrix(nrow = 0, ncol = 2, dimnames = list(NULL, c("from","to")))
    
    for (j in 2:length(ord)) {
      child <- ord[j]
      possible_parents <- ord[1:(j - 1)]
      
      cand <- possible_parents[runif(length(possible_parents)) < extra_edge_prob]
      
      # avoid duplicating the backbone parent
      cand <- setdiff(cand, ord[j - 1])
      
      room <- max_parents - parent_count[child]
      if (room <= 0) next
      if (length(cand) > room) cand <- sample(cand, room)
      
      if (length(cand) > 0) {
        arcs_extra <- rbind(arcs_extra, cbind(from = cand, to = rep(child, length(cand))))
        parent_count[child] <- parent_count[child] + length(cand)
      }
    }
    
    arcs_all <- unique(rbind(arcs_backbone, arcs_extra))
    
    dag <- empty.graph(nodes)
    arcs(dag) <- arcs_all
    
    if (is_connected_undirected(dag)) {
      return(list(dag = dag, order = ord))
    }
  }
  
  stop("Could not generate a connected DAG within max_tries. Try increasing max_tries or extra_edge_prob.")
}

# ---------- main: generate 5 DAGs with random size 3..5 from A..E ----------
generate_dags_3to5 <- function(n_dags = 5,
                               node_pool = c("A","B","C","D","E"),
                               size_range = 3:5,
                               max_parents = 2,
                               extra_edge_prob = 0.35,
                               seed_master = 20260204) {
  set.seed(seed_master)
  
  dags <- vector("list", n_dags)
  
  for (k in seq_len(n_dags)) {
    p <- sample(size_range, 1)                 # 3,4, or 5
    nodes_k <- sort(sample(node_pool, p))      # choose subset
    
    # per-DAG seed so it???s reproducible and independent
    seed_k <- seed_master + 1000 + k
    
    dags[[k]] <- make_connected_random_dag_on_nodes(
      nodes = nodes_k,
      max_parents = max_parents,
      extra_edge_prob = extra_edge_prob,
      seed = seed_k
    )
  }
  
  names(dags) <- paste0("DAG_", seq_len(n_dags))
  dags
}

dags <- generate_dags_3to5(
  n_dags = 5,
  node_pool = c("A","B","C","D","E"),
  size_range = 3:5,
  max_parents = 2,
  extra_edge_prob = 0.35,
  seed_master = 20260204
)

# ---------- print each DAG's node-set + arcs ----------
for (nm in names(dags)) {
  cat("\n", nm, " | nodes:", paste(nodes(dags[[nm]]$dag), collapse = ","), "\n", sep = "")
  print(arcs(dags[[nm]]$dag))
}

# ---------- plot ----------
plot_dags <- function(dags_list) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  par(mfrow = c(2, 3), mar = c(1,1,3,1))
  can_graphviz <- requireNamespace("Rgraphviz", quietly = TRUE)
  
  for (nm in names(dags_list)) {
    dag <- dags_list[[nm]]$dag
    if (can_graphviz) {
      bnlearn::graphviz.plot(dag, main = paste0(nm, " (p=", length(nodes(dag)), ")"))
    } else {
      if (!requireNamespace("igraph", quietly = TRUE)) {
        plot.new()
        title(main = paste0(nm, "\nInstall Rgraphviz or igraph"))
        next
      }
      ig <- igraph::graph_from_data_frame(bnlearn::arcs(dag), directed = TRUE, vertices = bnlearn::nodes(dag))
      plot(ig, main = paste0(nm, " (p=", length(nodes(dag)), ")"))
    }
  }
}

plot_dags(dags)
