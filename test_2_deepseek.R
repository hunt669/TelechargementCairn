# Nettoyer l'environnement et libérer la mémoire
rm(list = ls()); gc()

# Charger les bibliothèques
library(dplyr)
library(tidyr)
library(influential)
library(fastnet)
library(igraph)
library(doParallel)
library(parallel)
library(progress)

# Définir les fonctions
min_max_norm <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

clustered_seeding <- function(seeds, g, seeds_needed) {
  possible_seeds <- unique(unlist(sapply(seeds, function(x) neighbors(g, x, mode = "total"))))
  seeds_to_add <- possible_seeds[!possible_seeds %in% seeds]
  need_seeds <- seeds_needed > length(seeds_to_add)
  if (need_seeds) {
    return(possible_seeds)
  } else {
    final_seeds <- c(seeds, sample(seeds_to_add, seeds_needed))
    return(final_seeds)
  }
}

get_simple_centralities <- function(g) {
  centrality_df <- data.frame(
    seed = as.numeric(V(g)),
    degree = as.numeric(degree(g)),
    betweenness = as.numeric(betweenness(g)),
    eigen = as.numeric(eigen_centrality(g)$vector)
  )
  percolation <- collective.influence(graph = g, vertices = V(g), mode = "all", d = 3)
  centrality_df$percolation <- as.numeric(percolation)
  return(centrality_df)
}

get_complex <- function(seed, N, g, gmat, thresholds, num_seeds_to_add) {
  gmat_simulation <- matrix(nrow = N, ncol = N, 0)
  num_seeds_i <- num_seeds_to_add[seed]
  
  if (num_seeds_i > 0) {
    seeds <- as.numeric(neighbors(g, seed, mode = "total")) 
    num_seeds_needed <- num_seeds_i - length(seeds)
    need_seeds <- num_seeds_needed > 0 
    
    if (need_seeds) {
      seeds <- clustered_seeding(seeds, g, num_seeds_needed)
    } else if (length(seeds) > 1) {
      seeds <- c(seed, sample(seeds, num_seeds_i))
    } else {
      seeds <- c(seed, seeds)
    }
  } else {
    seeds <- seed
  }
  
  activated <- logical(N)
  activated[seeds] <- TRUE
  gmat_simulation[seeds, ] <- 1
  gmat_simulation[, seeds] <- 1
  gmat_simulation_run <- gmat * gmat_simulation
  
  spread <- TRUE
  while (spread) {
    influence <- colSums(gmat_simulation_run)
    influence_activated <- influence >= thresholds
    t <- which(activated)
    t_1 <- which(influence_activated)
    t_full <- union(t, t_1)
    spread <- length(t_full) > length(t)
    activated[t_full] <- TRUE
    adopters <- which(activated)
    gmat_simulation[adopters, ] <- 1
    gmat_simulation[, adopters] <- 1
    gmat_simulation_run <- gmat * gmat_simulation
  }
  
  num_adopters <- sum(activated)
  complex_g <- graph_from_adjacency_matrix(gmat_simulation_run)
  all_distances <- distances(complex_g, seed, V(complex_g), mode = "out")
  all_distances[all_distances == Inf] <- NA
  PLci <- mean(all_distances, na.rm = TRUE)
  
  return(c(seed, N, num_seeds_i, num_adopters, PLci))
}

# Charger et préparer les données
N <- 3739
g <- read.csv("C:/R_test_/r_test_1.csv", sep = ";")
g <- as.matrix(g) 
g <- graph_from_adjacency_matrix(g)
gmat <- as.matrix(as_adjacency_matrix(g))

# Configurer les paramètres de seuil
T_dist <- "homo"
T_type <- "abs"
thresholds <- replicate(N, 3)

if (T_type == "frac") {
  thresholds <- round(sapply(1:N, function(x) {
    thresholds[x] * length(neighbors(g, x, mode = "total"))
  }))
}

num_seeds_to_add <- thresholds - 1
num_seeds_to_add[num_seeds_to_add < 0] <- 0
thresholds[thresholds <= 0] <- 1

# Calculer les centralités simples
cat("Calcul des centralités simples...\n")
simple_centralities_df <- get_simple_centralities(g)
cat("Centralités simples calculées!\n")

# Initialiser le suivi de progression
pb <- progress_bar$new(
  format = "[:bar] :percent | Temps écoulé: :elapsed | Temps restant: :eta",
  total = N, 
  clear = FALSE, 
  width = 80)

# Initialiser le chronomètre
start_time <- Sys.time()

# Exécuter le modèle avec suivi de progression
model_output_list <- vector("list", N)

for (x in 1:N) {
  model_output_list[[x]] <- get_complex(x, N, g, gmat, thresholds, num_seeds_to_add)
  pb$tick()
  
  # Sauvegarde intermédiaire toutes les 100 itérations
  if (x %% 100 == 0 || x == N) {
    # Construire le dataframe temporaire
    temp_df <- as.data.frame(do.call(rbind, model_output_list[1:x]))
    colnames(temp_df) <- c("seed", "N", "num_neigh_seeds", "num_adopters", "PLci")
    temp_df$PLci_norm <- min_max_norm(temp_df$PLci)
    
    # Fusionner avec les centralités
    temp_full <- merge(temp_df, simple_centralities_df, by = "seed")
    
    # Sauvegarder
    saveRDS(temp_full, paste0("temp_results_", x, ".rds"))
    cat("\nSauvegarde intermédiaire à l'itération", x, "\n")
  }
}

# Post-traitement final
model_output_df <- as.data.frame(do.call(rbind, model_output_list))
colnames(model_output_df) <- c("seed", "N", "num_neigh_seeds", "num_adopters", "PLci")
model_output_df$PLci_norm <- min_max_norm(model_output_df$PLci)
model_output_full <- merge(model_output_df, simple_centralities_df, by = "seed")

# Sauvegarde finale
saveRDS(model_output_full, "final_results.rds")
write.csv(model_output_full, "final_results.csv", row.names = FALSE)

# Calcul du temps total
end_time <- Sys.time()
execution_time <- difftime(end_time, start_time, units = "mins")

# Rapport final
cat("\n\nCALCUL TERMINÉ!\n")
cat("Temps d'exécution total:", round(execution_time, 1), "minutes\n")
cat("Fichiers de résultats:\n")
cat("- final_results.csv (format CSV)\n")
cat("- final_results.rds (format R optimisé)\n")
cat("Aperçu des résultats:\n")
print(head(model_output_full))

# Notification sonore
if (Sys.info()['sysname'] == "Windows") {
  tryCatch({
    for (i in 1:3) {
      system("rundll32 user32.dll,MessageBeep -1")
      Sys.sleep(0.5)
    }
  }, error = function(e) NULL)
} else {
  tryCatch({
    system("tput bel")
  }, error = function(e) NULL)
}
