#### Function to estimate genetic offset between provenance and garden 
provgar_offset <- function(RDA, K, env_garden, env_provenance, weights = NULL){
  
  # Combine and scale env data for garden and 
  var_env <- rbind(env_garden, env_provenance)
  var_env <- apply(var_env, 2, scale)
  row.names(var_env) <- c("Garden", row.names(env_provenance))
  
  # Predicting population adaptive index based on RDA axes
  AI <- list()
  for(i in 1:K){
    AI_K <- apply(var_env, 1, function(x) sum(x * RDA$CCA$biplot[,i]))
    names(AI_K) <- row.names(var_env)
    AI[[i]] <- AI_K
    names(AI)[i] <- paste0("RDA", as.character(i))
  }
  
  # Adaptive Indice (AI) along each axis
  taball <- as.data.frame(do.call(cbind, AI))
  
  # Euclidean distance between common garden and provenances RDA scores, weighted or not
  if(!is.null(weights)){
    weights <- (RDA$CCA$eig/sum(RDA$CCA$eig))[1:K]
    taball <- do.call(cbind, lapply(1:K, function(x) taball[,x]*weights[x]))
    eucloffset <- as.matrix(dist(taball, method = "euclidean"))[,1]
  }
  else {
    eucloffset <- as.matrix(dist(taball, method = "euclidean"))[,1]
  }
  
  # Return the euclidean distances 
  return(eucloffset)
}