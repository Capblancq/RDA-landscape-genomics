#### Function to predict genomic offset from a RDA model
genomic_offset <- function(RDA, K, env_pres, env_fut, range = NULL, method = "loadings", scale_env, center_env){
  
  # Mask with the range if supplied
  if(!is.null(range)){
    env_pres <- mask(env_pres, range)
    env_fut <- mask(env_fut, range)
  }
  
  # Formatting and scaling environmental rasters for projection
  var_env_proj_pres <- as.data.frame(scale(rasterToPoints(env_pres[[row.names(RDA$CCA$biplot)]])[,-c(1,2)], center_env[row.names(RDA$CCA$biplot)], scale_env[row.names(RDA$CCA$biplot)]))
  var_env_proj_fut <- as.data.frame(scale(rasterToPoints(env_fut[[row.names(RDA$CCA$biplot)]])[,-c(1,2)], center_env[row.names(RDA$CCA$biplot)], scale_env[row.names(RDA$CCA$biplot)]))

  # Predicting pixels genetic component based on the loadings of the variables
  if(method == "loadings"){
    # Projection for each RDA axis
    Proj_pres <- list()
    Proj_fut <- list()
    Proj_offset <- list()
    for(i in 1:K){
      # Current climates
      ras_pres <- env_pres[[1]]
      ras_pres[!is.na(ras_pres)] <- as.vector(apply(var_env_proj_pres[,names(RDA$CCA$biplot[,i])], 1, function(x) sum( x * RDA$CCA$biplot[,i])))
      names(ras_pres) <- paste0("RDA_pres_", as.character(i))
      Proj_pres[[i]] <- ras_pres
      names(Proj_pres)[i] <- paste0("RDA", as.character(i))
      # Future climates
      ras_fut <- env_fut[[1]]
      ras_fut[!is.na(ras_fut)] <- as.vector(apply(var_env_proj_fut[,names(RDA$CCA$biplot[,i])], 1, function(x) sum( x * RDA$CCA$biplot[,i])))
      Proj_fut[[i]] <- ras_fut
      names(ras_fut) <- paste0("RDA_fut_", as.character(i))
      names(Proj_fut)[i] <- paste0("RDA", as.character(i))
      # Single axis genetic offset 
      Proj_offset[[i]] <- abs(Proj_pres[[i]] - Proj_fut[[i]])
      names(Proj_offset)[i] <- paste0("RDA", as.character(i))
    }
  }
  
  # Predicting pixels genetic component based on predict.RDA
  if(method == "predict"){ 
    # Prediction with the RDA model and both set of envionments 
    pred_pres <- predict(RDA, var_env_proj_pres[,-c(1,2)], type = "lc")
    pred_fut <- predict(RDA, var_env_proj_fut[,-c(1,2)], type = "lc")
    # List format
    Proj_offset <- list()    
    Proj_pres <- list()
    Proj_fut <- list()
    for(i in 1:K){
      # Current climates
      ras_pres <- rasterFromXYZ(data.frame(var_env_proj_pres[,c(1,2)], Z = as.vector(pred_pres[,i])), crs = crs(env_pres))
      names(ras_pres) <- paste0("RDA_pres_", as.character(i))
      Proj_pres[[i]] <- ras_pres
      names(Proj_pres)[i] <- paste0("RDA", as.character(i))
      # Future climates
      ras_fut <- rasterFromXYZ(data.frame(var_env_proj_pres[,c(1,2)], Z = as.vector(pred_fut[,i])), crs = crs(env_pres))
      names(ras_fut) <- paste0("RDA_fut_", as.character(i))
      Proj_fut[[i]] <- ras_fut
      names(Proj_fut)[i] <- paste0("RDA", as.character(i))
      # Single axis genetic offset 
      Proj_offset[[i]] <- abs(Proj_pres[[i]] - Proj_fut[[i]])
      names(Proj_offset)[i] <- paste0("RDA", as.character(i))
    }
  }
  
  # Weights based on axis eigen values
  weights <- RDA$CCA$eig/sum(RDA$CCA$eig)
  
  # Weighing the current and future adaptive indices based on the eigen values of the associated axes
  Proj_offset_pres <- do.call(cbind, lapply(1:K, function(x) rasterToPoints(Proj_pres[[x]])[,-c(1,2)]))
  Proj_offset_pres <- as.data.frame(do.call(cbind, lapply(1:K, function(x) Proj_offset_pres[,x]*weights[x])))
  Proj_offset_fut <- do.call(cbind, lapply(1:K, function(x) rasterToPoints(Proj_fut[[x]])[,-c(1,2)]))
  Proj_offset_fut <- as.data.frame(do.call(cbind, lapply(1:K, function(x) Proj_offset_fut[,x]*weights[x])))
  
  # Predict a global genetic offset, incorporating the K first axes weighted by their eigen values
  ras <- Proj_offset[[1]]
  ras[!is.na(ras)] <- unlist(lapply(1:nrow(Proj_offset_pres), function(x) dist(rbind(Proj_offset_pres[x,], Proj_offset_fut[x,]), method = "euclidean")))
  names(ras) <- "Global_offset"
  Proj_offset_global <- ras
  
  # Return projections for current and future climates for each RDA axis, prediction of genetic offset for each RDA axis and a global genetic offset 
  return(list(Proj_pres = Proj_pres, Proj_fut = Proj_fut, Proj_offset = Proj_offset, Proj_offset_global = Proj_offset_global, weights = weights[1:K]))
}
