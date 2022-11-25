###################################
#            RONA                 #
###################################

go_rona <- function(Y, X, X.pred, snps.set){
  
  nb_var <- ncol(X)
  n <- nrow(Y)

  mod_lm <- lm(as.matrix(Y[,snps.set]) ~ ., data = data.frame(X)) 
  sm <- summary(mod_lm)
  B <- sapply(sm, FUN = function(x) x$coeff[1:(nb_var + 1), 1])
  
  print(dim(B))
  
  
  X <-  cbind(rep(1.0, n), X)
  X.pred <- cbind(rep(1.0, n), X.pred)
  print(dim(X))
  

  Y.fit <- as.matrix(X) %*% as.matrix(B)
  Y.pred <- as.matrix(X.pred) %*% as.matrix(B)
  
  allele_frequency_shift <- abs(Y.fit - Y.pred)
  
  return(rowMeans(allele_frequency_shift))
  
}

go_rona_conf <- function(X, X.pred, lfmm.obj, snps.set){
  B = lfmm.obj@B
  
  Y.fit <- as.matrix(X) %*% t(B)
  Y.pred <- as.matrix(X.pred) %*% t(B)
  
  allele_frequency_shift <- abs(Y.fit - Y.pred)
  
  return(rowMeans(allele_frequency_shift[,snps.set]))
}

go_rona_pc <- function(Y, X, X.pred, pc, snps.set){
  
  d <- ncol(X)
  lm.obj <- lm(Y ~ cbind(X, pc))
  smr <- summary(lm.obj)
  coeff <- c()
  for (sm in smr){
    coeff <- cbind(coeff, sm$coefficients[,1])
  }
  
  B <- t(coeff[2:(d+1),])

  
  Y.fit <- as.matrix(X) %*% t(B)
  Y.pred <- as.matrix(X.pred) %*% t(B)
  
  allele_frequency_shift <- abs(Y.fit - Y.pred)
  
  return(rowMeans(allele_frequency_shift[,snps.set]))
}


###################################
#              RDA                #
###################################


go_rda <- function(Y, X, X.pred, snps.set){
  
  nb_var <- ncol(X)
  n <- nrow(Y)

  mod_lm <- lm(as.matrix(Y[,snps.set]) ~ ., data = data.frame(X)) 
  sm <- summary(mod_lm)
  B <- sapply(sm, FUN = function(x) x$coeff[1:(nb_var + 1), 1])
  
 
  X <-  as.matrix(cbind(rep(1.0, n), X))
  X.pred <-  as.matrix(cbind(rep(1.0, n), X.pred))

  
  pc = prcomp(X %*% B)
  
  proj.x = predict(pc, X %*% B)
  proj.xpred = predict(pc, X.pred %*% B)
  
  rda.go = rowSums((proj.x - proj.xpred)[,1:ncol(X)]^2)/ncol(Y[,snps.set])
  
  return(rda.go)
}

go_rda_conf <- function(Y, X, X.pred, lfmm.obj, snps.set){
  
  nb_var <- ncol(X)
  n <- nrow(Y)
  
  U = lfmm.obj@U
  V = lfmm.obj@V
  B = lfmm.obj@B
  
  X <- as.matrix(X)
  X.pred <- as.matrix(X.pred)
  

  Y.fit <- X %*% t(B) + U %*% t(V)
  Y.pred <- X.pred %*% t(B) + U %*% t(V)

  Y.fit <- Y.fit[,snps.set]
  Y.pred <- Y.pred[,snps.set]
  pc = prcomp(Y.fit)
  
  # calcul de l'offset sur les ncol(X) + ncol(U) =  4 + K dimensions libres
  
  proj.x = predict(pc, Y.fit)
  proj.xpred = predict(pc, Y.pred)
  
  rda.go = rowSums((proj.x - proj.xpred)[,1:(ncol(X)+ncol(U))]^2)/length(snps.set)
  
  return(rda.go)
}

go_rda_pc <- function(Y, X, X.pred, pc_value, snps.set){
  
  d <- ncol(X)

  lm.obj <- lm(Y ~ cbind(X, pc_value))
  smr <- summary(lm.obj)
  coeff <- c()
  for (sm in smr){
    coeff <- cbind(coeff, sm$coefficients[,1])
  }
  
  B <- t(coeff[2:(d+1),])
  

  X <- as.matrix(X)
  X.pred <- as.matrix(X.pred)
  
  
  Y.fit <- X %*% t(B) 
  Y.pred <- X.pred %*% t(B) 
  
  Y.fit <- Y.fit[,snps.set]
  Y.pred <- Y.pred[,snps.set]
  pc = prcomp(Y.fit)
  
  # calcul de l'offset sur les ncol(X) + ncol(U) =  4 + K dimensions libres
  
  proj.x = predict(pc, Y.fit)
  proj.xpred = predict(pc, Y.pred)
  
  rda.go = rowSums((proj.x - proj.xpred)[,1:(ncol(X)+ncol(pc_value))]^2)/length(snps.set)
  
  return(rda.go)
}


###################################
#              GF                 #
###################################




run_gf <- function(Y, X, causal_set, confounding_var=c()){
  # GF requires a dataframe object
  # We create two lists containing names for predictor (X) and names for output (OUT)
  # These lists will be used as column names in the data frame
  nb_env_var = dim(X)[2]
  var_name <- sapply(seq(1,nb_env_var), function(x) paste("X",x, sep=""))
  nb_causal <- length(causal_set)
  output_name <- sapply(seq(1,nb_causal), function(x) paste("OUT",x, sep=""))
  
  # This block of code check whether the user has specified variable as confounding var
  # If it's the case we will add these to the column names.
  nb_confound <- dim(confounding_var)[2]
  if (is.null(nb_confound)){
    confound_name <- c()
    df_gf <- data.frame(Y[,causal_set], X)
  }else{
    confound_name <- sapply(seq(1,nb_confound), function(x) paste("U",x, sep=""))
    df_gf <- data.frame(Y[,causal_set], X, confounding_var)
  }
  
  # In GF the predictor are the environmental variables and the possible confounding variables that are specified
  # The difference between environmental variables and confounding variables is that confounding variables will not
  # change in the new environmental condition
  pred_name <- c(var_name, confound_name)
  colnames(df_gf) <- c(output_name, pred_name)
  
  # This command will allow us to obtain a gradient forest object
  gf <- gradientForest(data=df_gf, predictor.vars=pred_name, response.vars = output_name, ntree=500)
  return(list(gf=gf, pred_name = var_name))
}

gf_pred <- function(gf, X, X_pred, pred_name){
  
  # We create data frame for current and new environmental conditions
  df_cur_var <- data.frame(X)
  df_fut_var <- data.frame(X_pred)
  colnames(df_cur_var) <- pred_name
  colnames(df_fut_var) <- pred_name
  
  # We obtain the cumulative importance value from GF for current and new environmental conditions
  currentcumimp <- predict(gf, df_cur_var)
  futurcumimp <- predict(gf, df_fut_var)
  
  # We compute the euclidean distance between current and new cumulative importance
  nb_ind <- nrow(X)
  genetic_offset <- c()
  for (i in seq(1,nb_ind)){
    genetic_offset <- c(genetic_offset, l2norm(futurcumimp[i,], currentcumimp[i,]))
  }
  
  return(genetic_offset)
}

go_gf <- function(Y, X, X_pred, causal_set, confounding_var = NULL){
  
  
  gf <- run_gf(Y, X, causal_set, confounding_var)
  genoffset <- gf_pred(gf$gf, X, X_pred, gf$pred_name)
  
  
  return(list(go=genoffset, varimp=gf$gf$overall.imp))
}


###############################
#   CLIMATE DISTANCE          #
###############################

causal_distance <- function(X.causal.cur, X.causal.pred){
  nb_var <- nrow(X.causal.cur)
  
  causal_distance <- c()
  for (i in seq(1,nb_var)){
    optimum_var1 <- X.causal.pred[i,1]
    phenotype_var1 <- X.causal.cur[i,1]
    optimum_var2 <- X.causal.pred[i,2]
    phenotype_var2 <- X.causal.cur[i,2]
    sigma_K <- 0.5
    causal_distance <- c(causal_distance, (1/2)*((optimum_var1 - phenotype_var1)^2/sigma_K^2 + (optimum_var2 - phenotype_var2)^2/sigma_K^2))
  }
  
  return(causal_distance)
}

euclidean_distance <- function(X, X.pred){
  nb_var <- nrow(X)
  
  euclidean.distance <- c()
  for (i in seq(1,nb_var)){
    euclidean.distance <- c(euclidean.distance, dist(rbind(X[i,], X.pred[i,]), method="euclidean"))
  }
  
  return(euclidean.distance)
}
