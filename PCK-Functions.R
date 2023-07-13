require(gstat)
require(proxy)
require(sf)
require(pbapply)
require(utils)

# Calculate Poisson semivariogram
pck.variogram <- function(data, cases = 'Y' , pop = 'N', maxd = NULL, nbins = 15 , pop.rate = 10000){
  
  #Get coordinates
  coords <- st_coordinates(st_centroid(data))/1000
  rownames(coords) <- counties$FIPS
  
  # Calculate number of observations
  data =  data %>% st_drop_geometry()
  No = nrow(data)
  # Select counts and population size
  y = as.numeric(data[,cases])
  n = as.numeric(data[,pop])
  
  
  # Create indexes to use in building of empirical semivariograms
  idx1 = rep(1:(No - 1), times = (No - 1):1)
  idx2 = unlist(sapply(1:(No - 1), function(i) (1:No)[-(1:i)]))
  
  # Reshape distances for unique pairs of points
  d = c(dist(coords, diag = T))
  
  # Calculate difference of unique pairs of points
  diff = (pop.rate * (y[idx1] / n[idx1] - y[idx2] / n[idx2])) ^ 2
  
  # Calculate weighted term
  w = (n[idx1] * n[idx2]) / (n[idx1] + n[idx2])
  # bias term
  bias = sum(y) / sum(n) * pop.rate
  
  # Determine maximum distance of estimator
  if (is.null(maxd)) {
    maxd = max(d) / 2
  }
  # Determine distances of tolerance regions for empirical semivariogram estimator
  bindist = seq(0, maxd, len = nbins + 1)
  # Determine classification cuts for unique pairs of points by distance
  dcuts = cut(d, breaks = bindist)
  
  
  # Calculate population-weighted empirical semivariogram
  np = unlist(lapply(split(d, dcuts), length), use.names = FALSE)
  middle = unlist(lapply(split(d, dcuts), median), use.names = FALSE)
  wbinds = split(w , dcuts)
  rbinds = split((w * diff - bias), dcuts)
  semivariance = unlist(lapply(1:nbins, function(x) {
    sum(rbinds[[x]]) / (2 * sum(wbinds[[x]]))
  }), use.names = FALSE)
  # Plot population-weighted semivariogram
  #plot(1:nbins, semivariance)
  
  # Transform varioram to object for gstat (easy manipulation of variogram)
  var <- data.frame(as.numeric(np) , as.numeric(middle),semivariance,rep(0, nbins),
                    rep(0, nbins),rep(as.factor("var1"), nbins))
  names(var) <-  c("np", "dist", "gamma", "dir.hor", "dir.ver", "id")
  var <- var[!is.na(var$gamma),]
  rownames(var) <- NULL
  class(var) <- c("gstatVariogram", "data.frame")
  
  return(var)
}

# Calculate  Poisson corss-semivariogram
pck.crossvariogram <-function(data, cases.a = 'Ya', cases.b = 'Yb' , cases.ab = 'Yab', pop = 'N', maxd = NULL, nbins = 15 , pop.rate = 10000){
  
  #Get coordinates
  coords <- st_coordinates(st_centroid(data))/1000
  rownames(coords) <- counties$FIPS
  
  # Global parameters
  data =  data %>% st_drop_geometry()
  No = nrow(data)
  
  # Define population sizes
  ya <- as.numeric(data[,cases.a])
  yb <- as.numeric(data[,cases.b])
  yab <- as.numeric(data[,cases.ab])
  # Define process
  n <- as.numeric(data[, pop])
  
  # Define total mean
  Mab <- sum(yab) / sum(n)
  
  # Create indexes to use in building of empirical semivariograms
  idx1 = rep(1:No, times = No:1)
  idx2 = c(1,unlist(sapply(1:No, function(i) (1:No)[-(1:(i - 1))])))
  
  # Reshape distances for unique pairs of points
  dist = as.matrix(dist(coords))
  d = c(dist[lower.tri(dist, diag = T)])
  
  # Calculate difference of unique pairs of points (centered and standardized)
  diff = ( pop.rate * pop.rate *((ya[idx1] / n[idx1]) - (ya[idx2] / n[idx2]))*((yb[idx1] / n[idx1]) - (yb[idx2] / n[idx2])))
  
  # Calculate weighted term
  w = (n[idx1] * n[idx2]) / (n[idx2] + n[idx1])
  # bias term
  bias = Mab * pop.rate
  
  # Determine maximum distance of estimator
  if (is.null(maxd)) {
    maxd = max(d) / 2
  }
  # Determine distances of tolerance regions for empirical semivariogram estimator
  bindist = seq(0, maxd, len = nbins + 1)
  # Determine classification cuts for unique pairs of points by distance
  dcuts = cut(d, breaks = bindist)
  
  
  # Calculate population-weighted empirical semivariogram
  np = unlist(lapply(split(d, dcuts), length), use.names = FALSE)
  middle = unlist(lapply(split(d, dcuts), mean), use.names = FALSE)
  wbinds = split(w , dcuts)
  rbinds = split((w * diff - bias), dcuts)
  semivariance = unlist(lapply(1:nbins, function(x) {
    sum(rbinds[[x]]) / (2 * sum(wbinds[[x]]))
  }), use.names = FALSE)
  
  if(all(semivariance < 0) ){
    stop("semivariance values cannot be negative. Check observations.")
  }
  
  # Transform varioram to object for gstat (easy manipulation of variogram)
  crossvar <- data.frame(as.numeric(np) , as.numeric(middle),semivariance,rep(0, nbins),
                         rep(0, nbins),rep(as.factor("var1.var2"), nbins))
  names(crossvar) <-  c("np", "dist", "gamma", "dir.hor", "dir.ver", "id")
  crossvar <- crossvar[!is.na(crossvar$gamma),]
  rownames(crossvar) <- NULL
  class(crossvar) <- c("gstatVariogram", "data.frame")
  return(crossvar)
}


# ---- Covariance function
Cexp =  function(h,a,b){b * exp(-h / a)}
Csph = function(h,a,b){ifelse(h <= a, b * (1-1.5*(h/a)+0.5*(h/a)^3), 0)}

#--- Select cov model
covariance.model <- function(model){
  if(model == "Exp"){
    cov = Cexp
  }else if(model == "Sph"){
    cov = Csph
  }else{
    stop("incorrect covariance model")
  }
  return(cov)
}


# Adjust linear model of corregionalization
lmc.poisson.cokrige <- function(var.a, var.b, crossvar.ab, data, cases.a = 'Ya', cases.b = 'Yb', pop = 'N', var.params){
  
  # Change id auxilliary variable
  var.b$id = as.factor("var2")
  
  # Integrate variograms
  variograms <- rbind(crossvar.ab, var.b, var.a)
  
  # Get coordinates
  coords <- st_coordinates(st_centroid(data))/1000
  
  # Gstat object
  g <- gstat(id = "var1", formula = c/p ~ 1, loc = ~ x + y, data = data.frame(x = coords[,1] ,y = coords[,2], c = data[,cases.a], p = data[,pop]))
  g <- gstat(g,id = "var2", formula = c/p ~ 1, loc = ~ x + y, data = data.frame(x = coords[,1] ,y = coords[,2], c = data[,cases.b], p = data[,pop]))
  g <- gstat(g, id = c("var1","var2"), model = var.params)
  fitted.lmc <- fit.lmc(variograms, g ,model = var.params, fit.method = 2, fit.ranges = F)
  print(plot(variograms,fitted.lmc))
  return(fitted.lmc)
}


# --- Poisson Cokriging
poisson.cokrige <- function(data.a, data.b, coords.pred, cases.a = 'Ya', cases.b = 'Yb' , cases.ab = 'Yab', pop = 'N', lmc, pop.rate = 10000){
  
  # Get coordinates
  coords.a = st_coordinates(st_centroid(data.a)) / 1000
  coords.b = st_coordinates(st_centroid(data.b)) / 1000
  
  data.a$x = coords.a[, 1]
  data.a$y = coords.a[, 2]
  
  data.b$x = coords.b[, 1]
  data.b$y = coords.b[, 2]
  
  # Remove geometry
  data.a =  data.a %>% st_drop_geometry()
  data.b =  data.b %>% st_drop_geometry()
  
  # Get covariance functions
  cov.funs = list()
  cov.funs$c1 = covariance.model(as.character(lmc$model$var1$model[1]))
  cov.funs$c2 = covariance.model(as.character(lmc$model$var1$model[2]))
  
  #Get covariance parameters
  params = list()
  params$c1 =  list(
    sill.a =  lmc$model$var1$psill[1]           ,
    range.a =  lmc$model$var1$range[1],
    sill.b =  lmc$model$var2$psill[1]           ,
    range.b =  lmc$model$var2$range[1],
    sill.ab =  lmc$model$var1.var2$psill[1]      ,
    range.ab =  lmc$model$var1.var2$range[1]
  )
  params$c2 =  list(
    sill.a =  lmc$model$var1$psill[2]           ,
    range.a =  lmc$model$var1$range[2],
    sill.b =  lmc$model$var2$psill[2]           ,
    range.b =  lmc$model$var2$range[2],
    sill.ab =  lmc$model$var1.var2$psill[2]     ,
    range.ab =  lmc$model$var1.var2$range[2]
  )
  
  if (nrow(coords.a) == nrow(coords.pred)) {
    message("----------------- Smoothing rates ----------------")
    
    
    fold = 1:nrow(data.a)
    preds = data.frame(matrix(as.numeric(NA), nrow(data.a), 6))
    names(preds) = c('rate.pred','rate.var','observed','residual','zscore', 'fold')
    
    folds = sort(unique(fold))
    progress_bar = txtProgressBar(min = 0,max = length(folds),style = 3,char  = "+")
    
    #Validation
    for (i in folds) {
      
      sel = which(fold == i)
      data.cv = data.a[-sel, ]
      data.target = data.a[sel,]
      
      coords.a = cbind(data.cv$x, data.cv$y)
      coords.xo =  cbind(data.target$x, data.target$y)
      
      # distances
      distMa = as.matrix(dist(coords.a))
      distMb = as.matrix(dist(coords.b))
      distMab = as.matrix(proxy::dist(coords.a, coords.b))
      distMba = as.matrix(proxy::dist(coords.b, coords.a))
      
      # data sizes
      size.a = nrow(data.cv)
      size.b = nrow(data.b)
      
      # joint risk
      yab = data.a[, cases.ab]
      
      # ---- Error term W_{lk}
      Waa <- diag((sum(data.cv[, cases.a]) / sum(data.cv[, pop]) * pop.rate) / data.cv[, pop])
      Wbb <- diag((sum(data.b[, cases.b]) / sum(data.b[, pop]) * pop.rate) / data.b[, pop])
      
      indexab <- which(distMab == 0, arr.ind = TRUE)[, 1]
      Wab <- matrix(0, dim(distMab)[1], dim(distMab)[2])
      Wab[which(distMab == 0, arr.ind = TRUE)] <- (sum(yab) / sum(data.cv[, pop]) * pop.rate) / data.cv[, pop][indexab]
      Wba = t(Wab)
      
      # ---- Covariance matrices with W_{lk}
      Caa <- (cov.funs$c1(distMa,  params$c1$range.a,  params$c1$sill.a) + cov.funs$c2(distMa,  params$c2$range.a,  params$c2$sill.a)) + Waa
      Cab <- (cov.funs$c1(distMab,  params$c1$range.ab,  params$c1$sill.ab) + cov.funs$c2(distMab,  params$c2$range.ab,  params$c2$sill.ab))  + Wab
      Cba <- (cov.funs$c1(distMba,  params$c1$range.ab,  params$c1$sill.ab) + cov.funs$c2(distMba,  params$c2$range.ab,  params$c2$sill.ab))  + Wba
      Cbb <- (cov.funs$c1(distMb,  params$c1$range.b,  params$c1$sill.b) + cov.funs$c2(distMb,  params$c2$range.b,  params$c2$sill.b))   + Wbb
      
      # ----- Unbiasedness constrains
      Ctotal <- rbind(cbind(Caa,Cab),cbind(Cba,Cbb))
      Ctotal <- cbind(Ctotal, c(rep(1, (size.a)), rep(0, size.b)), c(rep(0, size.a), rep(1, size.b)))
      Ctotal <- rbind(Ctotal, c(rep(1, size.a), rep(0, size.b), 0, 0), c(rep(0, size.a), rep(1, size.b), 0, 0))
      
      # ---- Inverse matrix
      Cinv <- solve(Ctotal)
      
      #---- Prediction
      pred.target <-
        as.data.frame(t(
          poisson.cokrige.pred(
            coords.xo,
            xta = coords.a,
            xtb = coords.b,
            Cinv = Cinv,
            cov.funs,
            params,
            data.cv,
            data.b,
            cases.a,
            cases.b,
            pop,
            pop.rate
          )
        ))
      
      preds[[1]][sel] = pred.target[[3]]
      preds[[2]][sel] = pred.target[[4]]
      preds[[3]][sel] = data.target[, cases.a] / data.target[, pop] * pop.rate
      preds[[4]][sel] = as.numeric(pred.target[[3]] - preds[[3]][sel])
      preds[[5]][sel] = as.numeric(preds[[4]][sel] / sqrt(pred.target[[4]]))
      preds[[6]][sel] = fold[i]
      
      setTxtProgressBar(progress_bar, value = i)
    }
    
    close(progress_bar)
    
  } else{
    
    message("----------------- Predicting rates ----------------")
    
    
    # distances
    distMa = as.matrix(dist(coords.a))
    distMb = as.matrix(dist(coords.b))
    distMab = as.matrix(proxy::dist(coords.a, coords.b))
    distMba = as.matrix(proxy::dist(coords.b, coords.a))
    
    # data sizes
    size.a = nrow(data.a)
    size.b = nrow(data.b)
    
    # joint risk
    yab = data.a[, cases.ab]
    
    # ---- Error term W_{lk}
    Waa <- diag((sum(data.a[, cases.a]) / sum(data.a[, pop]) * pop.rate) / data.a[, pop])
    Wbb <- diag((sum(data.b[, cases.b]) / sum(data.b[, pop]) * pop.rate) / data.b[, pop])
    
    indexab <- which(distMab == 0, arr.ind = TRUE)[, 1]
    Wab <- matrix(0, dim(distMab)[1], dim(distMab)[2])
    Wab[which(distMab == 0, arr.ind = TRUE)] <- (sum(yab) / sum(data.a[, pop]) * pop.rate) / data.a[, pop][indexab]
    Wba = t(Wab)
    
    # ---- Covariance matrices with W_{lk}
    Caa <- (cov.funs$c1(distMa,  params$c1$range.a,  params$c1$sill.a) + cov.funs$c2(distMa,  params$c2$range.a,  params$c2$sill.a)) + Waa
    Cab <- (cov.funs$c1(distMab,  params$c1$range.ab,  params$c1$sill.ab) + cov.funs$c2(distMab,  params$c2$range.ab,  params$c2$sill.ab))  + Wab
    Cba <- (cov.funs$c1(distMba,  params$c1$range.ab,  params$c1$sill.ab) + cov.funs$c2(distMba,  params$c2$range.ab,  params$c2$sill.ab))  + Wba
    Cbb <- (cov.funs$c1(distMb,  params$c1$range.b,  params$c1$sill.b) + cov.funs$c2(distMb,  params$c2$range.b,  params$c2$sill.b))   + Wbb
    
    # ----- Unbiasedness constrains
    Ctotal <- rbind(cbind(Caa,Cab),cbind(Cba,Cbb))
    Ctotal <- cbind(Ctotal, c(rep(1, (size.a)), rep(0, size.b)), c(rep(0, size.a), rep(1, size.b)))
    Ctotal <- rbind(Ctotal, c(rep(1, size.a), rep(0, size.b), 0, 0), c(rep(0, size.a), rep(1, size.b), 0, 0))
    
    # ---- Inverse matrix
    Cinv <- solve(Ctotal)
    
    #---- Prediction
    # predict
    preds <- as.data.frame(t(pbapply(coords.pred, 1, poisson.cokrige.pred, coords.a, coords.b, Cinv, cov.funs, 
                                     params, data.a, data.b, cases.a, cases.b, pop, pop.rate)))
    
  }
  return(preds)
}


# --- Poisson Cokriging helper function
poisson.cokrige.pred <- function(xao, xta, xtb, Cinv, cov.funs, params, data.cv, data.b, cases.a, cases.b, pop, pop.rate){

  # ----- Predict multiple locations
  ya <- data.cv[, cases.a]
  yb <- data.b[, cases.b]
  na <- data.cv[, pop]
  nb <- data.b[, pop]
  
  distXoa <- unlist(apply(xta,1, distFun, x0 = xao))
  distXob <- unlist(apply(xtb,1, distFun, x0 = xao))
  
  #--- Covariances to prediction
  covaao <-  (cov.funs$c1(distXoa,  params$c1$range.a,  params$c1$sill.a) + cov.funs$c2(distXoa,  params$c2$range.a,  params$c2$sill.a))
  covbao <-  (cov.funs$c1(distXob,  params$c1$range.ab,  params$c1$sill.ab) + cov.funs$c2(distXob,  params$c2$range.ab,  params$c2$sill.ab))
  coTotal <- c(covaao, covbao, 1, 0)
  
  # --- Calculate weights
  lambda <- Cinv%*%coTotal
  
  # predict
  pred.ra  = sum(lambda[1:(length(lambda)-2)] * c(ya/na, yb/nb) * pop.rate)
  var.ra = (cov.funs$c1(0,  params$c1$range.a,  params$c1$sill.a) + cov.funs$c2(0,  params$c2$range.a,  params$c2$sill.a)) - sum(lambda*coTotal)
  
  predictions = c(c(xao[1]), c(xao[2]) , pred.ra, var.ra/10)
  names(predictions) <- c('x','y','rate.pred','rate.var')
  return(predictions)
}

#---- Distance function
distFun <- function(xi, x0){sqrt((x0[1]- xi[1])^2 + (x0[2]- xi[2])^2)}

