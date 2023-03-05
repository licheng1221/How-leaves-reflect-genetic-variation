### Customized functions ###

# Calculate coefficient of variation
cal_cv <- function(data){
  CV = t(as.data.frame(apply(data,2,sd) / apply(data,2,mean)))
  return(CV)
}

# Output wavelengths that show (adjusted) p-value < 0.05
sigwavs <- function(pvalues){
  p = t(pvalues)
  wavs = rownames(p)[p<0.05]
  return(wavs)
}

# Model selection AZ
best_model_az <- function(wavelength){
  variable = az[,wavelength-350+1]
  M1 <- lm(variable ~ GeneGroup + Ctime + Day, data = az)
  M2 <- lm(variable ~ GeneGroup + Ctime + Day + Leaf_Num, data = az)
  M3 <- lm(variable ~ GeneGroup + Ctime + Batch + Leaf_Num, data = az)
  M4 <- lm(variable ~ GeneGroup + Ctime + Day + Batch + Leaf_Num, data = az)
  return(AIC(M1,M2,M3,M4))
}

# Stepwise selection AZ
step_interact_az <- function(wavelength){
  variable = az[,1877-350+1]
  M = lm(variable ~ GeneGroup + Ctime + Batch + Leaf_Num, data = az)
  new_M = step(M, scope = . ~ .^2, direction = 'forward')
  return(formula(new_M))
}

# Model selection Jena
best_model_jena <- function(wavelength){
  variable = jena[,wavelength-350+1]
  M1 <- lm(variable ~ GeneGroup, data = jena)
  M2 <- lm(variable ~ GeneGroup + Ctime, data = jena)
  return(AIC(M1,M2))
}

