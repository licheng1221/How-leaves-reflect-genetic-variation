### Customized functions ###

# LOF outliers detector
detect_outliers <- function(rawdata, spec_start_col,k) {
  # Define the types
  types <- c("WR", "WRL", "BR","BRL")
  # Sepc data
  df <- rawdata[,spec_start_col:(spec_start_col+2150)]
  # Initialize a vector to hold the filenames of the outliers
  outlier_filenames <- c()
  # Loop over the types
  for (type in types) {
    print(paste("Processing:", type))
    # Filter the data for the specific type
    rawdata_type <- rawdata[rawdata$type == type,]
    df_type <- df[rawdata$type == type,]
    # Calculate LOF scores
    lof_scores <- DescTools::LOF(df_type, k = k)
    # Identify outliers
    outlier_lof <- lof_scores > 2
    # Get the outlier data
    outliers <- rawdata_type[outlier_lof,]
    
    # Check if there are any outliers
    if (nrow(outliers) > 0) {
      # Get the unique Plant_IDs of the outliers
      outlier_plants <- unique(outliers$Plant_ID)
      # Check the conditions for removing plants
      for (Plant_ID in outlier_plants) {
        # Get the outlier data for this plant
        plant_outliers <- outliers[outliers$Plant_ID == Plant_ID,]
        # Get the filename and Plant_ID of the outlier
        filename <- plant_outliers$filename
        # Extract the number(s) before ".asd" in the filename
        number <- as.numeric(stringr::str_extract(filename, "(\\d{5})(?=\\.asd)"))
        # If there are more than one outlier scans, remove the plant
        if(nrow(plant_outliers) > 1){
          if(nrow(plant_outliers) == 2 & number[1] %% 5 == 0){
            print(paste("Outlier scan detected in", type, ": filename =", filename[2], ", Plant_ID =", Plant_ID,". Not first scan, set to NAs"))
                  # Set the values of the outlier to NAs
                  rawdata[rawdata$filename == filename[2],spec_start_col:(spec_start_col+2150)] <- NA
          } else{
          rawdata <- rawdata[rawdata$Plant_ID != Plant_ID,]
          print(paste("Removed plant", Plant_ID, "from", type, "due to multiple outlier scans:", filename))}
        } else {
          # Check if the outlier is the first scan
          if (number %% 5 == 0) {
            print(paste("Outlier scan detected in", type, ": filename =", filename, ", Plant_ID =", Plant_ID,". First scan, will be removed"))
          } else {
            print(paste("Outlier scan detected in", type, ": filename =", filename, ", Plant_ID =", Plant_ID,". Not first scan, set to NAs"))
            # Set the values of the outlier to NAs
            rawdata[rawdata$filename == filename,spec_start_col:(spec_start_col+2150)] <- NA
            # Add the filename to the outlier_filenames vector
            #outlier_filenames <- c(outlier_filenames, filename)
          }
        }
      }
    } else {
      print(paste("No outliers detected in", type))
    }
  }
  #return(list(rawdata = rawdata, outlier_filenames = outlier_filenames))
  return(rawdata)
}

# Calculated Reflectance and Absolute uncertainty 
calculate_CR_AU <- function(df){
  AU = matrix(, nrow = nrow(df)/20, 2151)
  CR = matrix(, nrow = nrow(df)/20, 2151)
  n = 1
  for (i in seq(1,nrow(df),20)){
    for (j in 1:ncol(df)){
      WR = mean(df[(i+1):(i+4),j],na.rm=TRUE)
      WRL = mean(df[(i+6):(i+9),j],na.rm=TRUE)
      BR = mean(df[(i+11):(i+14),j],na.rm=TRUE)
      BRL = mean(df[(i+16):(i+19),j],na.rm=TRUE)
      STD_WR = sd(df[(i+1):(i+4),j],na.rm=TRUE)
      STD_WRL = sd(df[(i+6):(i+9),j],na.rm=TRUE)
      STD_BR = sd(df[(i+11):(i+14),j],na.rm=TRUE)
      STD_BRL = sd(df[(i+16):(i+19),j],na.rm=TRUE)
      CR[n,j] = (BRL*WR - WRL*BR)/(WR-BR)
      AU[n,j] = sqrt((BR*(WRL-BRL)/(WR-BR)^2)^2 * (STD_WR/2)^2 + (BR/(WR-BR))^2 * (STD_WRL/2)^2 + (WR*(WRL-BRL)/(WR-BR)^2)^2 * (STD_BR/2)^2 + (WR/(WR-BR))^2 * (STD_BRL/2)^2)
    }
    n = n+1
  }
  return(list("CR"= CR, "AU" = AU))
}

# Calculate coefficient of variation
cal_cv <- function(data){
  CV = t(as.data.frame(apply(data,2,sd) / apply(data,2,mean)))
  return(CV)
}

# Output wavelengths that show (adjusted) p-value < 0.05
sigwavs <- function(pvalues){
  p = t(pvalues)
  wavs = as.numeric(rownames(p)[p<0.05])
  # If wavs is empty, return "None"
  if(length(wavs) == 0) {
    return("None")
  }
  # Find continuous ranges
  ranges = list()
  start = wavs[1]
  for (i in 2:length(wavs)) {
    if (wavs[i] != wavs[i-1] + 1) {
      ranges = c(ranges, paste(start, "~", wavs[i-1]))
      start = wavs[i]
    }
  }
  ranges = c(ranges, paste(start, "~", wavs[length(wavs)]))
  return(ranges)
}

# Convert a hexadecimal color to RGB, and make it transparent
RGB <- function(hex_color){
  rgb_values <- col2rgb(hex_color)
  rgb_normalized <- as.numeric(rgb_values) / 255
  RGB <- rgb(rgb_normalized[1], rgb_normalized[2], rgb_normalized[3], alpha = 0.5)
  return(RGB)
}

# Model selection AZ
best_model_az <- function(wavelength){
  variable = az[,wavelength-350+1]
  M1 <- lm(variable ~ GeneGroup + Ctime + Day, data = az)
  M2 <- lm(variable ~ GeneGroup + Ctime + Day + Leaf_Num, data = az)
  M3 <- lm(variable ~ GeneGroup + Ctime + Batch + Leaf_Num, data = az)
  M4 <- lm(variable ~ GeneGroup + Ctime + Day + Batch + Leaf_Num, data = az)
  #M5 <- lmer(variable ~ Ctime + Day + Batch + Leaf_Num + (1|GeneGroup), data = az)
  return(AIC(M1,M2,M3,M4))
}

# Stepwise selection AZ
step_interact_az <- function(wavelength){
  variable = az[,wavelength-350+1]
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


