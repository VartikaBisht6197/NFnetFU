#Function to plot missing value ggplot
#Input : data frame as matrix 
#Output : 1 and 0 list
missing_val_plot <- function(df){
  new_df <- melt(df)
  miss_val <- c()
  rounded <- round(new_df['value'][[1]], digits = 6)
  for (i in 1:length(rounded)) {
    if( rounded[i] == 0){
      miss_val[i] = 0
    }
    else{
      miss_val[i] = 1
    }
  }
  return(miss_val)
}