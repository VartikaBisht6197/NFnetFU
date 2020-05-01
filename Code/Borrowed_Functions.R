#Function to chage a list of float numbers to n decimal places float number
#Input : list , int (decimal places)
#Output : list
in_2_dec <- function(list_given,n){
  j = 1
  for (i in list_given){
    list_given[j] = formatC(i, digits = n, format = "f")
    j = j + 1
  }
  return(list_given)
}

#Create labels for data using factors
#Input : list - num[]
#Output : list - num[]
label_creator <- function(Data_Set_1){
  label_dat <- as.vector(Data_Set_1[2])
  factors <- factor(label_dat$Class)
  label_dat <- as.numeric(factors)
  return(label_dat)}


#Function to create start-end-edge node network data frame
#Input : data frame
#Output : data frame
Network_start_end <- function(df_need,manhattan){
  inDist <- dist(df_need, method = manhattan, diag = FALSE, upper = FALSE, p = 2)
  if (class(inDist) != "dist") stop("wrong input type")
  A <- attr(inDist, "Size")
  B <- if (is.null(attr(inDist, "Labels"))) sequence(A) else attr(inDist, "Labels")
  if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
  if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
  d = data.frame(
    row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
    col = rep(B[-length(B)], (length(B)-1):1),
    value = as.vector(inDist))
  return(d)
}



apply_cosine_similarity <- function(df){
  cos.sim <- function(df, ix) 
  {
    A = df[ix[1],]
    B = df[ix[2],]
    return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
  }   
  n <- nrow(df) 
  cmb <- expand.grid(i=1:n, j=1:n) 
  C <- matrix(apply(cmb,1,function(cmb){ cos.sim(df, cmb) }),n,n)
  C
}

