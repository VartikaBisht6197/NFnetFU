#Initial penalty calculator
#Input : data frame(features) ,data frame(estimate as factor)
#Output : data frame, weights, penalty , initial betas
cal_penalty <- function(data,label_dat){
  
  #Scale the data
  n <- nrow(data)
  data_mean <- colMeans(data)
  data_norm <- sqrt(n-1)*apply(data,2,sd)
  label_dat_mean <- mean(label_dat)
  data <- scale(data, center = data_mean, scale = data_norm)
  
  
  #fit a logistic regression model to the data to calculate the weights
  lm.fit <- glm(label_dat ~ data)
  beta.init <- coef(lm.fit)[-1]
  w  <- abs(beta.init) 

  #Scale the data by 1/weight[]
  data_scale <- scale(data, center=FALSE, scale=1/w)  
  
  #Calculate the initial coefficient
  lasso.fit <- cv.glmnet(data_scale, label_dat, family = "gaussian", alpha = 0, standardize = FALSE, nfolds = 5)
  beta <- predict(lasso.fit, data_scale, type="coefficients", s="lambda.min")[-1]
  beta_initial <- beta
  
  #calculate penalties
  beta <- beta * w
  beta <- matrix(beta, nrow=1)
  data_mean <- matrix(data_mean)
  b0 <- label_dat_mean -  beta%*%data_mean 
  coef <- cbind(b0, beta)
  
  #Name the colums in weight matrix
  names(w) <- colnames(data)
  
  return_val <- list(data,w,coef,beta_initial)
  
  return(return_val)
}
