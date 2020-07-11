## Input : Numeric Labels and Microbiome Abundance Data

#ANFIS preprocessing of data
data_train <- cbind(data1,label_dat)
range.data.input <- matrix(0,nrow = 2,ncol = dim(data_train)[2])
for (i in 1:dim(data_train)[2]) {
  range.data.input[1,i] <- min(data_train[,i])
  range.data.input[2,i] <- max(data_train[,i])
}

#Correlation plot for original dataset
tiff("Correlation plot for original dataset.tiff", width = 10, height = 10, units = 'in', res = 300)
par(cex = 0.7)
corrplot(cor(data1), type = "upper")
dev.off()

#ANFIS
microbiome_ANFIS <- frbs.learn(data_train, range.data = range.data.input, method.type = c("ANFIS"),control = list())
print("ANFIS DONE!")

#Rules Matrix Dimesions
rules_mat <- microbiome_ANFIS[["rule.data.num"]]
colnames(rules_mat) <- microbiome_ANFIS[["colnames.var"]]
rules_r <- dim(rules_mat)[1]
rules_c <- dim(rules_mat)[2]

#Save Labels
label_dat <- scale(rules_mat[,rules_c])
colnames(label_dat) <- colnames(rules_mat)[rules_c]
write.csv(label_dat, "Labels After ANFIS.csv")
print("New labels have been assigned!")

#Save Integer format of rules
rules_int <- rules_mat[,1:(rules_c-1)]
write.csv(rules_int,"Ruled Based Matrix.csv")
print("Rule based matrix is saved!")

#Correlation plot for rule based dataset
tiff("Correlation plot for rule based dataset.tiff", width = 10, height = 10, units = 'in', res = 300)
par(cex = 0.7)
corrplot(cor(rules_int), type = "upper")
dev.off()

#Significantly Correlated : No color = (p value cut off 0.05)
tiff("p value Correlation plot for rule based dataset.tiff", width = 10, height = 10, units = 'in', res = 300)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
par(cex = 0.4)
p.mat <- rcorr(rules_int, type ="spearman")
corrplot(p.mat$r, method = "color", col = col(200),number.cex = .7,
         type = "upper", addCoef.col = "black",tl.col = "black",
         p.mat = p.mat$P, sig.level = 0.05, insig = "blank", tl.srt = 90, diag = TRUE)
dev.off()

#Zero Standard Deviation
sdev_0 <- c()
for (i in 1:dim(rules_int)[2]) {
  if(sd(rules_int[,i]) < 0.1){
    sdev_0 <- c(sdev_0,i)
  }
}
if(length(sdev_0)>0)
{
  print(sprintf("%s features had Zero-Standard-Deviation",length(sdev_0)))
  write.csv(colnames(rules_int)[sdev_0], "Zero-Standard-Deviation features.csv")
  print("Zero-Standard-Deviation features saved.")
  
  #New Rule Based Matrix without Zero-Standard-Deviation Features
  rules_int <- subset(rules_int, select = -sdev_0)
  rules_c <- dim(rules_int)[2]
  write.csv(rules_int, "Ruled Based Matrix without Zero-Standard-Deviation Features.csv")
  print("New Ruled Based Matrix without Zero-Standard-Deviation Features saved")
  
  } else {
  sprintf("%s features had Zero-Standard-Deviation",length(sdev_0)) }

#Scale and Save the Integer format of rules
scaled_rules_int <- scale(rules_int)
write.csv(scaled_rules_int, "Scaled Ruled Based Matrix.csv")
print("Scaled Ruled Based Matrix saved")

## Output : Rule Based Matrix , Scaled Rule Based Matrix and Labels

