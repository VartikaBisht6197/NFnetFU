data_train = cbind(data1,label_dat)

#Learn the data to give rules
microbiome_ANFIS <- frbs.learn(data_train, range.data = NULL, method.type = c("ANFIS"),control = list())

#Save Integer format of rules
rules_int <- microbiome_ANFIS[["rule.data.num"]]
rules_int <- rules_int[ 1:dim(rules_int)[1] , 1:(dim(rules_int)[2]-1) ]
colnames(rules_int) <- microbiome_ANFIS[["colnames.var"]][1:dim(rules_int)[2]]
write.csv(rules_int, "rules_int.csv")

#Scale and Save the Integer format of rules
scaled_rules_int <- scale(microbiome_ANFIS[["rule.data.num"]])
scaled_rules_int <- scaled_rules_int[ 1:dim(rules_int)[1] , 1:dim(rules_int)[2] ]
colnames(scaled_rules_int) <- colnames(rules_int)
write.csv(scaled_rules_int, "scaled_rules_int.csv")

#Correlation plot for integer rules
jpeg("Correlation plot for integer rules.jpg", width = 1000, height = 1000)
corrplot(cor(rules_int), type = "upper")
dev.off()
#Correlation plot for original Data
jpeg("Correlation plot for original Data.jpg", width = 1000, height = 1000)
corrplot(cor(data1), type = "upper")
dev.off()

#Check smilarity/Correlation between new "rule-based" data and original data
cormat_TSK <- round(cor(rules_int),2)
cormat_data <- round(cor(data1),2)
print(paste0("Similarity between new rule-based data and original data : ", format( cor(c(cormat_TSK), c(cormat_data)), digits = 3, format = 'f' ) ) )
