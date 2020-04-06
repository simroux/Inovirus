# library(glmnet)
library(randomForest)
library(methods)

args <- commandArgs(trailingOnly=TRUE)
in_file <- args[1]
out_file <- args[2]
# search()
file_model <- "InovirusModel.RData"
if (args[3] != ""){
	file_model <- args[3]
}
load(file=file_model) 



# in_file <- "Tmp_input.csv"
# out_file <- "Tmp_output.csv"
print(in_file)
print(out_file)

data<-read.csv(in_file)
summary(data)
result <- data.frame(data[,1:2])
## for glm
# x<-data[,c(-1,-2)]
# summary(as.matrix(x))
# print(chosen_lambda)
# result$pred <-predict(glmmod,newx=as.matrix(x),s=chosen_lambda,type="response")
## for rf
pred_rf<-predict(fit, data, type="prob")
result$pred<-as.vector(pred_rf[,2])

write.csv(as.matrix(result),file=out_file)
