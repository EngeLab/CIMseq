#H2O multiclass

library(h2o)
h2o.init(nthreads=-1, max_mem_size="2G")
h2o.removeAll()
D = h2o.importFile(path = normalizePath("covtype.full.csv"))
h2o.summary(D)

data = h2o.splitFrame(D,ratios=c(.7,.15),destination_frames = c("train","test","valid"))
names(data) <- c("Train","Test","Valid")

y = "Cover_Type"
x = names(data$Train)
x = x[-which(x==y)]

m1 = h2o.glm(training_frame = data$Train, validation_frame = data$Valid, x = x, y = y,family='multinomial',solver='L_BFGS')
h2o.confusionMatrix(m1, valid=TRUE)


m2 = h2o.glm(training_frame = data$Train, validation_frame = data$Valid, x = x, y = y,family='multinomial',solver='L_BFGS', lambda = 0)
h2o.confusionMatrix(m2, valid=FALSE) # get confusion matrix in the training data
h2o.confusionMatrix(m2, valid=TRUE)  # get confusion matrix in the validation data





library('sp.scRNAseq')
library(h2o)
h2o.init(nthreads=1, max_mem_size="7G")
h2o.removeAll()

##try combining single cells to doublets
sampleType <- getData(expData, "sampleType")
counts <- getData(expData, "counts")
counts.log <- getData(expData, "counts.log")
scGroups <- spUnsupervised(expData)
mclust <- getData(scGroups, "mclust")
classification <- mclust$classification

dbl <- as.data.frame(counts[ ,sampleType=="Doublet"])
sng <- as.data.frame(counts[ ,sampleType=="Singlet"])
sngG1 <- sng[ ,classification == 2]
sngG2 <- sng[ ,classification == 11]

combined <- data.frame(
one = rowMeans(cbind(sngG1[,1], sngG2[,1])),
two = rowMeans(cbind(sngG1[,2], sngG2[,2])),
three = rowMeans(cbind(sngG1[,3], sngG2[,3])),
four = rowMeans(cbind(sngG1[,4], sngG2[,4])),
five = rowMeans(cbind(sngG1[,5], sngG2[,5])),
six = rowMeans(cbind(sngG1[,6], sngG2[,6])),
seven = rowMeans(cbind(sngG1[,7], sngG2[,7])),
eight = rowMeans(cbind(sngG1[,8], sngG2[,8])),
nine = rowMeans(cbind(sngG1[,9], sngG2[,9])),
ten = rowMeans(cbind(sngG1[,10], sngG2[,10])),
eleven = rowMeans(cbind(sngG1[,11], sngG2[,12])),
twelve = rowMeans(cbind(sngG1[,12], sngG2[,12])),
thirteen = rowMeans(cbind(sngG1[,13], sngG2[,13])),
row.names = rownames(sng)
)

sng <- as.data.frame(t(sng))
combined <- as.data.frame(t(combined))
sng$class <- as.factor(classification)
sng.hex <- as.h2o(sng)
combined.hex <- as.h2o(combined)

data = h2o.splitFrame(sng.hex,ratios=c(.75), destination_frames = c("train","valid"), seed=11)
names(data) <- c("Train", "Valid")

y = "class"
x = names(data$Train)
x = x[-which(x==y)]

m1 = h2o.gbm(
training_frame = data$Train,
validation_frame = data$Valid,
x = x,
y = y,
ntrees = 500,
min_rows = 10,
learn_rate = 0.001,
balance_classes=TRUE,
seed=11,
distribution="multinomial",
keep_cross_validation_predictions = TRUE,
keep_cross_validation_fold_assignment = TRUE,
score_each_iteration = TRUE,
nfolds=10
)

m1@model$training_metrics
h2o.confusionMatrix(m1, valid=FALSE) # get confusion matrix in the training data
h2o.confusionMatrix(m1, valid=TRUE)  # get confusion matrix in the validation data

fit = h2o.predict(object = m1, newdata = data$Test)
pre <- as.data.frame(fit)
real <- as.data.frame(data$Test)$class
pre$real <- real
pred <- pre[ ,c('real', 'predict', paste('p',1:12, sep=""))]

library(ggplot2)
library(reshape2)
m <- melt(pred, id.vars=c("real", "predict"))
ggplot(m, aes(x=predict, y=value, fill=variable))+
geom_bar(stat="identity", position=position_dodge(width=1))+
facet_grid(~real, scales="free")

fit2 <- h2o.predict(object = m1, newdata = combined.hex )
pre2 <- as.data.frame(fit2)
pre2$n <- 1:nrow(pre2)
m <- melt(pre2, id.vars=c("n", "predict"))

ggplot(m, aes(x=factor(n), y=value, fill=variable))+
geom_bar(stat="identity", position=position_dodge(width=1))+
facet_grid(. ~ n, scales="free")
