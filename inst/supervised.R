
spSupervised <- function(scGroups) {
    tsne <- getData(scGroups, "tsne")
    mclust <- getData(scGroups, "mclust")
    classification <- mclust$classification
    tsne <- as.data.frame(tsne)
    tsne$class <- classification
    
    library(caret)
    inTrain <- createDataPartition(y = tsne$class, p = .75, list = FALSE)
    training <- tsne[inTrain, ]
    testing <- tsne[-inTrain, ]
    testClass <- testing[ ,1:2]
    testing$class <- NULL
    
    set.seed(11)
    
    mod1 <- train(rownames(training) ~ .,
        data = training,
        method = "nnet",
        trControl = trainControl(method = "CV", number = 10)
    )
    
    set.seed(11)

    mod2 <- train(rownames(training) ~ .,
        data = training,
        method = "rbfDDA",
        trControl = trainControl(method = "CV", number = 10)
    )
    
    grid <-  expand.grid(
        hidden_dropout = 0,
        visible_dropout = 0,
        layer1 = 1:5,
        layer2 = 1:3,
        layer3 = 1:3
    )
    
    set.seed(11)
    
    mod3 <- train(rownames(training) ~ .,
        data = training,
        method = "dnn",
        trControl = trainControl(method = "CV", number = 10),
        tuneGrid = grid
    )
    
    mx.set.seed(0)
    
    grid <-  expand.grid(
        hidden_dropout = 0,
        visible_dropout = 0,
        layer1 = 1:5,
        layer2 = 1:5,
        layer3 = 1:5
    )
    
    traMxnet <- data.matrix(training)
    
    mod4 <- train(rownames(training) ~ .,
        data = training,
        method = "mx.mlp",
        trControl = trainControl(method = "CV", number = 10),
        tuneGrid = grid
    )
    
}

#mxnet
#C50 library(C50)
#R CMD javareconf
#C4.5 library(RWeka)

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

##try combining single cells to doublets
sampleType <- getData(expData, "sampleType")
counts <- getData(expData, "counts")
counts.log <- getData(expData, "counts.log")
scGroups <- spUnsupervised(expData)
mclust <- getData(scGroups, "mclust")
classification <- mclust$classification


sngG1 <- counts.log[ ,sampleType=="Singlet" & classification == 2]
sngG2 <- counts.log[ ,sampleType=="Singlet" & classification == 11]

data <- data.frame(
    one = rowMeans(cbind(sng[,1], sng[,2])),
    two = rowMeans(cbind(sng[,3], sng[,4])),
    three = rowMeans(cbind(sng[,5], sng[,6])),
    four = rowMeans(cbind(sng[,7], sng[,8])),
    five = rowMeans(cbind(sng[,9], sng[,10])),
    six = rowMeans(cbind(sng[,11], sng[,12]))
)


counts.log <- getData(expData, "counts.log")
sampleType <- getData(expData, "sampleType")
sng <- counts.log[ ,sampleType == "Singlet"]
dbl <- counts.log[ , sampleType == "Doublet"]


sng <- as.data.frame(t(sng))
sng$class <- as.factor(classification)
sng.hex <- as.h2o(sng)

data = h2o.splitFrame(sng.hex,ratios=c(.7,.15), destination_frames = c("train","test","valid"))
names(data) <- c("Train", "Test", "Valid")

y = "class"
x = names(data$Train)
x = x[-which(x==y)]

m1 = h2o.glm(training_frame = data$Train, validation_frame = data$Valid, x = x, y = y, family='multinomial', solver='L_BFGS')

h2o.confusionMatrix(m1, valid=FALSE) # get confusion matrix in the training data
h2o.confusionMatrix(m1, valid=TRUE)  # get confusion matrix in the validation data

fit = h2o.predict(object = m1, newdata = data$Test)
as.data.frame(fit)$predict
as.data.frame(data$Test)$class

dbl <- as.data.frame(t(dbl))
dbl.hex <- as.h2o(dbl)
fit2 = h2o.predict(object = m1, newdata = dbl.hex)
as.data.frame(fit2)$predict





