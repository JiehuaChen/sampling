# Ensemble regression predictions of site-level control yields (Yc) 
# and treatment response ratio indices (SRI) in Malawi.
# Malawi LREP response trial data (courtesy of LREP & Todd Benson)
# LREP data documentation at: https://www.dropbox.com/s/4qbxnz4mdl92pdv/Malawi%20area-specific%20fertilizer%20recs%20report.pdf?dl=0
# Data pre-processing with: https://github.com/mgwalsh/TRM/blob/master/MW_LREP_SI.R
# M.Walsh, J.Chen & A.Verlinden, November 2014

# Required packages
install.packages(c("downloader","raster","rgdal","MASS","rpart","randomForest"),dependencies=TRUE)
require(downloader)
require(raster)
require(rgdal)
require(MASS)
require(rpart)
require(randomForest)

# Data downloads ----------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("Data", showWarnings=F)
dat_dir <- "./Data"

# Site index data download to "./Data"
download("https://www.dropbox.com/s/o9588q2wci8mtiv/MW_Site_Indices.csv?dl=0", "./Data/MW_Site_Indices.csv", mode="wb")
mwsite <- read.table(paste(dat_dir, "/MW_Site_Indices.csv", sep=""), header=T, sep=",")

# Malawi grids download to "./Data" (~7.6 Mb)
download("https://www.dropbox.com/s/54di5f37yp30bz4/MW_grids.zip?dl=0", "./Data/MW_grids.zip", mode="wb")
unzip("./Data/MW_grids.zip", exdir="./Data", overwrite=T)
glist <- list.files(path="./Data", pattern="tif", full.names=T)
mwgrid <- stack(glist)

# Overlay grids & generate dataframes -------------------------------------
coordinates(mwsite) <- ~Easting+Northing
projection(mwsite) <- projection(mwgrid)
exgrid <- extract(mwgrid, mwsite)
Yc <- mwsite$Yc
SRI <- mwsite$SRI
ycdat <- data.frame(cbind(Yc, exgrid))
ycdat <- na.omit(ycdat)
srdat <- data.frame(cbind(SRI, exgrid))
srdat <- na.omit(srdat)

# sorting the data by their fpar values, so that in CV, we can be sure to produce outer-site prediction error
combinedat <- cbind(ycdat, SRI = srdat[, "SRI"])[order(ycdat[,"fPARs"]), ]
Xdat <- ycdat[,-1]

GetPredictions <- function(target, traindata, testX){ 
    # Regression models -------------------------------------------------------
    # Stepwise main effects GLM's
    if(target=="Yc"){
        Yc.glm <- glm(MyTarget ~ ., family=gaussian(link="log"), data=traindata)
        predicts_test <- exp(predict.glm(Yc.glm, testX))
    }
    if(target=="SRI"){
        Yc.glm <- glm(MyTarget ~ .,  data=traindata)
        predicts_test <- predict.glm(Yc.glm, testX)
    }
    # Regression trees
    # Control yield predictions (Yc)
    Yc.rt <- rpart(MyTarget ~ ., data=traindata)
    predicts_test <- cbind(predicts_test, predict(Yc.rt, testX))        

    # Random forests (no tuning default)
    # Control yield predictions (Yc)
    Yc.rf <- randomForest(MyTarget~., data = traindata, importance=T, proximity=T)
    predicts_test <- cbind(predicts_test, predict(Yc.rf, testX))
    return(predicts_test)
}

# number of cross-validation
nf <- 5
nrep <- 3 
num_wgt <- 3
wgt <- matrix(runif(num_wgt*nrep,0,1),nrep,num_wgt) 
wgt <- apply(wgt,1,function(x) x/sum(x))

TARGETS <- c("Yc", "SRI")
target_wgt_error <- new.env()

for(target in TARGETS){
    cat('working on: ', target, '\n')
    wgt_error <- NULL

    for(w in 1:ncol(wgt)){
        weight <- wgt[, w]
        cv_pred_error <- NULL

        for(cv in 1:nf){
            cat('CV',cv,':-------------------------------------\n')
            nd <- floor(nrow(combinedat)/nf)
            testindexes <- ((cv-1)*nd+1):(cv*nd)
            #Retriving train and test data frames from original data set
            trainX <- Xdat[-(testindexes),]
            testX <- Xdat[(testindexes),]
            tstTRG <- combinedat[(testindexes),target]
            MyTarget <- combinedat[-(testindexes), target]
            traindata <- data.frame(MyTarget, trainX)

            predicts_test <- GetPredictions(target, traindata, testX)
            # calculate ensemble predictions
            weighted_pred <-  predicts_test%*%weight
            cv_pred_error <- c(cv_pred_error,sqrt(sum((weighted_pred-tstTRG)^2)/length(nrow(testX))))
            cat('CV error: ', sqrt(sum((weighted.pred-tstTRG)^2)/length(nrow(testX))), '\n')
        }
        wgt_error <- c(wgt_error, mean(cv_pred_error))
    }
    target_wgt_error[[target]] <- wgt_error
}

target_wgt_error <- as.list(target_wgt_error)

# get the weights with minial average CV errors for each target
min_wgt<- wgt[, sapply(target_wgt_error, which.min) ]

# predictons on the whole grid
grid_values <- as.data.frame(getValues((mwgrid)))
for(k in 1:length(names(target_wgt_error))){
    target <- names(target_wgt_error)[k]
    cat("Predicting ", target, "\n")
    
    MyTarget <- combinedat[, target]
    traindata <-  data.frame(MyTarget, Xdat)

    predicts_test <- GetPredictions(target, traindata, grid_values)
    ensemble_predicts <- predicts_test%*%min_wgt[,k]
    ensemble_predicts <- predicts_test%*%min_wgt[,k]
    predict_grid_1k <- SpatialPointsDataFrame(coords = coordinates(mwgrid), data = data.frame(ensemble_predicts = ensemble_predicts))
    gridded(predict_grid_1k)<-TRUE
    writeGDAL (
               dataset=predict_grid_1k["ensemble_predicts"],
               fname=paste("ensemble_",target,".tif", sep=""),
               drivername= "GTiff",
               type="Float32")

}
    

