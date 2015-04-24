#!/usr/bin/env Rscript

## Running the Random Forest to classify transcripts between coding and non coding
## Need the randomForest and ROCR library

#########################
##### PREPROCESSING #####
#########################
## Include the tools package to get file with no extension
library(tools)

## Checking the number of arguments
args <- commandArgs(TRUE)

## Set generic values
codFile  <- args[1]
nonFile  <- args[2]
testFile <- args[3]
outFile  <- args[4]
outVar   <- paste(file_path_sans_ext(outFile), "_varImpPlot.png", sep="")
outStats <- paste(file_path_sans_ext(outFile), "_stats.txt", sep="")

## If number of argument == 4 then the threshold is not set, need to run 10-cross fold-validation
## and the packages randomForest and ROCR is needed
if(length(args) == 4)
    {
        thres            <- NULL
        outROC           <- paste(file_path_sans_ext(outFile), "_ROC.png", sep="")
        list.of.packages <- c("ROCR","randomForest")
    } else if(length(args) == 5)
## If number of argument == 5 then the threshold is set
## and only the library randomForest is needed
    {
        thres            <- args[5]
        list.of.packages <- c("randomForest")
    } else
## Else the number of arguments is not good
    {
        cat("Error: the number of argument pass to codpot_randomforest.r is wrong:\nUsing 10-fold cross-validation:\tcodpot_randomforest.r inCoding inNonCoding testFile outFile\nUsing a pre-defined threshold:\tcodpot_randomforest.r inCoding inNonCoding testFile outFile threshold\nQuit\n")
        quit()
    }


## Checking for the install of these packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) != 0)
    {
        cat("Please wait during the installation of the R packages (only done once): ", new.packages, ".\n", sep="")
        install.packages(new.packages, repos="http://cran.r-project.org/",  dependencies = TRUE)
        cat("R packages: ", new.packages, " installed.\n", sep="")
        }
## Loading library
for(pack in list.of.packages)
    {
        suppressMessages(library(pack, quietly=TRUE, verbose=FALSE, character.only=TRUE))
    }

###########################
##### CODE START HERE #####
###########################

## Read files
codMat  <- read.table(file=codFile,  header=TRUE, sep="\t")
nonMat  <- read.table(file=nonFile,  header=TRUE, sep="\t")
testMat <- read.table(file=testFile, header=TRUE, sep="\t")

## Modification of codMat and nonMat to add a label value (0: non coding; 1: coding)
codMat <- cbind(codMat, label=rep(x=1, length.out=nrow(codMat)))
nonMat <- cbind(nonMat, label=rep(x=0, length.out=nrow(nonMat)))

## Define the matrix -dat- with codMat and nonMat and random order
## and write the full matrix with label
## Write step
dat    <- rbind(codMat, nonMat)
outDat <- paste(file_path_sans_ext(outFile), "_learningData.txt", sep="")
write.table(x=dat, file=outDat, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

## Random step
randomOrder <- sample(nrow(codMat)+nrow(nonMat))
dat         <- dat[randomOrder,]
dat.nameID  <- 1
dat.featID  <- (2:(ncol(dat)-1))
dat.labelID <- ncol(dat)


if(is.null(thres))
    {
        cat("\tNo threshold: running 10-fold cross-validation on learning set to set a threshold.\n")

        ## variables counting
        number_row   <- nrow(dat)
        chunk        <- list()
        output       <- list()
        models       <- list()
        nb_cross_val <- 10

        ## Progress bar
        cat("\t10-fold cross-validation progress:\n")
        progress <- txtProgressBar(1, nb_cross_val, style=3)
        setTxtProgressBar(progress, 0)

        ## Split in 'nb_cross_val' fold cross validation
        for (n in 1:nb_cross_val)
            {
                ## split the dat in 'nb_cross_val' chunks
                chunk[[n]] <-  seq(as.integer((n-1)*number_row/nb_cross_val)+1,as.integer(n*number_row/nb_cross_val))

                ## Train the random forest model with (nb_cross_val-1) chunks and predict the value for the test dat set
                models[[n]] <- randomForest(    x=dat[-chunk[[n]], dat.featID],    y=as.factor(dat[-chunk[[n]], dat.labelID]),
                                            xtest=dat[chunk[[n]], dat.featID], ytest=as.factor(dat[chunk[[n]], dat.labelID]),
                                            ntree=500)

                ## Output results in list output
                output[[n]] <- as.data.frame(cbind(dat[chunk[[n]], dat.featID], "Label"=dat[chunk[[n]], dat.labelID], "Prob"=models[[n]]$test$votes[,2], "Pred"=models[[n]]$test$predicted))

                setTxtProgressBar(progress, n)
            }
        cat("\n")

        ## Extract a list of coding probabilities and label
        allRes <- sapply(seq(1:nb_cross_val), function(i){output[[i]]$Prob})
        allLab <- sapply(seq(1:nb_cross_val), function(i){output[[i]]$Label})

        ## Generate the predict values using ROCR
        pred <- prediction(allRes, allLab)

        ## Get sensitivity and specificity
        S <- performance(pred,measure="sens")
        P <- performance(pred,measure="spec")

        ## Apply a function to get the mean on the 10 cutoffs that maximize the sens and spec (or minimize absolute difference : Thanks Oliver Sander)
        mean_cutoff <- mean(sapply(1:length(pred@predictions), function(i) { S@x.values[[i]][which.min(abs(S@y.values[[i]]-P@y.values[[i]]))] } ))
        mean_Sn     <- mean(sapply(1:length(pred@predictions), function(i) { P@y.values[[i]][which.min(abs(S@y.values[[i]]-P@y.values[[i]]))] } ))

######################################
##### BEGIN THE PLOT OF THE ROCR #####
######################################
        ## plot curve
        png(outROC, h=800, w=800)
        par(cex.axis=1.2, cex.lab=1.2)

        ymin=0.5
        ymax=1
        plot(S,col="blue",lty=3,ylab="Performance",xlab="Coding Probability Cutoff",ylim=c(ymin,ymax),cex.axis=1.2,cex.label=1.2, main="Two-Graph ROC curves")
        plot(S,lwd=2,avg="vertical",add=TRUE,col="blue")
        plot(P,col="red",lty=3, add=TRUE,)
        plot(P,lwd=2,avg="vertical",add=TRUE,col="red")

        ## Sn
        abline(h=mean_Sn,lty="dashed",lwd=0.5)
        text(x=0, y = mean_Sn, labels = round(mean_Sn, digits = 3), cex=1.5 )

        ## Cutoffs
        abline(v=mean_cutoff,lty="dashed",lwd=0.5)
        text(x=mean_cutoff, y = ymin, labels = round(mean_cutoff, digits = 3), cex=1.5)

        ## Legend
        legend("right",col=c("blue","red"),lwd=2,legend=c("Sensitivity","Specificity"))

        tt <- dev.off()
        cat("\tROCR curve plot in '", outROC, "'.\n", sep="")

####################################
##### END THE PLOT OF THE ROCR #####
####################################
        thres <- mean_cutoff

        cat("\t10-fold cross-validation step is finish. Best threshold found: '", thres, "'.\n", sep="")
    }
## END of: if(thres==NULL)


## Make the random forest model
cat("\tMaking random forest model on '", codFile, "' and '", nonFile ,"' and apply it on '", testFile, "'.\n", sep="")

## RF model
dat.rf <- randomForest(x=dat[,dat.featID], y=as.factor(dat[,dat.labelID]), ntree=500)

## Prediction on learning data
dat.rf.learn <- predict(dat.rf, dat[,dat.featID], cutoff=c(1-thres, thres))

## Prediction on test data
dat.rf.test <- predict(dat.rf, testMat[,dat.featID], cutoff=c(1-thres, thres))

## RF on learning for stats
## dat.rf.learn <- randomForest(x=dat[,dat.featID], y=as.factor(dat[,dat.labelID]),
##                        xtest=dat[,dat.featID],
##                        ntree=50)
## RF on test
## dat.rf <- randomForest(x=dat[,dat.featID], y=as.factor(dat[,dat.labelID]),
##                        xtest=testMat[,dat.featID],
##                        ntree=50)

## Get sensitivity, specificity, accuracy and precision on learning set using the model and the threshold
cat("\tPrinting stats in '", outStats, "' the sensitivity, specificity, precision and accuracy obtain on learning data.\n", sep="")
## res <- rep(0, length.out=nrow(dat))
## res[dat.rf.learn$test$votes[,2]>=thres] <- 1

## cont <- table(data=dat[,dat.labelID], prediction=res)
cont <- table(data=dat[,dat.labelID], prediction=dat.rf.learn)
tp   <- cont[2,2]
fp   <- cont[1,2]
tn   <- cont[1,1]
fn   <- cont[2,1]
sen  <- tp / (tp+fn)
spe  <- tn / (fp+tn)
pre  <- tp / (tp+fp)
acc  <- (tp+tn) / sum(cont)

cont <- cbind(c("true negative","false negative","false positive","true positive"), as.data.frame(cont)[,3])
write.table(x="Number of true/false positives/negatives:", file=outStats, quote=FALSE, sep="", col.names=FALSE, row.names=FALSE)
write.table(x=cont, file=outStats, append=TRUE, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
write.table(x="\nMetric values:", file=outStats, append=TRUE, sep="", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(x=cbind(c("Sensitivity","Specificity","Precision","Accuracy"), c(sen,spe,pre,acc)), file=outStats, append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

## Obtain the coding label prediction on test data
## res <- rep(0, length.out=nrow(testMat))
## res[dat.rf$test$votes[,2]>=thres] <- 1

## Write the output
cat("\tWrite the coding label for '", testFile, "' in '", outFile, "'.\n", sep="")
## write.table(x=cbind(testMat, label=res), file=outFile, quote=FALSE, sep="\t", row.names=FALSE)
write.table(x=cbind(testMat, label=dat.rf.test), file=outFile, quote=FALSE, sep="\t", row.names=FALSE)

## Write the plot for variable importance
cat("\tPlot the variable importance as measured by a random forest in '", outStats, "'.\n", sep="")
png(outVar, h=800, w=800)
varImpPlot(dat.rf)
tt <- dev.off()
