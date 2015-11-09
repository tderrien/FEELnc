#!/usr/bin/env Rscript

## Running the Random Forest to classify transcripts between coding and non coding
## Need the randomForest and ROCR library

######################################
##### Version for two thresholds #####
######################################

#########################
##### PREPROCESSING #####
#########################
## Include the tools package to get file with no extension
library(tools)

## Checking the number of arguments
args <- commandArgs(TRUE)

## Set generic values
codFile    <- args[1]
nonFile    <- args[2]
testFile   <- args[3]
outFile    <- args[4]
numberT    <- as.numeric(args[5])
seed       <- as.numeric(args[6])
thresSpeM  <- as.numeric(args[7])
thresSpeL  <- as.numeric(args[8])
outSummary <- paste(file_path_sans_ext(outFile), "_summary.txt", sep="")
outVar     <- paste(file_path_sans_ext(outFile), "_varImpPlot.png", sep="")
outROC     <- paste(file_path_sans_ext(outFile), "_TGROC.png", sep="")
outStats   <- paste(file_path_sans_ext(outFile), "_statsLearn_CrossValidation.txt", sep="")
list.of.packages <- c("ROCR","randomForest")

if(length(args) != 8)
    {
        cat("Error: the number of argument pass to codpot_randomforest_2thres.r is wrong:\ncodpot_randomforest_2thres.r inCoding inNonCoding testFile outFile nTree seed specificityThresholdMrna specificityThresholdLncrna\nQuit\n")
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


## Fix seed
set.seed(seed)

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

## Random ordering of input matrix (fixe by seed)
randomOrder  <- sample(nrow(codMat)+nrow(nonMat))
dat          <- dat[randomOrder,]
dat.nameID   <- 1
dat.featID   <- (2:(ncol(dat)-1))
dat.labelID  <- ncol(dat)
test.nameID  <- 1
test.featID  <- (2:ncol(testMat))

## variables counting
number_row   <- nrow(dat)
chunk        <- list()
output       <- list()
models       <- list()
nb_cross_val <- 10
models.votes <- list()

## Progress bar
cat("\tRunning ", nb_cross_val, "-fold cross-validation on learning:\n", sep="")
progress <- txtProgressBar(0, nb_cross_val, style=3)
setTxtProgressBar(progress, 0)


## Split in 'nb_cross_val' fold cross validation
for (n in 1:nb_cross_val)
    {
        ## split the dat in 'nb_cross_val' chunks
        chunk[[n]] <-  seq(as.integer((n-1)*number_row/nb_cross_val)+1,as.integer(n*number_row/nb_cross_val))

        ## Get the minimum number of non coding and coding values
        minVal <- min(table(as.factor(dat[-chunk[[n]], dat.labelID])))

        ## Train the random forest model with (nb_cross_val-1) chunks and predict the value for the test dat set
        ## with an equal number of lncRNAs and mRNAs in each tree
        models[[n]] <- randomForest(    x=dat[-chunk[[n]], dat.featID], y=as.factor(dat[-chunk[[n]], dat.labelID]),
                                    ntree=numberT, sampsize=c("0"=(minVal), "1"=(minVal)))
        models.votes[[n]] <- predict(models[[n]], dat[chunk[[n]], dat.featID], type="vote")

        ## Output results in list output
        output[[n]] <- as.data.frame(cbind(dat[chunk[[n]], dat.featID], "Label"=dat[chunk[[n]], dat.labelID], "Prob"=models.votes[[n]][,2]))

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


## If the cutting between the spe and sens is higher than the giving specificity threshold, then use classical threshold
if(mean_Sn >= thresSpeM & mean_Sn >= thresSpeL)
    {
        cat("WARNING: the value where sensitivity equal specicifity: '", mean_Sn,"' is greater than the specificity threshold: mRNA: '", thresSpeM, "'; lncRNA: '", thresSpeL, "'. Use the best value.\n", sep="")

        cutoffThresSpeM <- mean_cutoff
        cutoffThresSpeL <- mean_cutoff
        thresSpeM       <- mean_Sn
        thresSpeL       <- mean_Sn
    } else {
        ## Get the threshold for each model with a mRNA specificity >= thresSpeM
        cutoffThresSpeM <- mean( sapply(1:length(pred@predictions), function(i) { min(P@x.values[[i]][which(P@y.values[[i]]>=thresSpeM)]) } ) )
        ## Same for thresSpeL
        cutoffThresSpeL <- mean( sapply(1:length(pred@predictions), function(i) { max(S@x.values[[i]][which(S@y.values[[i]]>=thresSpeL)]) } ) )

        cat("\t10-fold cross-validation step is finish. Best thresholds found: mRNA: '", cutoffThresSpeM, "'; lncRNA: '", cutoffThresSpeL,"'.\n", sep="")
    }

if (cutoffThresSpeM <= cutoffThresSpeL) {
    cat("WARNING: the threshold obtain for mRNA '", cutoffThresSpeM, "' is lesser than the the one for lncRNA '", cutoffThresSpeL, "'. Use the threshold '", mean_cutoff, "' where sensitivity equal to specicifity for mRNA: '", mean_Sn,"'.\n", sep="")

    cutoffThresSpeM <- mean_cutoff
    cutoffThresSpeL <- mean_cutoff
    thresSpeM       <- mean_Sn
    thresSpeL       <- mean_Sn
    }


## Print in outStats the sensitivity, specificity, accuracy and precision
cat("\tPrinting stats found with the 10-fold cross-validation in '", outStats, "'.\n", sep="")

res <- matrix(0, ncol=10, nrow=(nb_cross_val+1), dimnames=list(c((1:nb_cross_val),"mean"),c("sen","spe","pre","acc","tp","tn","fp","fn","TUCp_mRNA","TUCp_lncRNA")))
for(i in 1:nb_cross_val)
    {
        mod.lab                                         <- rep(-1, length.out=nrow(dat[chunk[[i]],]))
        mod.lab[models.votes[[i]][,2]>=cutoffThresSpeM] <- 1
        mod.lab[models.votes[[i]][,2]<=cutoffThresSpeL] <- 0

        cont2 <- table(dat[chunk[[i]],dat.labelID], mod.lab)
        if(all(dim(cont2) == c(2,2)))
            {
                cont  <- cont2
                tm    <- 0
                tl    <- 0
            } else {
                cont  <- cont2[,-1]
                tm    <- cont2[2,1]
                tl    <- cont2[1,1]
            }

        tp <- cont[2,2]
        fp <- cont[1,2]
        tn <- cont[1,1]
        fn <- cont[2,1]

        sen <- tp / (tp+fn)
        spe <- tn / (fp+tn)
        pre <- tp / (tp+fp)
        acc <- (tp+tn) / sum(cont)

        res[i,] <- c(sen,spe,pre,acc,tp,tn,fp,fn,tm,tl)
    }
res[nrow(res),] <- colMeans(res[-nrow(res),])
res             <- cbind(mod=rownames(res), round(res, digits=2))

write.table(x=res, file=outStats, quote=FALSE, sep="\t", row.names=FALSE)


######################################
##### BEGIN THE PLOT OF THE ROCR #####
######################################
## plot curve
cat("\tTwo-graphs ROCR curves in '", outROC, "'.\n", sep="")
png(outROC, h=800, w=800)
par(cex.axis=1.2, cex.lab=1.2)

ymin=0.5
ymax=1
plot(S,col="blue",lty=3,ylab="Performance",xlab="Coding Probability Cutoff",ylim=c(ymin,ymax),cex.axis=1.2,cex.label=1.2, main="Two-Graph ROC curves")
plot(S,lwd=2,avg="vertical",add=TRUE,col="blue")
plot(P,col="red",lty=3, add=TRUE)
plot(P,lwd=2,avg="vertical",add=TRUE,col="red")

## Specificity line (red: mRNA; blue: lncRNA)
abline(h=thresSpeM,lty="dashed",lwd=1.5, col="red")
abline(h=thresSpeL,lty="dashed",lwd=1.5, col="blue")
text(x=0, y = thresSpeM, col="red",  labels = round(thresSpeM, digits = 3), cex=1.5, pos=3 )
text(x=0, y = thresSpeL, col="blue", labels = round(thresSpeL, digits = 3), cex=1.5, pos=1 )

## Cutoffs (red: mRNA; blue: lncRNA)
abline(v=cutoffThresSpeM,lty="dashed",lwd=1.5, col="red")
abline(v=cutoffThresSpeL,lty="dashed",lwd=1.5, col="blue")
text(x=cutoffThresSpeM, y = ymin, labels = round(cutoffThresSpeM, digits = 3), cex=1.5, pos=4, col="red")
text(x=cutoffThresSpeL, y = ymin, labels = round(cutoffThresSpeL, digits = 3), cex=1.5, pos=2, col="blue")

## Legend
legend("right",col=c("blue","red"),lwd=2,legend=c("mRNA sensitivity","mRNA specificity"))

tt <- dev.off()

####################################
##### END THE PLOT OF THE ROCR #####
####################################


## Make the random forest model
cat("\tMaking random forest model on '", basename(codFile), "' and '", basename(nonFile), "' and apply it to '", basename(testFile), "'.\n", sep="")

## Get the minimum number of non coding and coding values
minVal <- min(table(as.factor(dat[, dat.labelID])))
## RF model
dat.rf <- randomForest(x=dat[,dat.featID], y=as.factor(dat[,dat.labelID]), ntree=numberT, sampsize=c("0"=(minVal), "1"=(minVal)))

## Prediction on test data
dat.rf.test.votes                                   <- predict(dat.rf, testMat[,test.featID], type="vote")
dat.rf.test                                         <- rep(-1, length.out=nrow(testMat))
dat.rf.test[dat.rf.test.votes[,2]>=cutoffThresSpeM] <- 1
dat.rf.test[dat.rf.test.votes[,2]<cutoffThresSpeL]  <- 0


## Write the output
cat("\tWrite the coding label for '", basename(testFile), "' in '", outFile, "'.\n", sep="")
write.table(x=cbind(testMat, coding_potential=dat.rf.test.votes[,2], label=dat.rf.test), file=outFile, quote=FALSE, sep="\t", row.names=FALSE)

## Write the summary file
## If there is only one cutoff
if(cutoffThresSpeM==cutoffThresSpeL)
{
	nbtuc	=	0
	nblnc	=	table(dat.rf.test)[1]
	nbmrna	=	table(dat.rf.test)[2]
} else {
	nbtuc	=	table(dat.rf.test)[1]
	nblnc	=	table(dat.rf.test)[2]
	nbmrna	=	table(dat.rf.test)[3]
}
cat("# Summary file:\n-With_mRNA_cutoff:\t",cutoffThresSpeM,"\n-With_lncRNA_cutoff:\t",cutoffThresSpeL," \n-Nb_TUCPs:\t",nbtuc,"\n-Nb_lncRNAs:\t",nblnc,"\n-Nb_mRNAs:\t",nbmrna,"\n", file=outSummary, sep = "")



## Write the plot for variable importance
cat("\tPlot the variable importance as measured by random forest in '", outStats, "'.\n", sep="")
png(outVar, h=800, w=800)
varImpPlot(dat.rf)
tt <- dev.off()
