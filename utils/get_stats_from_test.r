#!/usr/bin/env Rscript

## Get stats from a test data

#########################
##### PREPROCESSING #####
#########################
## Checking the number of arguments
args <- commandArgs(TRUE)


## If number of argument != 2, then error
if(length(args)!=2)
{
    cat("Error: the number of argument is not 2... Quit\n")
    quit()
}

## Set generic values
testFile <- args[1]
rfFile   <- args[2]

## Include the tools package to get file with no extension
library(tools)

matTest <- read.table(file=testFile, header=TRUE, sep="\t")
matRf   <- read.table(file=rfFile,   header=TRUE, sep="\t")

labTest <- matTest[,ncol(matTest)]
labRf   <- matRf[,ncol(matRf)]

cont <- table(test = labTest, randomForest=labRf)
tp   <- cont[2,2]
fp   <- cont[1,2]
tn   <- cont[1,1]
fn   <- cont[2,1]

sen <- tp / (tp+fn)
spe <- tn / (fp+tn)
pre <- tp / (tp+fp)
acc <- (tp+tn) / sum(cont)

cat("Contingency:\n", cont, "\n", sep="")
cat("Stats:\n\tsensibility\t", sen, "\n\tspecificity\t", spe, "\n\tprecision\t", pre, "\n\taccuracy\t", acc, "\n", sep="")
