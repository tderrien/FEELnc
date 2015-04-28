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

matTest <- read.table(file=testFile, header=TRUE, sep="\t")
matRf   <- read.table(file=rfFile,   header=TRUE, sep="\t")

rownames(matTest) <- matTest[,1]
rownames(matRf)   <- matRf[,1]

test.namesSort <- sort(as.character(matTest[,1]))
rf.namesSort   <- sort(as.character(matRf[,1]))

matTest <- matTest[test.namesSort,]
matRf   <- matRf[rf.namesSort,]

checkIn <- rownames(matTest) %in% rownames(matRf)

labTest <- matTest[checkIn,ncol(matTest)]
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

cat("Contingency:\n")
print(cont)
cat("Stats:\nsensibility\t", sen, "\nspecificity\t", spe, "\nprecision\t", pre, "\naccuracy\t", acc, "\n", sep="")
