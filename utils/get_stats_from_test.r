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

cont <- table(labTest, labRf, dnn=list(basename(testFile), basename(rfFile)))
if(all(dim(cont) == c(2,2)))
    {
        tp    <- cont[2,2]
        fp    <- cont[1,2]
        tn    <- cont[1,1]
        fn    <- cont[2,1]
        contS <- sum(cont)
    } else if(all(dim(cont) == c(2,3))) {
        tp    <- cont[2,3]
        fp    <- cont[1,3]
        tn    <- cont[1,2]
        fn    <- cont[2,2]
        contS <- sum(cont[,(2:3)])
    } else {
        cat("Error in the contingency table... Quit\n")
        quit()
    }


sen  <- round( (tp / (tp+fn)), digits=3)
spe  <- round( (tn / (fp+tn)), digits=3)
pre  <- round( (tp / (tp+fp)), digits=3)
acc  <- round( ((tp+tn) / contS), digits=3)
fsc  <- round( (2*tp) / (2*tp+fp+fn), digits=3)
mcc  <- round( ((tp*tn)-(fp*fn)) / ( sqrt(tp+fp) * sqrt(tp+fn) * sqrt(tn+fp) * sqrt(tn+fn)  ), digits=3)

cat("Contingency (-1: TUCp; 0: lncRNA; 1: mRNA):\n")
print(cont)
cat("Stats (on mRNA):\nsensibility\t", sen, "\nspecificity\t", spe, "\nprecision\t", pre, "\naccuracy\t", acc, "\nf1-score\t", fsc, "\n", phi-mcc\t", mcc, "\n",sep="")
