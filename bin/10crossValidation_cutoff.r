# Rscript to compute optimal cutoff
# adpated from Liguo Wang : http://dldcc-web.brc.bcm.edu/lilab/liguow/CGI/cpat/_build/html/index.html

# Install ROCR package if missing
# tx to : http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
# see also : http://statistics.berkeley.edu/computing/R-packages
# install in pwd by default
# Note:  Better is that user types install.packages('ROCR') in a R session
pwd <- Sys.getenv("PWD");
.libPaths(c(pwd,.libPaths()))
list.of.packages <- c("ROCR")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos="http://cran.r-project.org/",  dependencies = TRUE)

# Load lib
library('ROCR')

# CPAT output file (lncRNA and mRNA) with labels in argument
args <- commandArgs(trailingOnly = TRUE)
cpatinfile <-args[1]


data=read.table(file=cpatinfile,header=T,sep="\t")
attach(data)

# number of cross validation
nb_cross_val = 10

# variables counting 
number_row = nrow(data)
chunk <-list()
output <-list()

# Split in 'nb_cross_val' fold cross validation
for (n in 1:nb_cross_val){


	# split the data in 'nb_cross_val' chunks
	chunk[[n]] <-  seq(as.integer((n-1)*number_row/nb_cross_val)+1,as.integer(n*number_row/nb_cross_val))

	# extract data *without* the one chunk of interest
	vlabel = label[-chunk[[n]] ]
	vmrna = mRNA_size[-chunk[[n]] ]
	vorf = ORF_size[-chunk[[n]] ]
	vfickett = Fickett_score[-chunk[[n]] ]
	vhexamer = Hexamer_score[-chunk[[n]] ]
	
	# Train the logit model with 9 (or (nb_cross_val-1) chunks)
	mylogit <- glm(vlabel ~ vmrna + vorf + vfickett + vhexamer, family=binomial(link="logit"), na.action=na.pass)
	test <- data.frame(vmrna = mRNA_size[chunk[[n]] ], vorf = ORF_size[chunk[[n]] ],vfickett = Fickett_score[chunk[[n]] ], vhexamer = Hexamer_score[chunk[[n]] ], vlabel=label[chunk[[n]] ])
	
	# Predict value for the chunk
	test$prob <- predict(mylogit,newdata=test,type="response")
	# Output results in list output
	output[[n]] = as.data.frame(cbind("mRNA"=test$vmrna, "ORF"=test$vorf, "Fickett"=test$vfickett, "Hexamer" = test$vhexamer, "Label"=test$vlabel,"Prob"=test$prob))

}

# Extract a list of coding probabilities
Response = sapply(seq(1:nb_cross_val), function(i){output[[i]]$Prob})
# Extract a list of coding Label
Labls = sapply(seq(1:nb_cross_val), function(i){output[[i]]$Label})

# ROCR
ROCR_data = list(predictions=Response,Labels=Labls)
pred <- prediction(ROCR_data$predictions, ROCR_data$Labels)
#perf <- performance(pred,"auc")



png(paste(cpatinfile,".png",sep=""), h=800, w=800)
par(mfrow=c(2,2),mar=c(5,4,2,2),cex.axis=1.2, cex.lab=1.2)
#ROC curve
#pdf("Human_10fold.ROC.pdf")
perf <- performance(pred,"tpr","fpr")
plot(perf,col="blue",lty=3,xlab="1-Specificity",ylab="Sensitivity",ylim=c(0.7,1),xlim=c(0,0.3),main="",cex.axis=1.5,cex.label=1.5)	
plot(perf,lwd=2,avg="vertical",add=TRUE,col="red",xlab="1-specificity",ylab="sensitivity",main="",cex.axis=1.2,cex.label=1.2) 
#dev.off()

#precision
#pdf("Human_10fold.precision_vs_recall.pdf")
d=performance(pred,measure="prec", x.measure="rec")
plot(d,col="blue",lty=3,xlab="Recall (TPR)",ylab="Precision (PPV)",xlim=c(0.7,1),ylim=c(0.7,1),cex.axis=1.2,cex.label=1.2)
plot(d,lwd=2,avg="vertical",col="red",xlab="Recall (TPR)",ylab="Precision (PPV)",add=T,cex.axis=1.2,cex.label=1.2)
#dev.off()


#Accuracy
#pdf("Human_10fold.Accuracy.pdf")
perf <- performance(pred,"acc")
plot(perf,col="blue",lty=3,xlab="Coding probability cutoff",ylab="Accuracy",ylim=c(0.7,1),cex.axis=1.2,cex.label=1.2) 
plot(perf,lwd=2,avg="vertical",add=TRUE,col="red",cex.axis=1.2,cex.label=1.2) 
#dev.off()


#sensitivity vs specificity
S <- performance(pred,measure="sens")
P <- performance(pred,measure="spec")

# Apply a function to get the 10 cutoffs that maximize the sens and spec (or minimize absolute difference : Thanks Oliver Sander 
mean_cutoff=mean(sapply(1:length(pred@predictions), function(i) { S@x.values[[i]][which.min(abs(S@y.values[[i]]-P@y.values[[i]]))] } ))
mean_Sn=mean(sapply(1:length(pred@predictions), function(i) { P@y.values[[i]][which.min(abs(S@y.values[[i]]-P@y.values[[i]]))] } ))
# mean_sn == mean_Sp

#pdf("Human_10fold_sens_vs_spec.pdf")
ymin=0.5
ymax=1
plot(S,col="blue",lty=3,ylab="Performance",xlab="Coding Probability Cutoff",ylim=c(ymin,ymax),cex.axis=1.2,cex.label=1.2) 
plot(S,lwd=2,avg="vertical",add=TRUE,col="blue") 
plot(P,col="red",lty=3, add=TRUE,) 
plot(P,lwd=2,avg="vertical",add=TRUE,col="red") 

# Sn
abline(h=mean_Sn,lty="dashed",lwd=0.5)
text(x=0, y = mean_Sn, labels = round(mean_Sn, digits = 3), cex=1.5 ) 

# Cutoffs
abline(v=mean_cutoff,lty="dashed",lwd=0.5)
text(x=mean_cutoff, y = ymin, labels = round(mean_cutoff, digits = 3), cex=1.5)

# Legend
legend("right",col=c("blue","red"),lwd=2,legend=c("Sensitivity","Specificity"))

# Print Avg Cutoff and Sn/Sp
print (paste ("Cutoff_SnSp", mean_cutoff,mean_Sn, sep=" "))


dev.off()
