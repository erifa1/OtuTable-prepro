#!/usr/bin/Rscript

#Random sampling normalization 
##
##USAGE: 
#normalize_otu_table.R in.tabular.file min.abundance[numeric] out.file.name
#arg[1] must be tabular file, column are samples, row are features.
#arg[2] the minimum abundance wanted. Samples with total abundance lower than this value are discarded
#arg[3] Output file name.

#test args=c("laits_genus_table.csv", 6000, "laits_genus_table_norm.csv")


#setwd('')
args=commandArgs(trailingOnly=TRUE)
TABLE=as.character(args[1])
MinAb=as.numeric(args[2])
OUTname=as.character(args[3])
wd1=Sys.getenv('PWD')
if(length(args)<3){print("USAGE: normalize_otu_table.R in.tabular.file min.abundance[numeric] out.file.name")}else{


	otuTable <- as.matrix(read.table(paste(wd1,"/",TABLE,sep=""), sep = "\t", row.names = 1, h=TRUE))

	#Suppress low abandance samples
	SumSample=apply(otuTable, 2, sum)
	minAbondance=MinAb
	Tab2=otuTable[,SumSample>=MinAb]
	SampleSup=colnames(otuTable)[SumSample<MinAb]

	if (ncol(Tab2)==ncol(otuTable)){
		print("No sample discarded")} 
		else {
		print(paste("Sample with abundance < ",MinAb," :", sep=""))
		print(SampleSup)
	}

	cat("Processing....\n\n")

	A=Tab2
	#Random sampling normalisation
	OtuID=rownames(A)
	MINeff=min(apply(A[,2:ncol(A)],2, sum)); paste("min Abundance:",MINeff)
	stock=A
	for (j in 1:dim(stock)[2]){
		test=rep(as.character(OtuID),A[,j])
		sam1=sample(test,MINeff)
		for (i in 1:length(OtuID)){
			NOM=OtuID[i]
			COUNT=length(grep(paste("^",NOM,"$",sep=""),sam1))
			#print(c(i,COUNT))
			stock[i,j]=COUNT
		}
	}


	res=apply(stock[,1:dim(stock)[2]],2,sum)
	print("Total sum of each sample abundance:")
	print(res)

	final=data.frame(cbind(OtuID,stock))

	write.table(final, OUTname, quote=FALSE, sep="\t", row.names=FALSE)

}
