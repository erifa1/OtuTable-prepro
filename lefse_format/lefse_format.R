#!/usr/bin/Rscript --slave
#Rscript lefse_format.R lefse_format_input_table_example.csv lefse_format_mappingfile_example.csv 3 > log.Rout

#WARNING
#USE tabular file with exactly the same shape and column names. 7 tax rank in the seven first column, otuid in 8th column, counts from 9th column

#rm(list=ls())
args=commandArgs(trailingOnly=TRUE)

otu_table=as.character(args[1])
mapping_file=as.character(args[2])
index_fact_mapping_file=as.numeric(args[3])

if(length(args)<3){print("USAGE: Rscript lefse_format.R ABUNDANCE_TABLE.csv[tab separated] MAPPING_FILE.csv[tab separated] INDEX_FACTOR")}else{

	wd1=Sys.getenv('PWD')
	setwd(wd1)
	#setwd("/home/pri/Desktop/farm2cheese/lefse_input")
	A=read.table(otu_table, h=T)
	B=read.table(mapping_file, h=T)



	Btr=data.frame(t(B))

	ncol=dim(A)[2]

	#Domain counts (fungi, bacteria, archeae)
	Total=NULL
	for(i in 9:ncol){
		res=by(A[,i],as.factor(A[,1]),sum)
		Total=cbind(Total,res)
	}
	colnames(Total)=names(A)[9:ncol]

	#Sum counts for each tax group
	for(j in 2:7){
		print(j)
		lvl=as.factor(apply(A[,1:j], 1,paste,collapse="|")) #create taxonomy for each rank
		  
		stock=matrix(ncol=length(9:ncol),nrow=length(levels(lvl)))
		rownames(stock)=levels(lvl)
		for(i in 9:ncol){
			print(i)
			stock[,i-8]=as.numeric( by(A[,i],lvl,sum) ) #sum counts of each levels in the factor lvl
		}

		Total=rbind(Total,stock)	#concat with domain counts

	}

	#Finalizing lefse table
	Samples=colnames(Total)
	Fact=NULL
	for(i in 1:length(Samples)){
		Fact[i]=as.character(B[B$id==Samples[i],index_fact_mapping_file])
	}

	Final=rbind(Samples,Fact,Total)


	write.table(Final, "table_lefse.tsv",sep="\t",quote = FALSE, col.names=FALSE)

} #if
