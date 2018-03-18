#!/usr/bin/Rscript

#merge taxons of otu table without normalization. 
#TABLE="otu_table_example.tsv"; INDEX2=29; RANK=5; MINab=750; OUT1="test.tsv"

#INPUT
#Tabular formatted OTU table  as "otu_table_example.tsv"
#
#args
# 1: input OTU table tabular separated with "taxonomy" column in the last column, taxonomy separated with semicolon
# 2: index of the last sample column
# 3: taxonomic rank to output, i.e. 5 for family 6 for genus
# 4: min total abundance to discard low abundance sample
# 5: output name 


args=commandArgs(trailingOnly=TRUE)
TABLE=as.character(args[1])
INDEX2=as.numeric(args[2])
RANK=as.numeric(args[3])
MINab=as.numeric(args[4])
OUT1=as.character(args[5])

wd1=Sys.getenv('PWD')
if(length(args)<5){print("USAGE: taxon_table.R input_table.tsv[chr] index_last_sample[num] rank[num] min_abundance[num] name_outpu_table[chr]")
cat("\n\nINPUT\n#Tabular formatted OTU table  as otu_table_example.tsv
args
 1: input OTU table tabular separated with taxonomy column in the last column, taxonomy separated with semicolon
 2: index of the last sample column
 3: taxonomic rank to output, i.e. 5 for family 6 for genus
 4: min total abundance to discard low abundance sample
 5: output name \n")}else{

	#setwd("/home/pri/Desktop/farm2cheese/DATAS/00_contigs/00_fishflow/phyloseq")

	A=read.table(TABLE, sep="\t", h=TRUE)
	Asample=A[,2:INDEX2]	#84
	rownames(Asample)=A[,1]

	#Taxonomy extraction
	nsplt=strsplit(as.character(A$taxonomy), "[;]")
	namesTAB=as.character(do.call(rbind.data.frame, nsplt)[,RANK])
	

	MINeff=min(apply(Asample,1,sum))
	#Transpose
	B1=as.data.frame(Asample)

	#Sum multiple occurence of the same taxon
	stock1=NULL
	for (i in 1:dim(B1)[2]){
		#print(names(B1)[i])
		res=tapply(B1[,i],namesTAB,sum)
		stock1=cbind(stock1, res)
	}
	
	B=as.data.frame(stock1)
	colnames(B)=colnames(B1)

	#Discard sample with low abundance
	sum1=apply(B, 2, sum)
	C=B[,sum1>=MINab]
	if(ncol(B)==ncol(C)){
		print("No sample discarded")} else {
		print("Sample discarded:");print(names(B[,sum1<MINab]))
	}

	D=cbind.data.frame(taxon=row.names(C),C)

	#Output table
	write.table(D, OUT1, sep="\t", quote=FALSE, row.names=FALSE)
}
