#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(pegas))

filenames<-Sys.glob("*.fa")
d<-matrix(NA,ncol=7,nrow=(length(filenames)))
colnames(d)<-c("haplo","hom","het","seg.sites","theta","nuc.div","GC")
rownames(d)<-paste("loc",1:(length(filenames)),sep="")

for (f in (1:length(filenames))) {
	currFile<-filenames[f]
	dna<-read.dna(currFile,format="fasta")
	hap<-haplotype(dna)
	
	a<-array(NA,nrow(dna))
	m<-matrix(NA,nrow=(nrow(dna)/2),ncol=2)

	for (i in 1:(length(attr(hap,"index")))) { 
  		for (j in 1:(length(attr(hap,"index")[[i]]))) { 
    			a[attr(hap,"index")[[i]][j]]<-i ## insert haplotype number into array a at appropriate index
  		} 
	}

	m[,1]<-a[seq(1,(nrow(dna))-1,2)] ## odd indices, allele A
	m[,2]<-a[seq(2,(nrow(dna)),2)] ## even indices, allele B
	m<-t(apply(m,1,sort)) ## sort so that lower haplotype number always first

	## summarise m
	hom<-0
	het<-0
	for (r in (1:nrow(m))) {
  		if (m[r,1]==m[r,2]) {
		    hom<-hom+1
  		} else {
    			het<-het+1
  		}
	}
#	cat("File:",currFile,"\n")
#	cat("Number of haplotypes:",nrow(hap),"\n")
#	cat("Homozygotes:",hom,"\n")
#	cat("Heterozygotes:",het,"\n\n")

#	cat("Number segregating sites:",length(seg.sites(dna)),"\n")
#	cat("Average pairwise differences:",theta.s(length(seg.sites(dna)),nrow(dna)),"\n")
#	cat("Nucleotide diversity:",nuc.div(dna),"\n")
#	cat("%GC:",GC.content(dna),"\n\n")


	#print(i)
	d[f,1]<-(nrow(hap))
	d[f,2]<-hom
	d[f,3]<-het
	d[f,4]<-length(seg.sites(dna))
	d[f,5]<-theta.s(length(seg.sites(dna)),nrow(dna))
	d[f,6]<-nuc.div(dna)
	d[f,7]<-GC.content(dna)
}
mean<-colMeans(d,na.rm=TRUE)
d<-rbind(d,mean)
write.csv(d,file="bdevolver.stats")
print(d)
