## load option parsing lib
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(pegas))
option_list = list(
  make_option(c("-f","--fasta"), type="character", default="input.fasta", help="Input haplotype sequences in fasta format [default=%default]"),
  make_option(c("-s","--stats"), type = "character", default = "input.stats", help = "Input stats file [default=%default]")
#  make_option(c("-G","--graphic"), type = "character", default = "png", help = "Type of graphical output desired (can be png, jpeg, svg, pdf) [default=%default]")
);
opt_parser=OptionParser(option_list=option_list);
opt=parse_args(opt_parser);

## get dna and stats
cat("\nFasta:",opt$fasta,"\n")
cat("Stats:",opt$stats,"\n")
dna<-read.dna(opt$fasta,format="fasta")
stats<-read.table(opt$stats,head=TRUE)

## plot from stats
filename<-paste("theta.",opt$fasta,".svg",sep="",collapse="")
svg(filename)
plot(stats$gen,stats$theta.s,type="l",xlab="generations",ylab="average pairwise differences (theta)",main="",col="darkgrey")
garbage<-dev.off() ## prevents "null device" message going to STDOUT
filename<-paste("nuc_div.",opt$fasta,".svg",sep="",collapse="")
svg(filename)
plot(stats$gen,stats$nuc.div,type="l",xlab="generations",ylab="nucleotide diversity",main="",col="darkgrey")
garbage<-dev.off()
filename<-paste("haplotypes.",opt$fasta,".svg",sep="",collapse="")
svg(filename)
plot(stats$gen,stats$haplo,type="l",xlab="generations",ylab="haplotypes",main="",col="darkgrey")
garbage<-dev.off()

## make haplonet
hap<-haplotype(dna)
net<-haploNet(hap)

radii<-array(NA,nrow(hap)) ## get number of haps as NA array
for (i in 1:(nrow(hap))){ radii[i]<-length(attr(hap,"index")[[i]]) } ## populate with nseqs in each hap category

filename<-paste("net.",opt$fasta,".svg",sep="",collapse="")
svg(filename)
plot(net,size=radii,threshold=0) ## plot network without alternative mutational pathways drawn
garbage<-dev.off()

################################################################################
## analyse hap #
################

## generate NA array equal to number of alleles in sample
a<-array(NA,nrow(dna))

## generate NA matrix where rows == individual cols == genotype
m<-matrix(NA,nrow=(nrow(dna)/2),ncol=2)

for (i in 1:(length(attr(hap,"index")))) {
  for (j in 1:(length(attr(hap,"index")[[i]]))) {
    a[attr(hap,"index")[[i]][j]]<-i ## insert haplotype number into array a at appropriate index
  }
}

## then put genotype info into matrix m
m[,1]<-a[seq(1,(nrow(dna))-1,2)] ## odd indices, allele A
m[,2]<-a[seq(2,(nrow(dna)),2)] ## even indices, allele B
m<-t(apply(m,1,sort)) ## sort so that lower haplotype number always first

## summarise m
hom<-0
het<-0
for (r in (1:nrow(m))) {
  if (m[r,1]==m[r,2]) {
    #    cat("Ind",r,"is homozygote:",m[r,],"\n")
    hom<-hom+1
  } else {
    het<-het+1
  }
}
cat("Number of haplotypes:",nrow(hap),"\n")
cat("Homozygotes:",hom,"\n")
cat("Heterozygotes:",het,"\n")
cat("Number segregating sites:",length(seg.sites(dna)),"\n")
cat("Average pairwise differences:",theta.s(length(seg.sites(dna)),nrow(dna)),"\n")
cat("Nucleotide diversity:",nuc.div(dna),"\n")
cat("%GC:",GC.content(dna),"\n~~~\n")

##########################################################################################
## triplets is best done as sliding window (use bdevolver.triplets.R)

## make table that sums all genotypes
#tab<-matrix(0,nrow=nrow(hap),ncol=nrow(hap))
#for (i in (1:nrow(m))) {
#  tab[m[i,1],m[i,2]]<-tab[m[i,1],m[i,2]]+1
#}
## NB get table like
#     [,1] [,2] [,3] [,4]
#[1,]    3    2    1    1
#[2,]    2    1    1    0
#[3,]    2    1    0    1
#[4,]    3    0    0    2
# which shows the number of inds with each genotype, ie 3 inds with genotype 1|1, 2 with 1|2 etc.

## generate all possible combinations of 3 haplotypes
## cycles has rows = 1:3, haplotype number, cols = all possible combo's
#cycles<-combn(1:nrow(hap),3)
#colnames(cycles)<-rep("FALSE",ncol(cycles))
#for (j in (1:ncol(cycles))) {
#  count<-0
#  if (tab[cycles[1,j],cycles[2,j]]>0) count<-count+1
#  if (tab[cycles[1,j],cycles[3,j]]>0) count<-count+1
#  if (tab[cycles[2,j],cycles[3,j]]>0) count<-count+1
#  if (count==3) {colnames(cycles)[j]<-"TRUE"} else {colnames(cycles)[j]<-"FALSE"}
#}
## LOGIC
## For there to be a cycle or 'triplet', need 3 heterozygous inds
## with genotypes 1|2, 1|3, and 2|3
## in the table above, that looks like:
##      [,1] [,2] [,3]
## [1,]    0    1    1
## [2,]    1    0    1
## [3,]    1    1    0

##number of cycles
#triplets<-sum(as.logical(colnames(cycles)))
#cat("Number of triplet genotypes (a|b, a|c, b|c):",triplets,"\n")
#if (triplets>0) {
#  for (col in (1:ncol(cycles))) {
#    if (colnames(cycles)[col]==TRUE) {
#      alleles<-t(combn(cycles[,col],2))
#      cat("\tAlleles: ")
#      cat(cycles[,col],sep=", ")
#      cat("; found in:\n")
#      for (i in 1:nrow(m)) {
#        for (j in 1:nrow(hap)) {
#          if(m[i,1]==alleles[j,1] & m[i,2]==alleles[j,2]){
#            cat("\t\tInd_",i," (genotype ",sep="")
#            cat(m[i,],sep="|")
#            cat(")\n",sep="")
#          }
#        }
#      }
#    }
#  }
#}

##proportion of possible cycles observed
#triplets.prop<-sum(as.logical(colnames(cycles)))/ncol(cycles)
#cat("As a proportion of possible triplets:",triplets.prop,"\n")
