#!/usr/bin/env Rscript

###########################################################
## sliding window
#################

## load option parsing lib
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(pegas))
option_list = list(
	make_option(c("-f","--fasta"), type="character", help="Input haplotype sequences in fasta format [REQUIRED]"),
	make_option(c("-w","--window"), type="integer", default=10, help="Window size for sliding window [default=%default]"),
	make_option(c("-s","--step"), type="integer", default=1, help="Step size for sliding window [default=%default]"),
	make_option(c("-t","--strip"), action='store_true', default=FALSE, help="Strip monomorphic alignment columns [default=%default]"),
	make_option(c("-d","--decreasing"), action='store_true', default=FALSE, help="Use decreasing window approach, ie iterate thru windows of {n..1} where n=alignment length, breaking at the 1st instance of a haplotype trio is found (ie, the longest window) [default=%default]"),
	make_option(c("-q","--quitonfirst"), action='store_true', default=FALSE, help="Strip monomorphic alignment columns [default=%default]"),
	make_option(c("-v","--verbose"), type="integer", default=1, help="Level of verbosity (>1 prints all individuals in haplotype trio, plus their genotype) [default=%default]")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## get dna from fasta
dna<-read.dna(opt$fasta,format="fasta")
#dna<-read.dna("sim3.loc1.test.fna.triplet.fasta",format="fasta")
#opt$strip<-TRUE

if (opt$verbose>0) {
  cat("[START]\n\n")
  cat("Input filename:",opt$fasta,"\n")
  cat("Window size:",opt$window,"\n")
  cat("Step size:",opt$step,"\n")
  cat("Strip invariant sites:",opt$strip,"\n")
  cat("Decreasing window:",opt$decreasing,"\n")
}

## strip invariant sites in alignment
if (opt$strip==TRUE) {
  before<-length(dna[1,])
  if (opt$verbose>0) cat("\nStripping invariant sites from alignment... (length from ",before," to ",sep="")
  df<-as.data.frame(as.character(dna)) ## convert to df
  df<-df[sapply(df, function(x) length(unique(x))>1)] ## keep only aln cols with SNPs in them
  dna<-as.DNAbin(as.matrix(df)) ## convert back to DNAbin
  after<-length(dna[1,])
  if (opt$verbose>0) cat(after,")\n",sep="")
}

if (opt$decreasing==FALSE) {
  starts<-seq(1,length(dna[1,]),by=opt$step) ## gives start coordinates for all windows
  ends<-starts+(opt$window-1) ## end coords for all windows
  starts<-starts[ends<=length(dna[1,])] ## trims window coordinates to fit within seq coords
  ends<-ends[ends<=length(dna[1,])]
  if (opt$verbose>0) cat("Number of windows:",length(starts),"\n")
  
  #chunks.triplets<-list() ## initialise list to store triplet window alignments as DNAbin objs
  windows.triplets<-rep(0,length(starts)) ## initialise a vector of length == number windows
  
  for (window.index in (1:length(starts))) {
    ## get chunk and make hap from it
    chunk.dna<-dna[,starts[window.index]:ends[window.index]]
    if (opt$verbose>1) cat(starts[window.index],":",ends[window.index],"\n",sep="")
    chunk.hap<-haplotype(chunk.dna)
    
    if (nrow(chunk.hap)>=3) { ## no point testing windows with <3 haplotypes
      a<-array(NA,nrow(chunk.dna)) ## generate NA array equal to number of alleles in sample
      m<-matrix(NA,nrow=(nrow(chunk.dna)/2),ncol=2) ## generate NA matrix where rows == individual cols == genotype
      
      for (i in 1:(length(attr(chunk.hap,"index")))) { 
        for (j in 1:(length(attr(chunk.hap,"index")[[i]]))) { 
          a[attr(chunk.hap,"index")[[i]][j]]<-i ## insert haplotype number into array a at appropriate index
        } 
      }
      
      ## then put genotype info into matrix m
      m[,1]<-a[seq(1,(nrow(chunk.dna))-1,2)] ## odd indices, allele A
      m[,2]<-a[seq(2,(nrow(chunk.dna)),2)] ## even indices, allele B
      m<-t(apply(m,1,sort)) ## sort so that lower haplotype number always first
      
      ## make table that sums all genotypes
      tab<-matrix(0,nrow=nrow(chunk.hap),ncol=nrow(chunk.hap))
      for (i in (1:nrow(m))) {
        tab[m[i,1],m[i,2]]<-tab[m[i,1],m[i,2]]+1
      }
      
      ## generate all possible combinations of 3 haplotypes
      ## cycles has rows = 1:3, haplotype number, cols = all possible combo's 
      cycles<-combn(1:nrow(chunk.hap),3)
      colnames(cycles)<-rep("FALSE",ncol(cycles))
      for (j in (1:ncol(cycles))) {
        count<-0
        if (tab[cycles[1,j],cycles[2,j]]>0) count<-count+1
        if (tab[cycles[1,j],cycles[3,j]]>0) count<-count+1
        if (tab[cycles[2,j],cycles[3,j]]>0) count<-count+1
        if (count>=3) {colnames(cycles)[j]<-"TRUE"} else {colnames(cycles)[j]<-"FALSE"}  
      }
      
      triplets<-sum(as.logical(colnames(cycles)))
      if (triplets>0) {
        windows.triplets[window.index]<-1
        if (ends[window.index]>length(dna[1,])) { coords.end<-length(dna[1,]) } else { coords.end<-ends[window.index] } ## make end coordinate make sense RE length of sequence
        if (opt$verbose>1) cat("Triplet found in window ",window.index,", coords ",sep="")
        if (opt$verbose>1) cat(starts[window.index],":",coords.end,"\n",sep="")
        
        if (opt$verbose>1) {
          for (col in (1:ncol(cycles))) {
            if (colnames(cycles)[col]==TRUE) {
              alleles<-t(combn(cycles[,col],2)) ## this is a genotype table of triplet, showing which hap genotypes are contributing to the triplet
              cat("\tAlleles: ")
              cat(cycles[,col],sep=", ")
              cat("; found in:\n")
              for (i in 1:nrow(m)) {
                for (j in 1:nrow(alleles)) {
                  if(m[i,1]==alleles[j,1] & m[i,2]==alleles[j,2]){
                    cat("\t\tInd_",i," (genotype ",sep="")
                    cat(m[i,],sep="|")
                    cat(")\n",sep="")
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  ## get longest sequential run of haplotype trio windows
  if (sum(windows.triplets)>0) {
    cat("\nHaplotype trio detected!\n")
    cat("Number of windows with haplotype trio:",sum(windows.triplets),"\n")
    r<-rle(windows.triplets)
    runs<-r[[1]][r[[2]]==1]
    total<-sum(opt$window+(runs-opt$step))
    cat("Total haptrio sites (with current window):",total,"\n")
    cat("Length of longest congruent haptrio (with current window): ",max(opt$window+(runs-opt$step))," (sites)\n",sep="")
  } else {
    cat("\nNo haplotype trios detected\n")
  }
  
} else if (opt$decreasing==TRUE) {
  
  ###################################
  ## Do decreasing window analysis ##
  ###################################
  
  for (window in (length(dna[1,]):1)) {
    
    starts<-seq(1,length(dna[1,]),by=opt$step) ## gives start coordinates for all windows
    ends<-starts+(window-1) ## end coords for all windows
    starts<-starts[ends<=length(dna[1,])] ## trims window coordinates to fit within seq coords
    ends<-ends[ends<=length(dna[1,])]
    if (opt$verbose>0) {
      cat("\nWindow size:",window,"\n")
      cat("Number of windows:",length(starts))
    }
    
    windows.triplets<-rep(0,length(starts)) ## initialise a vector of length == number windows
    
    for (window.index in (1:length(starts))) {
      ## get chunk and make hap from it
      chunk.dna<-dna[,starts[window.index]:ends[window.index]]
      if (opt$verbose>1) cat(starts[window.index],":",ends[window.index],"\n",sep="")
      chunk.hap<-haplotype(chunk.dna)
      
      if (nrow(chunk.hap)>=3) { ## no point testing windows with <3 haplotypes
        a<-array(NA,nrow(chunk.dna)) ## generate NA array equal to number of alleles in sample
        m<-matrix(NA,nrow=(nrow(chunk.dna)/2),ncol=2) ## generate NA matrix where rows == individual cols == genotype
        
        for (i in 1:(length(attr(chunk.hap,"index")))) { 
          for (j in 1:(length(attr(chunk.hap,"index")[[i]]))) { 
            a[attr(chunk.hap,"index")[[i]][j]]<-i ## insert haplotype number into array a at appropriate index
          } 
        }
        
        ## then put genotype info into matrix m
        m[,1]<-a[seq(1,(nrow(chunk.dna))-1,2)] ## odd indices, allele A
        m[,2]<-a[seq(2,(nrow(chunk.dna)),2)] ## even indices, allele B
        m<-t(apply(m,1,sort)) ## sort so that lower haplotype number always first
        
        ## make table that sums all genotypes
        tab<-matrix(0,nrow=nrow(chunk.hap),ncol=nrow(chunk.hap))
        for (i in (1:nrow(m))) {
          tab[m[i,1],m[i,2]]<-tab[m[i,1],m[i,2]]+1
        }
        
        ## generate all possible combinations of 3 haplotypes
        ## cycles has rows = 1:3, haplotype number, cols = all possible combo's 
        cycles<-combn(1:nrow(chunk.hap),3)
        colnames(cycles)<-rep("FALSE",ncol(cycles))
        for (j in (1:ncol(cycles))) {
          count<-0
          if (tab[cycles[1,j],cycles[2,j]]>0) count<-count+1
          if (tab[cycles[1,j],cycles[3,j]]>0) count<-count+1
          if (tab[cycles[2,j],cycles[3,j]]>0) count<-count+1
          if (count>=3) {colnames(cycles)[j]<-"TRUE"} else {colnames(cycles)[j]<-"FALSE"}  
        }
        
        triplets<-sum(as.logical(colnames(cycles)))
        if (triplets>0) {
          windows.triplets[window.index]<-1
          if (ends[window.index]>length(dna[1,])) { coords.end<-length(dna[1,]) } else { coords.end<-ends[window.index] } ## make end coordinate make sense RE length of sequence
          if (opt$verbose>1) cat("Triplet found in window ",window.index,", coords ",sep="")
          if (opt$verbose>1) cat(starts[window.index],":",coords.end,"\n",sep="")
          
          if (opt$verbose>1) {
            for (col in (1:ncol(cycles))) {
              if (colnames(cycles)[col]==TRUE) {
                alleles<-t(combn(cycles[,col],2)) ## this is a genotype table of triplet, showing which hap genotypes are contributing to the triplet
                cat("\tAlleles: ")
                cat(cycles[,col],sep=", ")
                cat("; found in:\n")
                for (i in 1:nrow(m)) {
                  for (j in 1:nrow(alleles)) {
                    if(m[i,1]==alleles[j,1] & m[i,2]==alleles[j,2]){
                      cat("\t\tInd_",i," (genotype ",sep="")
                      cat(m[i,],sep="|")
                      cat(")\n",sep="")
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    
    ## get longest sequential run of haplotype trio windows
    if (sum(windows.triplets)>0) {
      cat("\nHaplotype trio detected!\n")
      cat("Number of windows with haplotype trio:",sum(windows.triplets),"\n")
      r<-rle(windows.triplets)
      runs<-r[[1]][r[[2]]==1]
      total<-sum(window+(runs-opt$step))
      cat("Total haptrio sites (with current window):",total,"\n")
      cat("Length of longest congruent haptrio (with current window): ",max(opt$window+(runs-opt$step))," (sites)\n",sep="")
    } else {
      cat("\nNo haplotype trios detected\n")
    }
    if (opt$quitonfirst==TRUE && sum(windows.triplets)>0) { break }
  }
}
cat("\n[END]\n\n")
