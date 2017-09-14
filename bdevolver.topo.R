#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(pegas))
suppressPackageStartupMessages(library(phangorn))
suppressPackageStartupMessages(library(optparse))
option_list = list(
  make_option(c("-i","--trees1"), type="character", default=NULL, help="Input multitree file in Newick format [REQUIRED]"),
  make_option(c("-p","--prefix"), type="character", default="out", help="Output filename written to ${PREFIX}.treedistcsv [default=%default]"),
  make_option(c("-s","--trees2"), type="character", default=NULL, help="Second multitree file in Newick format; if none provided then comparisons within trees1 are computed"),
  make_option(c("-v","--verbose"), type="integer", default=0, help="Verboseness (0=quiet, 1=verbose) [default=%default]")
); 
opt_parser=OptionParser(option_list=option_list);
opt=parse_args(opt_parser);

n<-0
df<-data.frame()
trees1<-read.tree(opt$trees1)
if (!(is.null(opt$trees2))) trees2<-read.tree(opt$trees2)
#trees1<-read.tree("mutation/trees/mut.phylo")
#trees2<-read.tree("conversion/trees/con.phylo")

if (!(exists("trees2"))) { ## do comparisons within trees1
  for (i in (1:length(trees1))) {
    for (j in (length(trees1):(1+n))) {
      if (i!=j) {
        if (opt$verbose > 0) cat(i,"vs",j,"\n")
#        if (opt$verbose > 0) print(treedist(trees1[[i]],trees1[[j]],check.labels=F))
        df<-rbind(df,cbind(paste(i,":",j,sep=""),rbind(treedist(trees1[[i]],trees1[[j]],check.labels=F))))
      }
    }
    n<-n+1
  }
} else { ## compare trees1 vs trees2
  for (i in (1:length(trees1))) {
    for (j in (length(trees2):1)) {
      if (opt$verbose > 0) cat(i,"vs",j,"\n")
#      if (opt$verbose > 0) print(treedist(trees1[[i]],trees1[[j]],check.labels=F))
      df<-rbind(df,cbind(paste(i,":",j,sep=""),rbind(treedist(trees1[[i]],trees1[[j]],check.labels=F))))
    }
  }
}
colnames(df)<-c("comp","symmetric.difference","branch.score.difference","path.difference","quadratic.path.difference")
filename<-paste(opt$prefix,".treedist.csv",sep="")
write.csv(df,file=filename,row.names=F)