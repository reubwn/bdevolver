################################################################################
## bdelloid_evolver #
#####################
## Script for simulating forward-in-time sequences under asexuality + mutation + drift + gene conversion + HGT

################################################################################
## get options #
################

## load option parsing lib
suppressPackageStartupMessages(library(optparse))
option_list = list(
	make_option(c('-N','--ninds'), type='integer', default=2, help='Number of individuals to simulate [default=%default]', metavar='integer'),
	make_option(c('-n','--nlocs'), type='integer', default=2, help='Number of genes per individual [default=%default]',metavar='integer'),
	make_option(c('-p','--ploidy'),type='integer', default=2,help='Ploidy level [default=%default] NOTE: only diploid supported atm',metavar='integer'),
	make_option(c('-l','--genelength'), type='integer', default=100, help='Gene length [default=%default]', metavar='integer'),
	make_option(c('-m','--mutation'), type='double', default=0.01, help='Probability of mutation occurring at a given nucleotide site, per generation [default=%default]', metavar='double'),
	make_option(c('-A','--freqA'), type='double', default=0.25, help='Expected frequency of A nucleotide [default=%default]', metavar='double'),
	make_option(c('-T','--freqT'), type='double', default=0.25, help='Expected frequency of T nucleotide [default=%default]', metavar='double'),
	make_option(c('-G','--freqG'), type='double', default=0.25, help='Expected frequency of G nucleotide [default=%default]', metavar='double'),
	make_option(c('-C','--freqC'), type='double', default=0.25, help='Expected frequency of C nucleotide [default=%default]', metavar='double'),
	make_option(c('-c','--conversion'), type='double', default=0.00, help='Probability of a gene conversion event initiating at a given nucleotide site, per generation [default=%default]', metavar='double'),
	make_option(c('--lamdaC1'),type='integer', default=30, help='First mean tractlength for gene conversion [default=%default]', metavar='integer'),
	make_option(c('--lamdaC2'),type='integer', default=200, help='Second mean tractlength for gene conversion if tract lengths are sampled from a bimodal distribution [default=%default]', metavar='integer'),
	make_option(c('--propC1'),type='double', default=1, help='Proportion of times conversion tractlength should be sampled from the first specified distribution [default=%default]', metavar='integer'),
	make_option(c('-t','--transfer'), type='double', default=0.00, help='Probability of a horizontal transfer event initiating at a given nucleotide site, per generation [default=%default]', metavar='double'),
	make_option(c('--lamdaT1'),type='integer', default=30, help='First mean tractlength for horizontal transfer [default=%default]', metavar='integer'),
	make_option(c('--lamdaT2'),type='integer', default=200, help='Second mean tractlength for horizontal transfer if tract lengths are sampled from a bimodal distribution [default=%default]', metavar='integer'),
	make_option(c('--propT1'),type='double', default=1, help='Proportion of times transfer tractlength should be sampled from the first specified distribution [default=%default]', metavar='integer'),
	make_option(c('-g','--ngens'), type='integer', default=100, help='Number of generations to evolve sequences for [default=%default]', metavar='integer'), ## throws weird bug with -g
	make_option(c('-r','--nreps'), type='integer', default=1, help='Number of replicate simulations to perform [default=%default]', metavar='integer'),
	make_option(c('--noseqs'), action='store_true', default=FALSE, help='Don\'t print fasta sequences at end of simulation [default=%default]'),
	make_option(c('--nostats'), action='store_true', default=FALSE, help='Don\'t print any stats per locus [default=%default]'),
	make_option(c('--parentage'), action='store_false', default=TRUE, help='Print parentage information to file [default=%default]'),
	make_option(c('-v','--verbose'), type='integer', default=1, help='Verboseness (0=quiet, 1=verbose, 2=vv, 3=vvv >3=debug) [default=%default]', metavar = 'integer')
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (opt$verbose>2) {
	library(pegas)
} else {
	suppressPackageStartupMessages(library(pegas))
}

## hard-code options for testing purposes ##
# opt$ninds<-20
# opt$ninds<-100
# opt$nlocs<-2
# opt$genelength<-100
# opt$mutation<-0.04
# opt$tractlength<-5
# opt$transfer<-0
# opt$ngens<-100
# opt$verbose<-1
# opt$parentage<-F

## evaluate parameters and convert if necessary
if(is.numeric(opt$mutation)){ mutation<-opt$mutation } else { mutation<-eval(parse(text=opt$mutation)) } ## mutation
if(is.numeric(opt$conversion)){ conversion<-opt$conversion } else { conversion<-eval(parse(text=opt$conversion)) } ## gene conversion
if(is.numeric(opt$transfer)){ transfer<-opt$transfer } else { transfer<-eval(parse(text=opt$transfer)) } ## HGT
if(is.numeric(opt$ngens)){ ngens<-opt$ngens } else { ngens<-eval(parse(text=opt$ngens)) } ## ngens
if(is.numeric(opt$genelength)){ genelength<-opt$genelength } else { genelength<-eval(parse(text=opt$genelength)) } ## genelength

## because 1 conversion/transfer event may affect multiple nucleotides, useful to calculate the
## per-nuc Pr to be converted/transferred by multiplying rates by average tractlengths
prob_conversion<-conversion*((opt$lamdaC1*opt$propC1)+(opt$lamdaC2*(1-opt$propC1)))
prob_transfer<-transfer*((opt$lamdaT1*opt$propT1)+(opt$lamdaT2*(1-opt$propT1)))
conversion.lengths<-array(NA,opt$ngens*opt$nlocs)
transfer.lengths<-array(NA,opt$ngens*opt$nlocs)

###################
## run simulation #
###################
cat('\n[START]\n\n')
for (reps in (1:opt$nreps)) { ## iterate thru nreps
	if (opt$verbose>0) cat('Replicate:',reps,'\n')
	for (m in (1:length(mutation))) { ## iterate thru values of mutation
		if (opt$verbose>0) cat('Mutation rate:',format(mutation[m],scientific=T,digits=4),'\n')
		for (c in (1:length(conversion))) { ## iterate thru values of rho.intra
			if (opt$verbose>0) cat('Conversion rate:',format(conversion[c],scientific=T,digits=4),'\n')
			if (opt$verbose>0) cat('Mean conversion tractlength:',format(((opt$lamdaC1*opt$propC1)+(opt$lamdaC2*(1-opt$propC1))),scientific=F,digits=4),'\n')
			if (opt$verbose>0) cat('Probability of nucleotide to be converted:',format(prob_conversion,scientific=T,digits=4),'\n')
			for (h in (1:length(transfer))) { ## iterate thru values of rho.inter
				if (opt$verbose>0) cat('Transfer rate:',format(transfer[h],scientific=T,digits=4),'\n')
				if (opt$verbose>0) cat('Mean transfer tractlength:',format(((opt$lamdaT1*opt$propT1)+(opt$lamdaT2*(1-opt$propT1))),scientific=F,digits=4),'\n')
				if (opt$verbose>0) cat('Probability of nucleotide to be transferred:',format(prob_transfer,scientific=T,digits=4),'\n')

				## make population of initial, identical sequences
				seq.matrix<-matrix(sample(c('A','T','G','C'),genelength,replace=T,prob=c(opt$freqA,opt$freqT,opt$freqG,opt$freqC)),nrow=opt$ploidy*opt$ninds,ncol=genelength,byrow=T) ## rows are haplotypes; columns are sequence position
				locus.array<-array(seq.matrix,dim=c(opt$ploidy*opt$ninds,genelength,opt$nlocs))
				dimnames(locus.array)[[1]]<-paste("ind_",rep(1:opt$ninds,each=2),'',rep(c('A','B'),opt$ninds),sep='')
				dimnames(locus.array)[[3]]<-paste("locus_",1:opt$nlocs,sep='')
				## NB
				## The object locus.array is an array of nucleotide matrices
				## so locus.array[,,,] prints all nucs for all loci
				## locus.array[1,,2] prints 1st haplotype for locus 2
				## locus.array[1,10,2] prints 10th nuc of 1st hap for locus 2 etc.

#				conversion.lengths<-array(NA,opt$ngens*opt$nlocs)
#				transfer.lengths<-array(NA,opt$ngens*opt$nlocs)
#				sample<-array(NA,opt$ngens)

				## begin ngen loop
				for (gen in (1:ngens)) {
					if (opt$verbose>0) cat('\tGeneration:',gen,'\n')

					## generate stats filename per locus
					filename.stats.list<-array(NA,opt$nlocs)

					for (loc in (1:opt$nlocs)) {
						if (opt$verbose>1) cat('\nLocus:',loc,'\n')

						###################
						## apply MUTATION #
						###################
						muts<-runif(opt$ninds*genelength*opt$ploidy)<mutation[m] ## logical vector of length opt$ninds*l*ploidy, dictating which base in that vector will be mutated
#						print(sum(muts))
#						sample[gen]<-sum(muts)
						if (sum(muts>0)) {
							if (opt$verbose>1) {
								mutation.tab<-which(matrix(muts,ncol=genelength),arr.ind=TRUE)
								mutation.tab<-cbind(mutation.tab,locus.array[,,loc][muts])
							}

							## do the mutation
							locus.array[,,loc][muts]<-sample(c('A','T','G','C'),sum(muts),replace=T,prob=c(opt$freqA,opt$freqT,opt$freqG,opt$freqC))

							if (opt$verbose>1) {
								mutation.tab<-cbind(mutation.tab,locus.array[,,loc][muts])
								colnames(mutation.tab)<-c('hap','base','from','to')
								cat('Mutation table:\n')
								print(data.frame(mutation.tab))
							}
						}

						##########################
						## apply GENE-CONVERSION #
						##########################
						## (ie intra-individual 'recombination')
						## this should in effect replace a tract within allele 1 of a given individual with the homologous tract from
						## the alternative allele
						## tract lenth sampled from Poisson with mean of opt$tractlength

						## rhos.intra is a logical vector of length equal to the total number of base-pairs across all individuals in population,
						## returning 'TRUE' only when random deviate from runif is < given rho-inter value from --conversions
						rhos.intra<-runif(opt$ninds*genelength*opt$ploidy)<conversion[c]
						if (sum(rhos.intra)>0) {

							## make indices of conversion event(s) (ie, which nucleotide should GC initiate from)
							conversion.tab<-which(matrix(rhos.intra,ncol=genelength),arr.ind=TRUE)
							## NB
							## conversion.tab returns something like:
							##       row col
							## [1,]    6   8
							## indicating the conversion should start at the 8th nucleotide of haplotype 6

							## iterate across conversion events
							for (r in (1:nrow(conversion.tab))) {
								## sample tract length from Poisson with mean opt$lamdaC1
							  ## allows for sampling from a bimodal distribution if opt$propC1 is set to <1
							  x<-runif(1)
							  if (x<=opt$propC1) {
							    conversion.tractlength<-rpois(1,opt$lamdaC1)
							    conversion.lengths[gen]<-conversion.tractlength
							  } else {
							    conversion.tractlength<-rpois(1,opt$lamdaC2)
							    conversion.lengths[gen]<-conversion.tractlength
							  }

								## get start-end coordinates for conversion
								conversion.start<-conversion.tab[r,2]
								if ((conversion.tab[r,2]+conversion.tractlength)>genelength) { ## if tractlength exceeds gene length, take the tract to the end of the gene
									conversion.end<-genelength
								} else {
									conversion.end<-conversion.tab[r,2]+(conversion.tractlength-1) ## -1 because of 1-offset indexing
								}

								## get donor-recipient info
								conversion.recipient<-locus.array[conversion.tab[r,1],,loc]
								conversion.recipient.name<-dimnames(locus.array)[[1]][conversion.tab[r,1]]
								if (conversion.tab[r,1] %% 2 == 0) { ## if recipient allele is even, donor tract should be from previous allele in locus.array
									conversion.donor<-locus.array[conversion.tab[r,1]-1,,loc]
									conversion.donor.name<-dimnames(locus.array)[[1]][conversion.tab[r,1]-1]
								} else { ## else if recipient is odd, donor should be from next allele
									conversion.donor<-locus.array[conversion.tab[r,1]+1,,loc]
									conversion.donor.name<-dimnames(locus.array)[[1]][conversion.tab[r,1]+1]
								}

								## do the conversion
								locus.array[conversion.recipient.name,conversion.start:conversion.end,loc]<-locus.array[conversion.donor.name,conversion.start:conversion.end,loc]

								## die if donor tract != recipient tract
								if (paste(locus.array[conversion.recipient.name,conversion.start:conversion.end,loc],sep='',collapse='') != paste(locus.array[conversion.donor.name,conversion.start:conversion.end,loc],sep='',collapse='')) {
									stop('\nDonor sequence not equal to recipient sequence after conversion...?\n\n',call. = F)
								}

								if (opt$verbose>1) {
									cat('Conversion event:','[donor ',conversion.donor.name,'][recipient ',conversion.recipient.name,'][length ',conversion.tractlength,'][coords ',conversion.start,'::',conversion.end,']\n',sep='')
								}
							}## end of gene conversion events loop
						}

						###############################
						## apply HORIZONTAL TRANSFER ##
						###############################
						## models between-individual recombination, BIR, as the non-reciprocal transfer of an allele to another individual
						## BIR is still occurring within the species boundary, so assume that horizontally transferred alleles are
						## integrating into the genome via HR
						## could add an option that ensures that the HGT transfer of an allele is always immediately followed by the conversion of the other allele, such that the HGT donor becomes homozygous for the transferred allele...?

						## rhos.inter is a logical vector of length equal to the total number of base-pairs across all individuals in population,
						## returning 'TRUE' only when random deviate from runif is < given rho-inter value from --hgt
						rhos.inter<-runif(opt$ninds*genelength*opt$ploidy)<transfer[h]
						if (sum(rhos.inter)>0) {
							transfer.tab<-which(matrix(rhos.inter,ncol=genelength),arr.ind=TRUE) ## table indicating index of recipient and position in gene where transfer begins

							for (r in (1:nrow(transfer.tab))) {
								transfer.donor.index<-sample(1:(opt$ninds*opt$ploidy),1) ## choose index of allele to be transferred

								## get tractlength of transferred element
								if (runif(1)<=opt$propT1) {
								  transfer.tractlength<-rpois(1,opt$lamdaT1)
								  transfer.lengths[gen]<-transfer.tractlength
								} else {
								  transfer.tractlength<-rpois(1,opt$lamdaT2)
								  transfer.lengths[gen]<-transfer.tractlength
								}
								print(transfer.tractlength)
								## get transfer coordinates
								transfer.start<-transfer.tab[r,2]
								if ((transfer.start+transfer.tractlength)>genelength) { ## if tractlength exceeds gene length, take the tract to the end of the gene
									transfer.end<-genelength
								} else {
									transfer.end<-transfer.tab[r,2]+(transfer.tractlength-1) ## -1 because of 1-offset indexing
								}

								## do the transfer!
								## AND ensure transfers to self are not allowed...
								if (transfer.tab[r,1]  %% 2 == 0) { ## if recipient allele index is even number
									if (transfer.tab[r,1] != transfer.donor.index && transfer.tab[r,1] != transfer.donor.index-1) {
										locus.array[transfer.tab[r,1],transfer.start:transfer.end,loc]<-locus.array[transfer.donor.index,transfer.start:transfer.end,loc]

										if (opt$verbose>1) {
											transfer.donor.name<-dimnames(locus.array)[[1]][transfer.donor.index]
											transfer.recipient.name<-dimnames(locus.array)[[1]][transfer.tab[r,1]]
											cat('transfer event:','[donor ',transfer.donor.name,'][recipient ',transfer.recipient.name,'][length ',transfer.tractlength,'][coords ',transfer.start,'::',transfer.end,']\n',sep='')
										}
									}
								} else { ## recipient allele index is odd number
									if (transfer.tab[r,1] != transfer.donor.index && transfer.tab[r,1] != transfer.donor.index+1) {
										locus.array[transfer.tab[r,1],transfer.start:transfer.end,loc]<-locus.array[transfer.donor.index,transfer.start:transfer.end,loc]

										if (opt$verbose>1) {
											transfer.donor.name<-dimnames(locus.array)[[1]][transfer.donor.index]
											transfer.recipient.name<-dimnames(locus.array)[[1]][transfer.tab[r,1]]
											cat('transfer event:','[donor ',transfer.donor.name,'][recipient ',transfer.recipient.name,'][length ',transfer.tractlength,'][coords ',transfer.start,'::',transfer.end,']\n',sep='')
										}
									}
								}
							}## end of transfer events loop
						}

						#################
						## output stats #
						#################
						## calculate metrics
						theta<-theta.s(length(seg.sites(as.DNAbin(tolower(locus.array[,,loc])))),nrow(as.DNAbin(locus.array[,,loc])))
						pi<-nuc.div(as.DNAbin(tolower(locus.array[,,loc])))
						GC<-GC.content(as.DNAbin(tolower(locus.array[,,loc])))
						hap<-haplotype(as.DNAbin(tolower(locus.array[,,loc])))

						## output
						filename.stats<-paste('bdevolver.','reps',format(reps),'.inds',opt$ninds,'.gens',format(ngens,scientific=T),'.mut',format(mutation[m],scientific=T,digits=4),'.con',format(conversion[c],scientific=T,digits=4),'.transfer',format(transfer[h],scientific=T,digits=4),'.loc',loc,'.out.stats',sep='')
						cat(loc,gen,theta,pi,GC,nrow(hap),'\n',sep='\t',file=filename.stats,append=TRUE)
						filename.stats.list[loc]<-filename.stats
					}## end of nloc loop

					################
					## apply DRIFT #
					################
					## random sampling of parents for next generation
					## select with replacement N individuals (opt$ninds*ploidy haplotypes) from locus.array
					if (gen<ngens) { ## don't do it on the final iteration, otherwise stats and seqs will have different final values
						indices<-sample(seq(from=1,to=(opt$ninds*opt$ploidy),by=2),size=opt$ninds,replace=T) ## randomly select ODD NUMBERED ONLY indices from 1 to number alleles; this will always be the 'A' allele of any given ind_*
						selected<-array(NA,opt$ninds*opt$ploidy) ## create dummy array
						selected[seq(1,(opt$ninds*opt$ploidy)-1,2)]<-indices ## insert odd-numbered indices into selected
						selected[seq(2,opt$ninds*opt$ploidy,2)]<-indices+1 ## for each A allele, insert corresponding B allele (always indice[i]+1)
						locus.array<-locus.array[selected,,] ## create new locus.array from selected parents
					}

					## rename alleles at end of each generation
					dimnames(locus.array)[[1]]<-paste("ind_",rep(1:opt$ninds,each=2),'',rep(c('A','B'),opt$ninds),sep='')

					## print parentage info if --parentage
					if (opt$parentage == FALSE) {
						filename.parents<-paste('bdevolver.','rep',format(reps),'.mut',format(mutation[m],scientific=T,digits=4),'.con',format(conversion[c],scientific=T,digits=4),'.transfer',format(transfer[h],scientific=T,digits=4),'.out.parents',sep='')
						cat(gen,dimnames(locus.array)[[1]],'\n',sep='\t',file=filename.parents,append=T)
					}
				}## end of ngen loop

				################
				## output seqs #
				################
				if(!(opt$noseqs)) {
					for (loc in (1:length(dimnames(locus.array)[[3]]))) { ## print one file per locus
						filename.fasta<-paste('bdevolver.','reps',format(reps),'.inds',opt$ninds,'.gens',format(ngens,scientific=T),'.mut',format(mutation[m],scientific=T,digits=4),'.con',format(conversion[c],scientific=T,digits=4),'.transfer',format(transfer[h],scientific=T,digits=4),'.loc',loc,'.out.fa',sep='')
						write.dna(locus.array[,,loc],file=filename.fasta,format="fasta") ## write seqs
						if (opt$verbose>0) cat('\nSeqs written to file:',filename.fasta) ## print info
					}
				}

				## if parents is specified
				if (opt$parentage==FALSE) cat('\n\tParents written to file:',filename.parents)
				## prepend header to stats file
				for (fii in (1:length(filename.stats.list))) {
					fi<-file(filename.stats.list[fii],"r+")
					lines<-readLines(fi)
					writeLines(c("loc\tgen\ttheta.s\tnuc.div\tgc\thaplo",lines),con=fi)
					close(fi)
				}

				## print hist of conversion lengths
				if (opt$conversion>0 && sum(conversion.lengths,na.rm=T)>0) {
				  svg("conversion.tractlengths.svg")
				  hist(conversion.lengths,breaks=(length(conversion.lengths)/10),col="darkred",border="darkred",xlab="tractlength (bp)")
				  garbage<-dev.off()
				}
				if (opt$transfer>0 && sum(transfer.lengths,na.rm=T)>0) {
				  svg("transfer.tractlengths.svg")
				  hist(transfer.lengths,breaks=(length(transfer.lengths)/10),col="darkblue",border="darkblue",xlab="tractlength (bp)")
				  garbage<-dev.off()
				}
			}## end of rho.inter loop
		}## end of rho.intra loop
	}## end of mutation loop
}## end of nreps loop
cat('\n\n[END]\n\n')
