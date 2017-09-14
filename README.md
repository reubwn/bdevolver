# bdevolver
simulating nucleotide evolution with mutation, drift, gene conversion and horizontal transfer

### libraries
Requires the libraries `optparse` and `pegas`:
```R
install.packages("optparse")
install.packages("pegas")
```

## bdevolver.simulation.R
Main simulation script. Outputs sequences in fasta format.

Type `bdevolver.simulation.R -h` list of options:

```R
Usage: /home/rnowell/scripts/bdevolver.simulation.R [options]


Options:
        -N INTEGER, --ninds=INTEGER
                Number of individuals to simulate [default=2]

        -n INTEGER, --nlocs=INTEGER
                Number of genes per individual [default=2]

        -p INTEGER, --ploidy=INTEGER
                Ploidy level [default=2] NOTE: only diploid supported atm

        -l INTEGER, --genelength=INTEGER
                Gene length [default=100]

        -m DOUBLE, --mutation=DOUBLE
                Probability of mutation occurring at a given nucleotide site, per generation [default=0.01]

        -A DOUBLE, --freqA=DOUBLE
                Expected frequency of A nucleotide [default=0.25]

        -T DOUBLE, --freqT=DOUBLE
                Expected frequency of T nucleotide [default=0.25]

        -G DOUBLE, --freqG=DOUBLE
                Expected frequency of G nucleotide [default=0.25]

        -C DOUBLE, --freqC=DOUBLE
                Expected frequency of C nucleotide [default=0.25]

        -c DOUBLE, --conversion=DOUBLE
                Probability of a gene conversion event initiating at a given nucleotide site, per generation [default=0]

        --lamdaC1=INTEGER
                First mean tractlength for gene conversion [default=30]

        --lamdaC2=INTEGER
                Second mean tractlength for gene conversion if tract lengths are sampled from a bimodal distribution [default=200]

        --propC1=INTEGER
                Proportion of times conversion tractlength should be sampled from the first specified distribution [default=1]

        -t DOUBLE, --transfer=DOUBLE
                Probability of a horizontal transfer event initiating at a given nucleotide site, per generation [default=0]

        --lamdaT1=INTEGER
                First mean tractlength for horizontal transfer [default=30]

        --lamdaT2=INTEGER
                Second mean tractlength for horizontal transfer if tract lengths are sampled from a bimodal distribution [default=200]

        --propT1=INTEGER
                Proportion of times transfer tractlength should be sampled from the first specified distribution [default=1]

        -g INTEGER, --ngens=INTEGER
                Number of generations to evolve sequences for [default=100]

        -r INTEGER, --nreps=INTEGER
                Number of replicate simulations to perform [default=1]

        --noseqs
                Don't print fasta sequences at end of simulation [default=FALSE]

        --nostats
                Don't print any stats per locus [default=FALSE]

        --parentage
                Print parentage information to file [default=TRUE]

        -v INTEGER, --verbose=INTEGER
                Verboseness (0=quiet, 1=verbose, 2=vv, 3=vvv >3=debug) [default=1]

        -h, --help
                Show this help message and exit
```

## bdevolver.triplets.R
Script for the detection of haplotype trios in a sequence alignment. Type `-h` for help.

## bdevolver.plot.R
Script that plots some time-series population genetics from the output of bdevolver.simulation.R. Type `-h` for help.

## bdevolver.stats.R
Outputs some stats in a table format.
