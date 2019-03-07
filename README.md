# HILDA


# Introduction

The R package `HiLDA` is developed under the Bayesian framework to allow statisticaly testing whether there is a change in the mutation burdnes of mutation sigantures between two groups. The mutation signature is defined based on the independent model proposed by Shiraishi's et al. 

# Paper

- Shiraishi et al. A simple model-based approach to inferring and visualizing cancer mutation signatures, bioRxiv, doi: [http://dx.doi.org/10.1101/019901](http://dx.doi.org/10.1101/019901).

# Input data

`HiLDA` is a package built on some basic functions from `pmsignature` including how to read the input data. Here is an example from `pmsignature` on the input data, *mutation features* are elements used for categorizing mutations such as: 
  
* 6 substitutions (C>A, C>G, C>T, T>A, T>C and T>G)
* 2 flanking bases (A, C, G and T)
* transcription direction.

## Mutation Position Format

sample1 chr1  100	A	C	

sample1	chr1	200	A	T	

sample1	chr2	100	G	T	

sample2	chr1	300	T	C	

sample3	chr3	400	T	C	
  
* The 1st column shows the name of samples 
* The 2nd column shows the name of chromosome 
* The 3rd column shows the coordinate in the chromosome
* The 4th column shows the reference base (A, C, G, or T).
* The 5th colum shows the alternate base (A, C, G, or T).


# Workflow 
## Install the packages
First, a few R packages such as `pmsignature`, `gtools`, `R2jags` have to be installed prior to using `HiLDA`. Currently, the easiest way for installing pmsignature is to use the package devtools. 

## Install `JAGS`
Download and install JAGS by following the instructions from http://mcmc-jags.sourceforge.net/

## Get input data
Read in the data by using the function from `pmsignature`. 

```r
G <- pmsignature::readMPFile(inputFile, numBases = 5, trDir = FALSE, 
                bs_genome = BSgenome.Hsapiens.UCSC.hg18::BSgenome.Hsapiens.UCSC.hg18,
                txdb_transcript = TxDb.Hsapiens.UCSC.hg18.knownGene::TxDb.Hsapiens.UCSC.hg18.knownGene)
```

Here, *inputFile* is the path for the input file. *numBases* is the number of flanking bases to consider including the central base (if you want to consider two 5' and 3' bases, then set 5). Also, you can add transcription direction information using *trDir*. *numSig* sets the number of mutation signatures estimated from the input data.  

```r
library(HiLDA)
inputFile <- system.file("data/sampleG.rdata", package = "HiLDA")
load(inputFile)
```

# Get signatures

```r
K <- 3
Param <- pmsignature::getPMSignature(G, K = K)
```

# Visualize the mutation signatures 

```r
hilda_plotSignature(Param)
```

# Generate the initial values and run the global test (Bayes Factor) and the local test 
```r
hilda_gloabl <- hilda_bayesfactor(inputG = G, inputParam = Param, refGroup = seq(1,20,2), n.iter = 2000)
hilda_local <- hilda_test(inputG = G, inputParam = Param, refGroup = seq(1,20,2), n.iter = 2000)
```

# Assess Convergence of MCMC chains
```r
hilda_rhat(hilda_local)
```

# Output the posterior distribution of the mean difference
```r
hilda_bayesfactor_result(hilda_gloabl)
```

# Output the posterior distribution of the mean difference
```r
hilda_posterior(hilda_local)
```


