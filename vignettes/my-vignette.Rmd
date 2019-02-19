---
title: "HiLDA: a package for testing the burdens of mutational signatures"
author: 
- name: Zhi Yang
  affiliation:
  - Department of Preventive Medicine, University of Southern California, Los Angeles, USA 
date: "`r Sys.Date()`"
vignette: |
  %\VignetteIndexEntry{An introduction to HiLDA}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
output:
  BiocStyle::html_document:
    toc_float: true
  BiocStyle::pdf_document: default
abstract: | 
  Instructions on using _HiLDA_ on testing the burdens of mutational signatures. 

---

```{r style, echo = FALSE, results = 'asis'}
library(BiocStyle)
```

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
## Install the package
First, a few R packages such as `pmsignature`, `gtools`, `R2jags` have to be installed prior to using `HiLDA`. Currently, the easiest way for installing pmsignature is to use the package devtools. 

## Install `JAGS`
Download and install JAGS by following the instructions from http://mcmc-jags.sourceforge.net/

## Get input data
Read in the data by using the function from `pmsignature`. 

```{r eval=FALSE}
G <- pmsignature::readMPFile(inputFile, numBases = 5, trDir = FALSE, 
                bs_genome = BSgenome.Hsapiens.UCSC.hg18::BSgenome.Hsapiens.UCSC.hg18,
                txdb_transcript = TxDb.Hsapiens.UCSC.hg18.knownGene::TxDb.Hsapiens.UCSC.hg18.knownGene)
```

Here, *inputFile* is the path for the input file. *numBases* is the number of flanking bases to consider including the central base (if you want to consider two 5' and 3' bases, then set 5). Also, you can add transcription direction information using *trDir*. *numSig* sets the number of mutation signatures estimated from the input data.  

```{r eval=FALSE}
library(HiLDA)
inputFile <- system.file("data/sampleG.rdata", package = "HiLDA")
```

# Get signatures

```{r eval=FALSE}
K <- 3
Param <- pmsignature::getPMSignature(G, K = K)
```

# Visualize the mutation signatures 

```{r eval=FALSE}
hilda.plotSignature(Param)
```


# Generate the initial values and run HiLDA test
```{r eval=FALSE}
hilda.result <- hilda.test(G, Param, refGroup = seq(1,20,2), n.iter = 2000)
```

# Assess Convergence of MCMC chains
```{r eval=FALSE}
hilda.rhat(hilda.result)
```

# Output the posterior distribution of the mean difference in mutational exposures
```{r eval=FALSE}
hilda.rhat(hilda.test)
```


