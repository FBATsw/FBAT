# Source Code
This repository contains the Family Based Association Test (FBAT) software source code. The FBAT software provides the implementation for several family-based association tests. 

## Update January 2025
The new code release contains methodology to perform statistical fine-mapping of FBAT association statistics in a genetic region of interest. This code is under development and will be updated accordingly.

## FBAT-TOOLKIT
The FBAT-Toolkit was developed in the Department of Biostatistics at the Harvard T.H. Chan School of Public Health. The FBAT software (Linux executable, 64 bit, region-based extension "gen_rv", development version) and example pedigree files can be found [here](https://sites.google.com/view/fbatwebpage). The user manual can be found [here](https://drive.google.com/file/d/1QYada0wegEbspwFPRlyv7g9hKNv7krmA/view).

## FBAT Overview
FBAT is an acronym for Family-Based Association Tests in genetic analyses. Family-based association designs, as opposed to case-control study designs, are particularly attractive, since they test for linkage as well as association, avoid spurious associations caused by admixture of populations, and are convenient for investigators interested in refining linkage findings in family samples.

The unified approach to family-based tests of association, introduced by Rabinowitz and Laird (2000) and Laird et al. (2000), builds on the original TDT method (Spielman et al., 1993) in which alleles transmitted to affected offspring are compared with the expected distribution of alleles among offspring. In particular, the method puts tests of different genetic models, tests of different sampling designs, tests involving different disease phenotypes, tests with missing parents, and tests of different null hypotheses, all in the same framework. Similar in spirit to a classical TDT test, the approach compares the genotype distribution observed in the ‘cases’ to its expected distribution under the null hypothesis, the null hypothesis being “no linkage and no association” or “no association, in the presence of linkage”. Here, the expected distribution is derived using Mendel’s law of segregation and conditioning on the sufficient statistics for any nuisance parameters under the null. Since conditioning eliminates all nuisance parameters, the technique avoids confounding due to model misspecification as well as admixture or population stratification (Rabinowitz and Laird, 2000; Lazzeroni and Lange, 1998).

In order to adapt these “classical” family-based association tests to even more complex scenarios such as multivariate or longitudinal data sources with either binary or quantitative traits, a broader class of conditional tests has been defined (refer to Laird and Lange, 2006). The FBAT rare variant statistic is the latest addition to these tests (De et al, 2011, Yip et al, 2011, Zhou et al, 2013).

## Implementation of Methods
These methods have all been implemented in the FBAT-toolkit, which consists of two packages: FBAT and PBAT. The software provides methods for a wide range of situations that arise in family-based association studies. It provides options to test linkage or association in the presence of linkage, using marker or haplotype data, single or multiple traits. PBAT can compute a variety of univariate, multivariate and time-to-onset statistics for nuclear families as well as for extended pedigrees. PBAT can also include covariates and gene/covariate-interactions in all computed FBAT-statistics. Further, PBAT can be used for pre- and post-study power calculations and construction of the most powerful test statistic. For situations in which multiple traits and markers are given, PBAT provides screening tools to sift through a large pool of traits and markers and to select the most ‘promising’ combination of traits and markers thereof, while at the same time handling the multiple testing problem. For further details on PBAT, see the PBAT webpage; the remainder of this manual will focus on the FBAT package.



