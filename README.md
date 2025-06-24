# Benchmarking of methods to analyze data derived from GBS-MeDIP

This repository contains scripts used for benchmarking various statistical and bioinformatic methods specifically designed for analyzing data obtained through Genotyping-by-Sequencing combined with Methylated DNA Immunoprecipitation (GBS-MeDIP).

## Overview

GBS-MeDIP is a cost-effective method developed to investigate DNA methylation. Due to its unique data characteristics, traditional methods for RNA-seq or MeDIP-seq analyses are not suitable. This repository provides scripts used to create a count matrix and to perform differential methylated analysis specifically for GBS-MeDIP data.

## Content Description

#### 1. Count Matrix Generation

-  `Preprocessing` Folder containing scripts for merging genomic coordinates and generating SAF files used in subsequent counting steps.

-  `MEDIPS_vs_featureCounts` Folder containing scripts R and shell scripts utilizing featureCounts software and MEDIPS R package to generate count matrices from SAF files and the script used for the comparision of both methods.

#### 2. Distribution Analysis

-  `Distribution_analysis` Folder containing an R script to evaluate the statistical distribution (normal, poisson, negative-binomial) fitting the GBS-MeDIP count data. Calculates and visualizes Akaike Information Criterion (AIC) values for the different tested distributions.

#### 3. Statistical Method Benchmarking

-  `Null_distribution` Folder containing scripts that evaluates the uniformity of p-value distributions to detect biases in the various statistical methods described in this article. In those scripts the reader will find for each of the evaluated method code to compute the False Positive Rate (FPR).

-  `Simulation` Folder containing script which simulates GBS-MeDIP data to compute True Positive Rates (TPR) and generates Receiver Operating Characteristic (ROC) curves to compare the accuracy of different methods.

## Recommended Pipeline

Based on our benchmarking results, we recommend the following pipeline for GBS-MeDIP data:

-  Count Matrix Generation: featureCounts

-  Differential Methylation Analysis: Mann-Whitney test

## Dependencies

R packages: ggplot2 v.3.5.1, fitdistrplus v.1.1-11, EdgeR v4.6.1, Limma v3.64.0, MEDIPS v1.52.0

Tools: featureCounts v2.0.3, samtools v.1.14

## Usage

Clone the repository and follow the script-specific instructions provided as comments within each script.
```
git clone https://github.com/Violeta-de-Anca/Benchmarking-of-methods-to-analyse-data-derived-from-GBS-MeDIP.git
```
## References

Detailed methodology, results, and discussion are presented in our published paper:

Violeta de Anca Prado, et al. "Benchmarking of methods to analyze data derived from GBS-MeDIP." (Details and DOI/link to be added upon publication).

## Authors

-  [Violeta de Anca Prado](https://orcid.org/0000-0003-1845-509X)

-  [Fábio Pértille](https://orcid.org/0000-0002-7214-9184)

-  [Pedro Sá](https://orcid.org/0000-0002-1588-6778)

-  [Marta Gòdia](https://orcid.org/0000-0002-0439-4014)

-  [Joëlle Rüegg](https://orcid.org/0000-0002-6580-9201)

-  [Carlos Guerrero Bosagna (Corresponding author)](https://orcid.org/0000-0003-1935-5875)

## Funding

This project was funded by the European Union (Grant Agreement No 101000236) and John Templeton Foundation (Grant ID 62167).

## Acknowledgements

We acknowledge support from the National Academic Infrastructure for Supercomputing in Sweden (NAISS), the John Templeton Foundation, the European Union, Svenska Forskningsrådet, and the help from all the attendants from the workshop on GBS-MeDIP data analysis in relation to the GEroNIMO project, especially Sonia Eynard and Gwendal Restoux.
