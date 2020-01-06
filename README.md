---
title: "Readme.Rmd"
author: "Hamid.Bolouri"
date: "12/12/2019"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# This reporsitory contains R scripts used to generate the results in "Integrative network modeling reveals mechanisms underlying T cell exhaustion", Scientific Reports, 2020.

## Description of files:

(Note: To maintain 1:1 mapping to our local files, the original file names are maintained here.)

## Make Cytoscape annotation files from infection expression data
mkCyEdgeAnnotTbl_July2017.R
mkCyNodeAnnot_July2017.R

## Make Cytoscape anotation files from cancer expression
mkCancerEdgeAnnotTbl_July2017.R
mkCancerNodeAnnotTbl_July2017.R

## Additional analysis of data from Andrea Scjietinger's lab
Schietinger_final.R
mkScheitingerFC.tbl.R
edgeScoreSchietinger.R

## Calculate concordance of TCE network edges to published data
concordanceMatching3_25May2018.R

## Example logic simulations
TCE2.0_sim.R # TCE net
simRandTex.R # randomized net
compSimNet2CyNet.R # Comparison

## Soft clustering of genomwide expression data
mfuzz_clustering.R

## Metabolic analysis
metabolicActivity.R

## Prediction of the effect of perturbing each gene in the TCE network
twoStateNet_sensitivityAnalysis3.R


