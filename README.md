# Using RNA-Seq Data to Explore the Transcriptional Differences Between Gout (GA) and Septic Arthritis (SA) in Blood

## Introduction
Patient with SA shares a similar clinical presentation with gout, which makes the diagnosis of SA challenging.
Delay and inadequate treatment can result in irreversible joint destruction of joint and even case-fatality.
As a result, it is crucial to find out the biomarker of these two diseases and can be used in blood diagnosis to distinguish SA from gout.\
In this project, the transcriptional level of 29814 genes in blood samples collected from Healthy group, Gout group and SA group were measured.
*The main aim of this report is to search for the genes that only have significant different between SA and healthy group also between gout and healthy group (target genes).*\
The influence of sex and neutrophil score on genes are also evaluated to clarify the relationship between expression level of genes and clinical measurements.

## Method
In this report, *t-test* is performed using R to compare the expression level of gene in different groups. Each group has 8 samples and the expression value of gene in each sample within a group will be collected and make a comparison with other group to get the p-value and log2fold change. *Pearson correlation test* is used to test the correlation. *ANOVA test* is performed to test the correlation of sex (factor) and gene expression. General Linear model is used to test the correlation of neutrophil score and gene expression. 
