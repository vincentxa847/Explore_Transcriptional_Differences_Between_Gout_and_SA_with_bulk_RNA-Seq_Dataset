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

## Result
Using the p.adj values in de_gout_vs_HC and de_sa_vs_HC to get the significant differential genes between Healthy (HC) and disease groups.
69 genes (de_gout_vs_HC_sig) show significant different between HC and gout (p.adj <0.05) and 13046 genes (de_sa_vs_HC_sig) show significant different between HC and SA (p.adj <0.05).
```
gene_data_sig_gout_vs_HC = subset(de_gout_vs_HC, p.adj < 0.05) 
gene_data_sig_sa_vs_HC = subset(de_sa_vs_HC, p.adj < 0.05)
ggplot(de_gout_vs_HC ,aes(x = log2Fold, y = -log(p.adj,10)))+geom_point(colour = "black")+geom_point(data = gene_data_sig_gout_vs_HC, colour = "red")
ggplot(de_sa_vs_HC ,aes(x = log2Fold, y = -log(p.adj,10)))+geom_point(colour = "black")+geom_point(data = gene_data_sig_sa_vs_HC, colour = "red")
```
![Figure1](https://github.com/vincentxa847/Statistics-for-Bioinformatics-Msc_course/assets/118545004/b258ee61-7506-4551-8366-d1c5ffcb4c00)
![Figure 2](https://github.com/vincentxa847/Statistics-for-Bioinformatics-Msc_course/assets/118545004/419321df-51e9-42ea-ae97-a9c73d102445)



