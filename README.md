# Using RNA-Seq Dataset to Explore the Transcriptional Differences Between Gout (GA) and Septic Arthritis (SA) in Blood

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
**69 genes (de_gout_vs_HC_sig)** show significant different between HC and gout (p.adj <0.05) and **13046 genes (de_sa_vs_HC_sig)** show significant different between HC and SA (p.adj <0.05).
```
gene_data_sig_gout_vs_HC = subset(de_gout_vs_HC, p.adj < 0.05) 
gene_data_sig_sa_vs_HC = subset(de_sa_vs_HC, p.adj < 0.05)
```

Another table is being created to store the t-test results and log2 fold change of gene expression for a different comparison within the DE group. For example, information between gout and HC group in de_sa_vs_HC_sig. The genes show significant up- and down- regulation in one disease group may also have significant change in another disease group. For example, the GATD3A in de_gout_vs_HC_sig group, which has the most dramatic change between HC and gout group (log2fold = 4.757419), also shows a notable increase between HC and SA group (log2fold = 4.323593). This makes GATD3A inappropriate as a biomarker to distinguish SA and gout. It highlights the importance of considering the gene expression between HC and both disease groups when selecting target genes. 

```
# Find out the genes that only significant express in sa but not in gout group# 8515 genes (Same way for genes that only significant express in sa but not in gout group)
# Before this step, merge the DE table with the expression table, and then use a subset function to filter the differentially expressed genes.
comparsion_gout_vs_HC_sig_in_sa_sig = as.data.frame(matrix(0,ncol = 2, nrow =nrow(expression_de_sa_vs_HC_sig))) 
names(comparsion_gout_vs_HC_sig_in_sa_sig) = c("log2fold","p")
row.names(comparsion_gout_vs_HC_sig_in_sa_sig) = row.names(expression_de_sa_vs_HC_sig)
for (row in 1:nrow(expression_de_sa_vs_HC_sig))
{
expression_sa_sig_HC = as.numeric(expression_de_sa_vs_HC_sig[row,c(1:9)]) # HC group 
expression_sa_sig_gout = as.numeric(expression_de_sa_vs_HC_sig[row,c(10:18)]) # Gout group
mean_expression_sa_sig_HC_row = mean(expression_sa_sig_HC)
mean_expression_sa_sig_gout_row = mean(expression_sa_sig_gout)
log2Fold = log2(mean_expression_sa_sig_gout_row)-log2(mean_expression_sa_sig_HC_row)
p_row = t.test(expression_sa_sig_HC,expression_sa_sig_gout)
p_row = p_row$p.value

comparsion_gout_vs_HC_sig_in_sa_sig[row,"log2fold"] = log2Fold
comparsion_gout_vs_HC_sig_in_sa_sig[row,"p"] = p_row
}

comparsion_gout_vs_HC_sig_in_sa_sig_p0.05 = subset(comparsion_gout_vs_HC_sig_in_sa_sig,p>0.05)
de_sa_vs_HC_only_sig_in_sa = de_sa_vs_HC[row.names(comparsion_gout_vs_HC_sig_in_sa_sig_p0.05),]
```
It turns out that **8515 genes have significant differential expression only between HC and SA groups (de_sa_vs_HC_only_sig_in_sa, target genes)** and **3 genes that have significant differential expression only between HC and gout groups (de_gout_vs_HC_only_sig_in_gout, target genes)**. For the 8515 genes, target genes were identified by narrowing down based on log2 fold change.

```
# Sort the table by log2Fold (positive and negative value for up- an down regulation) and select  the top 10 genes as target genes
de_sa_vs_HC_only_sig_in_sa_up = de_sa_vs_HC_only_sig_in_sa[order(de_sa_vs_HC_only_sig_in_sa$log2Fold, decreasing = TRUE),]
de_sa_vs_HC_only_sig_in_sa_up = de_sa_vs_HC_only_sig_in_sa_up[c(1:10),]
de_sa_vs_HC_only_sig_in_sa_down = de_sa_vs_HC_only_sig_in_sa[order(de_sa_vs_HC_only_sig_in_sa$log2Fold),]
de_sa_vs_HC_only_sig_in_sa_down = de_sa_vs_HC_only_sig_in_sa_down[c(1:10),] 
```

