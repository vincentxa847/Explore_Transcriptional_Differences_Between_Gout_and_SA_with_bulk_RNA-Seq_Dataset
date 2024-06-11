# Using RNA-Seq Dataset to Explore the Transcriptional Differences Between Gout (GA) and Septic Arthritis (SA) Blood Sample

## Introduction
Patient with SA shares a similar clinical presentation with gout, which makes the diagnosis of SA challenging.
Delay and inadequate treatment can result in irreversible joint destruction of joint and even case-fatality.
As a result, it is crucial to find out the biomarker of two diseases and can be used in blood diagnosis to distinguish SA from gout.\
In this project, the transcriptional levels of 29814 genes in blood samples collected from Healthy group, Gout group and SA group were measured.
* The main aim of this report is to identify the genes that show significant differences only between the SA and Healthy groups, as well as between the Gout and Healthy groups (target genes).*\
Additionally, the influence of sex and neutrophil score on gene expression is evaluated to clarify the relationship between gene expression levels and clinical measurements.

## Method
In this report, *t-test* is performed using R to compare gene expression levels across different groups. Each group consists of 8 samples and gene expression values from each sample within a group were collected  and compared with those from other groups to obtain p-values and log2 fold changes. *Pearson correlation test* was used to assess the correlation. *ANOVA test* was performed to examine the correlation between sex (factor) and gene expression. General Linear model was used to investigate the correlation between neutrophil score and gene expression. 

## Result
Using the p.adj values in de_gout_vs_HC and de_sa_vs_HC (Table of Differential Values) to get the significant differential genes between Healthy (HC) and disease groups.
**69 genes (de_gout_vs_HC_sig)** show significant different between HC and gout (p.adj <0.05) and **13046 genes (de_sa_vs_HC_sig)** show significant different between HC and SA (p.adj <0.05).
```
gene_data_sig_gout_vs_HC = subset(de_gout_vs_HC, p.adj < 0.05) 
gene_data_sig_sa_vs_HC = subset(de_sa_vs_HC, p.adj < 0.05)
```

Another table is being created to store the t-test results and log2 fold change of gene expression for a different comparison within the DE group. For example, information between gout and HC group in de_sa_vs_HC_sig. The genes show significant up- and down- regulation in one disease group may also have significant change in another disease group. For example, the GATD3A in de_gout_vs_HC_sig group, which has the most dramatic change between HC and gout group (log2fold = 4.757419), also shows a notable increase between HC and SA group (log2fold = 4.323593). This makes GATD3A inappropriate as a biomarker to distinguish SA and gout. It highlights the importance of considering the gene expression between HC and both disease groups when selecting target genes. 

```
# 8515 genes (Same process for genes that only significant express in sa but not in gout group)
# Before this step, merge the DE table with the expression table.
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
It turns out that **8515 genes have significant differential expression exclusively between HC and SA groups (de_sa_vs_HC_only_sig_in_sa, target genes)** and **3 genes that have significant differential expression exclusively between HC and gout groups (de_gout_vs_HC_only_sig_in_gout, target genes)**. For the 8515 genes, target genes were identified by narrowing down based on log2 fold change.

```
# Sort the table by log2Fold (positive and negative value for up- an down regulation) and select the top 10 genes as target genes
de_sa_vs_HC_only_sig_in_sa_up = de_sa_vs_HC_only_sig_in_sa[order(de_sa_vs_HC_only_sig_in_sa$log2Fold, decreasing = TRUE),]
de_sa_vs_HC_only_sig_in_sa_up = de_sa_vs_HC_only_sig_in_sa_up[c(1:10),]
de_sa_vs_HC_only_sig_in_sa_down = de_sa_vs_HC_only_sig_in_sa[order(de_sa_vs_HC_only_sig_in_sa$log2Fold),]
de_sa_vs_HC_only_sig_in_sa_down = de_sa_vs_HC_only_sig_in_sa_down[c(1:10),]
de_sa_vs_HC_only_sig_in_sa_upanddown = rbind(de_sa_vs_HC_only_sig_in_sa_up,de_sa_vs_HC_only_sig_in_sa_down)

ggplot(de_sa_vs_HC ,aes(x = log2Fold, y = -log(p.adj,10)))+ geom_point(colour = "black")+ geom_point(data = de_sa_vs_HC_only_sig_in_sa_upanddown, colour = "red")+ geom_text_repel(data = de_sa_vs_HC_only_sig_in_sa_upanddown,aes(label = symbol), vjust = -1, hjust = 1)

ggplot(de_gout_vs_HC ,aes(x = log2Fold, y = -log(p.adj,10)))+geom_point(colour = "black")+geom_point(data = de_gout_vs_HC_only_sig_in_gout, colour = "red")+geom_text_repel(data = de_gout_vs_HC_only_sig_in_gout,aes(label = symbol), vjust = -1, hjust = 1)

```
![de_sa_vs_HC_only_sig_in_sa_upanddown](https://github.com/vincentxa847/Statistics-for-Bioinformatics-Msc_course/assets/118545004/d1faf300-ce9d-4103-9331-286c19be26dc)\
*Genes showing significant differential expression exclusively between HC and SA groups, top 10 upregulated and downregulated genes were plotted.*\
![de_gout_vs_HC_only_sig_in_gout](https://github.com/vincentxa847/Statistics-for-Bioinformatics-Msc_course/assets/118545004/e95adfc4-82b3-4a29-880e-ad49b65e34ee)\
*Genes showing significant differential expression exclusively between between HC and gout.*

After selecting the target genes, the effect of clinical measurements is evaluated. From the de_gout_vs_HC_sig dataset, *HIPK2* was selected and plotted. The data clearly show that HIPK2 expression is positively correlated with neutrophil scores across the three groups. Pearson correlation test yields the correlation value of 0.6902029, indicating a positive correlation between HIPK2 expressions and neutrophil score. To further confirm this relationship, General Linear Model (GLM) was applied. The GLM variables for HIPK2 and neutrophil scores indicate an intercept of 529.04 and a slope of 99.42.

```
HIPK2 = data.frame(t(expression_de_gout_anno_sig["HIPK2",-c(28:30)]))
HIPK2$neutrophils= sample_information$NEUTROPHILS 
cor(x=HIPK2$HIPK2,y=HIPK2$neutrophils) # 0.6902029
model_HIPK2_NEU = lm(HIPK2$HIPK2~HIPK2$neutrophils)
summary(model_HIPK2_NEU)

ggplot(HIPK2, aes(x = HIPK2,y=neutrophils)) + geom_point()+ geom_smooth(method = "lm", se = FALSE)
```
![HIPK2_NEUTROPHILS](https://github.com/vincentxa847/Statistics-for-Bioinformatics-Msc_course/assets/118545004/f4880b4d-fce1-451d-959f-8b6f78ae020f)\
*Correlation between HIPK2 expression and neutrophil score*

However, there is no clear relationship between HIPK2 expression and sex. The p-value of the F-statistic from the ANOVA test is 0.858, indicating that the null hypothesis cannot be rejected.

```
HIPK2$sex= sample_information$SEX
model_HIPK2_sex = lm(HIPK2$HIPK2 ~ HIPK2$sex)
anova(model_HIPK2_sex)

ggplot(HIPK2,aes(y=HIPK2,x=sex))+ geom_boxplot()
```
![HIPK2_SEX](https://github.com/vincentxa847/Statistics-for-Bioinformatics-Msc_course/assets/118545004/b3af3b42-c1d7-418a-8708-b05c3073b61b)\
*Distribution of HIPK2 in male and female*

## Discussion
In this report, blood markers that can be used to distinguish gout and SA were identified. Three target genes in the gout group were found: *ASPH, IFT46, and SULT4A*. [Aspartate β-hydroxylase (ASPH) is a highly conserved enzyme that is widely expressed in proliferating placenta trophoblastic cells](https://pubmed.ncbi.nlm.nih.gov/33125119/). This enzyme is also a biomarker and potential drug target of cancer cells. [IFT46 is a core component of the intraflagellar transport machinery and is related to the formation of all cilia](https://pubmed.ncbi.nlm.nih.gov/25722189/). In the provided data, IFT46 shows significant downregulation in gout group (p.adj =0.006439270). Most research on IFT46 focuses on its function in cilia and therefore its role in human disease remains future investigation. SULT4A1 was identified in both human and rat brain and belongs to [sulphotransferase (SULT) gene family](https://pubmed.ncbi.nlm.nih.gov/10698717/). [It is associated with several neurologic disorders such as Phelan McDermid syndrome](https://pubmed.ncbi.nlm.nih.gov/18823757/). The data provided show a notable upregulation of SULT4A1 in the gout group, showing its role in inflammatory disease worth investigating.\
Out of the 8515 target genes in the gout group, the gene that exhibits the most significant difference is *AKR1B10*, which shows a notable upregulation. AKR1B10 is a human nicotinamide adenine dinucleotide phosphate (NADPH)-dependent reductase belonging to the aldo-keto reductase (AKR) 1B subfamily. [The AKR1B family is associated with diabetes and Inflammatory Disease](https://pubmed.ncbi.nlm.nih.gov/30362099/) but its mechanism in disease remains unclear and requires further research to confirm. While target genes in this dataset exhibit significant differences between control and disease groups, they also encompass a wide range of biological functions. The diversity makes it challenging to investigate their molecular mechanism in disease. 
