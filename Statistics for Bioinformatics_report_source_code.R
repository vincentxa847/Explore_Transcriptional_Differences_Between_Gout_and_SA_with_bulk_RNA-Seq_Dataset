#### Import the data ####
de_gout_vs_HC = read.table("../Data for Report Assessment-20221015/DE_GOUT_vs_HC.csv", header=TRUE, row.names=1, sep="\t")  
de_sa_vs_HC = read.table("../Data for Report Assessment-20221015/DE_SA_vs_HC.csv", header=TRUE, row.names=1, sep="\t") 
expression_table = read.table("../Data for Report Assessment-20221015/Expression_Table.csv", header=TRUE, row.names=1, sep="\t")  
sample_information =read.table("../Data for Report Assessment-20221015/Sample_Information.csv", header = TRUE, row.names = 1, sep = "\t")
annotations = read.table("../Data for Report Assessment-20221015/Annotations.csv", header =  TRUE, row.names = 1, sep = "\t")

#### Merge expression,annotation and de table. Then create table to subset the de genes of sa and gouts respectively ####
# merge expression and annotation table 
expression_annotation = merge(annotations, expression_table, by.x = 0,by.y = 0)
row.names(expression_annotation) = expression_annotation[,1]
names(expression_annotation)[1] = "GENE_ID"

# Get significant gene of gout,sa and HC
de_gout_vs_HC_sig = subset(de_gout_vs_HC, p.adj < 0.01) 
de_sa_vs_HC_sig = subset(de_sa_vs_HC, p.adj <0.01)
nrow(de_gout_vs_HC_sig) # 12 (p.adj < 0.01)
nrow(de_sa_vs_HC_sig) # 9894 (p.adj < 0.01)

# Create expression table of only significant gene of de_gout here # 69 (p.adj < 0.05) 
# (This table is very useful) 
expression_de_gout_merge = merge(expression_table, de_gout_vs_HC, by.x=0,by.y = 0)
expression_de_gout_anno_merge = merge(annotations,expression_de_gout_merge,by.x=0,by.y = 1)
# expression_de_gout_anno_merge use to index gene symbols and gene id (not cut the gene_id)
expression_de_gout_anno_sig = subset(expression_de_gout_anno_merge, p.adj<0.05)
gene_symbols = expression_de_gout_anno_sig[,"symbol"]
row.names(expression_de_gout_anno_sig) = gene_symbols
expression_de_gout_anno_sig = expression_de_gout_anno_sig[7:36]

# Create expression table of only significant gene de_sa here # 12999 (p.adj < 0.05)
expression_de_sa_merge = merge(expression_table, de_sa_vs_HC, by.x=0,by.y = 0)
expression_de_sa_anno_merge = merge(annotations,expression_de_sa_merge,by.x=0,by.y = 1) 
expression_de_sa_anno_sig = subset(expression_de_sa_anno_merge, p<0.05)
gene_symbols = expression_de_sa_anno_sig[,"symbol"]
row.names(expression_de_sa_anno_sig) = gene_symbols
expression_de_sa_anno_sig = expression_de_sa_anno_sig[7:36]

#### Cherry-pick some genes for discussion ####
# Create the histogram of the sample (by column) 
library(ggplot2)
ggplot(expression_table, aes(x=log10(HC_1)))+geom_histogram()

# Create a histogram of a gene,take IFT46 and ASPH gene for example (by row) and make a plot 
IFT46 = data.frame(t(expression_de_gout_anno_sig["IFT46",-c(28:30)]))
ggp = ggplot(IFT46,aes(x =IFT46)) +geom_density()
ggp # This shows the gene_expression value across the x-axis and their frequencies in y-axis. In this plot, most of the expression values of this gene are between 400-450 
ASPH = data.frame(t(expression_de_gout_anno_sig["ASPH",-c(28:30)]))
ggp = ggplot(ASPH,aes(x =ASPH)) +geom_density()
ggp # we can see that most of the expression have some degree of centrality

# Make correlation plot 
genes = c("SULT4A1","SNORA73B") # genes that significant higher in gout 
SULT4A1_SNORA73B = data.frame(t(expression_de_gout_anno_sig[genes,-c(28:30)])) # the way to select multiple gene in matrix
SULT4A1_SNORA73B = SULT4A1_SNORA73B[c(10:18),] # select the gout group here
ggp = ggplot(SULT4A1_SNORA73B,aes(x=SULT4A1,y=SNORA73B))+ geom_point()
ggp # the genes here have no clear correlation, even compare within gout group, maybe we don't have enough group
# Make a scatter plot, comparing patient with different situations (neutrophils score)
HIPK2 = data.frame(t(expression_de_gout_anno_sig["HIPK2",-c(28:30)]))
HIPK2$neutrophils= sample_information$NEUTROPHILS # select the target gene and add a sample information to see if neutrophils will affect the gene expression
ggp = ggplot(HIPK2,aes(x=HIPK2,y=neutrophils))+geom_point()
ggp # the HIPK2 gene and neutrophils score show correlation  (report in assignment)

# Make density plot, comparing patient with different situations (sex)
# KLHDC7A
KLHDC7A_spilt_by_sex = data.frame(t(expression_de_gout_anno_sig["KLHDC7A",-c(28:30)]))
KLHDC7A_spilt_by_sex$sex = sample_information$SEX
KLHDC7A_spilt_by_sex = KLHDC7A_spilt_by_sex[c(10:18),] # select the gout group here
ggp = ggplot(KLHDC7A_spilt_by_sex,aes(x=KLHDC7A,fill=sex))+geom_density()
ggp # KLHDC7A seems to affect by sex both in gout group and in whole group (report in assignment)
# Try another gene (LINC00327 here)
LINC00327_spilt_by_sex = data.frame(t(expression_de_gout_anno_sig["LINC00327",-c(28:30)]))
LINC00327_spilt_by_sex$sex = sample_information$SEX
LINC00327_spilt_by_sex =LINC00327_spilt_by_sex[c(10:18),] # select the gout group here
ggp = ggplot(LINC00327_spilt_by_sex,aes(x=LINC00327,fill=sex))+geom_density() # In gout group, female have two peak with a single peak at 60
ggp
# Make a boxplot spilt by sex
ggp = ggplot(KLHDC7A_spilt_by_sex, aes(x=sex,y=KLHDC7A))+ geom_boxplot()
ggp
ggp = ggplot(LINC00327_spilt_by_sex, aes(x=sex,y=LINC00327))+ geom_boxplot()
ggp

#### Using expression table to find genes that only show significantly differential expression in only one group ####
# Start to make the significant gene expression table in gout vs HC and sp vs HC 
de_gout_vs_HC_sig_id = row.names(de_gout_vs_HC_sig) # p<0.01
expression_de_gout_vs_HC_sig = expression_table[de_gout_vs_HC_sig_id,]
de_sa_vs_HC_sig_id = row.names(de_sa_vs_HC_sig)
expression_de_sa_vs_HC_sig = expression_table[de_sa_vs_HC_sig_id,]
# Compare expression between HC and sa base on de_sa_vs_HC (we can use this to compare within group to answer task1)
# take the second gene "ENSG00000115919" for example here
expression_sa_sig_HC_gene2 = as.numeric(expression_de_sa_vs_HC_sig[2,c(1:9)])
expression_sa_sig_SA_gene2 = as.numeric(expression_de_sa_vs_HC_sig[2,c(19:27)])
mean_expression_sa_sig_HC = mean(expression_sa_sig_HC_gene2)
mean_expression_sa_sig_SA = mean(expression_sa_sig_SA_gene2)
p = t.test(expression_sa_sig_HC_gene2,expression_sa_sig_SA_gene2)
p = p$p.value
p
# Compare expression between HC and gout base on de_sa_vs_HC (we can use this to compare within group to answer task1)
# take the second gene "ENSG00000115919" for example here
expression_sa_sig_HC_gene2 = as.numeric(expression_de_sa_vs_HC_sig[2,c(1:9)])
expression_sa_sig_gout_gene2 = as.numeric(expression_de_sa_vs_HC_sig[2,c(10:18)])
mean_expression_sa_sig_HC = mean(expression_sa_sig_HC_gene2)
mean_expression_sa_sig_gout = mean(expression_sa_sig_gout_gene2)
p = t.test(expression_sa_sig_HC_gene2,expression_sa_sig_gout_gene2)
p = p$p.value
# Create a table to store the results show the gout group expression level of HC_sa_sig_gene, 
# this can help select the gene that only sig in sa group 
comparsion_gout_vs_HC_sig_in_sa_sig = as.data.frame(matrix(0,ncol = 2, nrow =nrow(expression_de_sa_vs_HC_sig)))
names(comparsion_gout_vs_HC_sig_in_sa_sig) =c("log2fold","p")
row.names(comparsion_gout_vs_HC_sig_in_sa_sig) =row.names(expression_de_sa_vs_HC_sig)
for (row in 1:nrow(expression_de_sa_vs_HC_sig))
{
  expression_sa_sig_HC = as.numeric(expression_de_sa_vs_HC_sig[row,c(1:9)])
  expression_sa_sig_gout = as.numeric(expression_de_sa_vs_HC_sig[row,c(10:18)])
  mean_expression_sa_sig_HC_row = mean(expression_sa_sig_HC)
  mean_expression_sa_sig_gout_row = mean(expression_sa_sig_gout)
  log2Fold = log2(mean_expression_sa_sig_gout_row)-log2(mean_expression_sa_sig_HC_row)
  p_row = t.test(expression_sa_sig_HC,expression_sa_sig_gout)
  p_row = p_row$p.value
  
  comparsion_gout_vs_HC_sig_in_sa_sig[row,"log2fold"] = log2Fold
  comparsion_gout_vs_HC_sig_in_sa_sig[row,"p"] = p_row
}
# find out the gene that only significant express in sa group but not in gout # 8515 (target gene)
comparsion_gout_vs_HC_sig_in_sa_sig_p0.05 = subset(comparsion_gout_vs_HC_sig_in_sa_sig,p>0.05)
# get the value of those 8515 target genes from de_sa_vs_HC
de_sa_vs_HC_only_sig_in_sa = de_sa_vs_HC[row.names(comparsion_gout_vs_HC_sig_in_sa_sig_p0.05),]
# sort the table by p column (lowest will be first) in this case we have many target genes so we can use the most significant gene as target
de_sa_vs_HC_only_sig_in_sa = de_sa_vs_HC_only_sig_in_sa[order(de_sa_vs_HC_only_sig_in_sa$p),]
# Use the most significant gene as target to show the effect of sex on its expression , ENSG00000198074 AKR1B10 (just look at sa group)
target1 = data.frame(t(expression_table["ENSG00000198074",]))
target1$SEX = sample_information$SEX
target1_in_sa = target1[c(19:27),]
ggp = ggplot(target1_in_sa,aes(x = SEX,y=ENSG00000198074))+geom_boxplot()
ggp # sex have no significant effect on ENSG00000198074

# Create another table to store the results show that in HC_gout_sig_gene, the expression level of sa group 
# this can help select the gene that only sig in gout group 
comparsion_sa_vs_HC_sig_in_gout_sig = as.data.frame(matrix(0,ncol = 2, nrow =nrow(expression_de_gout_vs_HC_sig)))
names(comparsion_sa_vs_HC_sig_in_gout_sig) =c("log2fold","p")
row.names(comparsion_sa_vs_HC_sig_in_gout_sig) =row.names(expression_de_gout_vs_HC_sig)
for (row in 1:nrow(expression_de_gout_vs_HC_sig))
{
  expression_gout_sig_HC = as.numeric(expression_de_gout_vs_HC_sig[row,c(1:9)])
  expression_gout_sig_sa = as.numeric(expression_de_gout_vs_HC_sig[row,c(19:27)])
  mean_expression_gout_sig_HC_row = mean(expression_gout_sig_HC)
  mean_expression_gout_sig_sa_row = mean(expression_gout_sig_sa)
  log2Fold = log2(mean_expression_gout_sig_sa_row)-log2(mean_expression_gout_sig_HC_row)
  p_row = t.test(expression_gout_sig_HC,expression_gout_sig_sa)
  p_row = p_row$p.value
  
  comparsion_sa_vs_HC_sig_in_gout_sig[row,"log2fold"] = log2Fold
  comparsion_sa_vs_HC_sig_in_gout_sig[row,"p"] = p_row
}
# find out the gene that only significant express in gout group but not in sa # 3 (target gene)
comparsion_sa_vs_HC_sig_in_gout_sig_p0.05 = subset(comparsion_sa_vs_HC_sig_in_gout_sig,p>0.05) 
# get the value of those 3 target genes from de_gout_vs_HC
de_gout_vs_HC_only_sig_in_gout = de_gout_vs_HC[row.names(comparsion_sa_vs_HC_sig_in_gout_sig_p0.05),]
# sort the table by p column (lowest will be first) in this case we have 3 target genes so we can use the most significant gene as target
de_gout_vs_HC_only_sig_in_gout = de_gout_vs_HC_only_sig_in_gout[order(de_gout_vs_HC_only_sig_in_gout$p),]

# Use the most significant gene as target to show the effect of sex on its expression , ENSG00000198363 ASPH (just look at gout group)
target2 = data.frame(t(expression_table["ENSG00000198363",]))
target2$SEX = sample_information$SEX
target2_in_gout = target2[c(10:18),]
ggp = ggplot(target2_in_gout,aes(x = SEX,y=ENSG00000198363))+geom_boxplot()
ggp

#### Target genes data for discussion ####
# Significant gene only in sa group
target1 = data.frame(t(expression_table["ENSG00000198074",]))
target1$NEUTROPHILS = sample_information$NEUTROPHILS
target1$SEX = sample_information$SEX
target1$SEX = as.factor(target1$SEX)
target2 = data.frame(t(expression_table["ENSG00000115919",]))
target2$NEUTROPHILS = sample_information$NEUTROPHILS
target2$SEX = sample_information$SEX
target2$SEX = as.factor(target2$SEX)
target3 = data.frame(t(expression_table["ENSG00000124102",]))
target3$NEUTROPHILS = sample_information$NEUTROPHILS
target3$SEX = sample_information$SEX
target3$SEX = as.factor(target3$SEX)
target4 = data.frame(t(expression_table["ENSG00000188373",]))
target4$NEUTROPHILS = sample_information$NEUTROPHILS
target4$SEX = sample_information$SEX
target4$SEX = as.factor(target4$SEX)
target5 = data.frame(t(expression_table["ENSG00000135114",]))
target5$NEUTROPHILS = sample_information$NEUTROPHILS
target5$SEX = sample_information$SEX
target5$SEX = as.factor(target5$SEX)
# Significant gene only in gout group
target6 = data.frame(t(expression_table["ENSG00000198363",]))
target6$NEUTROPHILS = sample_information$NEUTROPHILS
target6$SEX = sample_information$SEX
target6$SEX = as.factor(target6$SEX)
target7 = data.frame(t(expression_table["ENSG00000118096",]))
target7$NEUTROPHILS = sample_information$NEUTROPHILS
target7$SEX = sample_information$SEX
target7$SEX = as.factor(target7$SEX)
target8 = data.frame(t(expression_table["ENSG00000130540",]))
target8$NEUTROPHILS = sample_information$NEUTROPHILS
target8$SEX = sample_information$SEX
target8$SEX = as.factor(target8$SEX)

#### Using density plot to see the effect of sex on target genes expression and the distribution of target genes in different groups ####
# Because we just have 30 sample groups, so use density plot rather boxplot is better 
ggp = ggplot(target2,aes(x=ENSG00000115919,fill=SEX)) + geom_density()
ggp
# Do some statistic test in target gene (testing)
mean(target8$NEUTROPHILS)
summary(target8[19:27,]) 
sd(target8[19:27,1])
# See the distribution of target gene in different group 
ggp = ggplot(target8[1:9,],aes(x=ENSG00000130540))+geom_density()
ggp # HC, most of them distribute in x = 5
ggp = ggplot(target8[10:18,],aes(x=ENSG00000130540))+geom_density()
ggp # gout, most of them distribute in about x=50, have a significant peak in about x= 550
ggp = ggplot(target8[19:27,],aes(x=ENSG00000130540))+geom_density()
ggp # gout most of them distribute in x = 20

ggp = ggplot(target1[1:9,],aes(x=ENSG00000198074))+geom_density()
ggp # HC, most of them distribute in x = 50, have a significant peak in about x= 75 and 100
ggp = ggplot(target1[10:18,],aes(x=ENSG00000198074))+geom_density()
ggp # gout, most of them distribute in about x=45
ggp = ggplot(target1[19:27,],aes(x=ENSG00000198074))+geom_density()
ggp # gout most of them distribute in x = 4500 , have a broad peak in about x = 6500

#### Correlation test to EVALUATE THE EFFECT OF NEUTROPHILS ####
cor(target8$ENSG00000130540,target8$NEUTROPHILS) # overall
cor(target8[1:9,]$ENSG00000130540,target8[1:9,]$NEUTROPHILS) # HC group
cor(target8[10:18,]$ENSG00000130540,target8[10:18,]$NEUTROPHILS) # gout group
cor(target8[19:27,]$ENSG00000130540,target8[19:27,]$NEUTROPHILS) # sa group
# plot the correlation
ggp = ggplot(target8,aes(x=ENSG00000130540,y=NEUTROPHILS))+geom_point()
ggp # no significant correlation ( for neutrophil score )
ggp = ggplot(target8,aes(y=ENSG00000130540,fill=SEX))+geom_boxplot()
ggp # no significant correlation ( for sex )
ggp = ggplot(target8[10:18,],aes(y=ENSG00000130540,fill=SEX))+geom_boxplot()
ggp # see the correlation within group
ggp = ggplot(target8[1:9,],aes(y=ENSG00000130540,fill=SEX))+geom_boxplot()
ggp # see the correlation within  target 8 HC group (have some different) # figure 1
ggp = ggplot(target6[10:18,],aes(y=ENSG00000198363,fill=SEX))+geom_boxplot()
ggp # see the correlation within  target 6 gout group (have some different) # figure 2
ggp = ggplot(target6,aes(y=ENSG00000198363,fill=SEX))+geom_boxplot() # figure 3
ggp # the correlation in target6 across group show an opposite trend compare to within gout group
 

#### GLM [for sex (factor), anova] ####
## target 6
model_target6_sex = lm(target6$ENSG00000198363 ~ target6$SEX)
anova(model_target6_sex)
# Response: target6$ENSG00000198363
# Df   Sum Sq Mean Sq F value Pr(>F)
# target6$SEX  1    10395   10395  0.0078 0.9303 # no significant relation across group. although in figure 2, it seems to have some relation between expression and sex
# Residuals   25 33323283 1332931
 
## target  8
model_target8_sex = lm(target8$ENSG00000130540 ~ target8$SEX)
anova(model_target8_sex)
# Response: target8$ENSG00000130540
# Df Sum Sq Mean Sq F value Pr(>F)
# target8$SEX  1  10338   10338  0.9959 0.3279 # no significant relation across group
# Residuals   25 259514   10380

# Response: target8[10:18, ]$ENSG00000130540
# Df Sum Sq Mean Sq F value Pr(>F)
# target8[10:18, ]$SEX  1  30450   30450  1.1371 0.3217# no significant relation within gout group. although in figure 1, it seems to have some relation between expression and sex 
# Residuals             7 187456   26780
  

# GLM [for NEUTROPHILS SCORE (variable), linear model]
# HIPK2 (ENSG00000064393, most significant gene in de_gout_HC) see the correlation in geom_point
HIPK2 = data.frame(t(expression_de_gout_anno_sig["HIPK2",-c(28:30)]))
HIPK2$neutrophils= sample_information$NEUTROPHILS
model_HIPK2_neu = lm(HIPK2$HIPK2 ~ HIPK2$neutrophils)

summary(model_HIPK2_neu)
# Residuals:
#  Min      1Q  Median      3Q     Max 
# -853.89 -434.55    4.61  413.67  581.59 

lm(formula = HIPK2$HIPK2 ~ HIPK2$neutrophils)
# Coefficients:
# (Intercept)         529.04     196.86   2.687   0.0126 *  # have correlation
#   HIPK2$neutrophils    99.42      20.85   4.769 6.78e-05 ***

# Visualize the model fit (for factor) TU6
ggp = ggplot(target8, aes(x = SEX, y = ENSG00000130540 , fill = SEX))+ geom_violin()+stat_summary(fun=mean)
ggp # figure4 from violin plot we can see that the samples are not a good distribution here
# this may be the reason why it look to have some different in boxblot but actually no relation by annova test

