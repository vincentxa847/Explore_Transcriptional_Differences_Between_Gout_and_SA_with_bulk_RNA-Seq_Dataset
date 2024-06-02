#### Import the data ####
de_gout_vs_HC = read.table("C:/Users/user/Desktop/Msc courses/Statistics for Bioinformatics/Data for Report Assessment-20221015/DE_GOUT_vs_HC.csv", header=TRUE, row.names=1, sep="\t")  
de_sa_vs_HC = read.table("C:/Users/user/Desktop/Msc courses/Statistics for Bioinformatics/Data for Report Assessment-20221015/DE_SA_vs_HC.csv", header=TRUE, row.names=1, sep="\t") 
expression_table = read.table("C:/Users/user/Desktop/Msc courses/Statistics for Bioinformatics/Data for Report Assessment-20221015/Expression_Table.csv", header=TRUE, row.names=1, sep="\t")  
sample_information =read.table("C:/Users/user/Desktop/Msc courses/Statistics for Bioinformatics/Data for Report Assessment-20221015/Sample_Information.csv", header = TRUE, row.names = 1, sep = "\t")
annotations = read.table("C:/Users/user/Desktop/Msc courses/Statistics for Bioinformatics/Data for Report Assessment-20221015/Annotations.csv", header =  TRUE, row.names = 1, sep = "\t")

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

# TU4 use the most significant gene as target to show the effect of sex on its expression , ENSG00000198363 ASPH (just look at gout group)
target2 = data.frame(t(expression_table["ENSG00000198363",]))
target2$SEX = sample_information$SEX
target2_in_gout = target2[c(10:18),]
ggp = ggplot(target2_in_gout,aes(x = SEX,y=ENSG00000198363))+geom_boxplot()
ggp
