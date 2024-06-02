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



