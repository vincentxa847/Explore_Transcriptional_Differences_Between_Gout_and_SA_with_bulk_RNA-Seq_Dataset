# Import the data
de_gout_vs_HC = read.table("C:/Users/user/Desktop/Msc courses/Statistics for Bioinformatics/Data for Report Assessment-20221015/DE_GOUT_vs_HC.csv", header=TRUE, row.names=1, sep="\t")  
de_sa_vs_HC = read.table("C:/Users/user/Desktop/Msc courses/Statistics for Bioinformatics/Data for Report Assessment-20221015/DE_SA_vs_HC.csv", header=TRUE, row.names=1, sep="\t") 
expression_table = read.table("C:/Users/user/Desktop/Msc courses/Statistics for Bioinformatics/Data for Report Assessment-20221015/Expression_Table.csv", header=TRUE, row.names=1, sep="\t")  
sample_information =read.table("C:/Users/user/Desktop/Msc courses/Statistics for Bioinformatics/Data for Report Assessment-20221015/Sample_Information.csv", header = TRUE, row.names = 1, sep = "\t")
annotations = read.table("C:/Users/user/Desktop/Msc courses/Statistics for Bioinformatics/Data for Report Assessment-20221015/Annotations.csv", header =  TRUE, row.names = 1, sep = "\t")

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




