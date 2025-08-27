# Essential packages
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManaer")

BiocManager::install(c("RDAVIDWebService", "BiocParallel", "AnnotationHub", "clusterProfiler", "pathview"))

install.packages("rJava")

library(rJava)
library(RDAVIDWebService)
library(BiocParallel)
library(AnnotationHub)
library(clusterProfiler)
library(pathview)

# Macaca fascicularis 유전체 주석 데이터베이스 로드
ah <- AnnotationHub()

query_results <- query(ah, pattern = c("Macaca fascicularis", "OrgDb"))
print(query_results)

mfa_db_id <- names(query_results)[1]
org.Mfa.eg.db <- ah[[mfa_db_id]]

print(keytypes(org.Mfa.eg.db))
## file import from my document
file_path <- "D:/xeno/kidney/pathview/FDR_Final_table_5F1.csv"
gene_5F1 <- read.csv(file_path,header=T)


## make the input files
FDR<-gene_5F1[,7]
names(FDR)<-gene_5F1[,1]
allgene <- names(FDR)[FDR >= 0]
siggene <- names(FDR)[FDR <= 0.05]

## Biological Id TranslatoR (bitr)
gene.df <- bitr(allgene, fromType = "SYMBOL", 
                toType = "ENTREZID",
                OrgDb = org.Mfa.eg.db) #★★★★★★★

colnames(mygene1)[1] <- "SYMBOL"
temp <- merge(gene.df, gene_5F1, by="SYMBOL", all=FALSE)
vector1 <- temp[,3]
names(vector1) <- temp[,2]


## Pathview for significant KEGGs

# Renal cell carcinoma
mcf05211 <- pathview(gene.data  = vector1,
                     pathway.id = "05211",
                     species    = "mcf",
                     limit      = list(2, cpd= 1),
                     bins       = list(gene= 18, cpd= 1),
                     low        = list(gene="#00498c", cpd="#00498c"),
                     mid        = list(gene="white", cpd="white"),
                     high       = list(gene="#990000", cpd="#990000"),
                     node.sum   = "max.abs")
