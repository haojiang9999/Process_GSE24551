# Process_GSE24551
For Xiang
#####R和Rstudio的安装
我用的[R](https://www.r-project.org/)版本是
```
R version 3.3.3 (2017-03-06) -- "Another Canoe"
```
我觉得3.4应该也没问题，太新的怕有packages不兼容
[Rstudio](https://www.rstudio.com/products/rstudio/#Desktop)下个最新的就行
####一些我用过的参考网址
http://www.bio-info-trainee.com/1399.html

http://www.bio-info-trainee.com/2087.html
https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r/geo/
http://www.bio-info-trainee.com/1085.html

https://www.r-bloggers.com/creating-annotated-data-frames-from-geo-with-the-geoquery-package/

#####首先通过[bioconductor](http://www.bioconductor.org/install/)安装需要的安装包（Installing）
```
source("https://bioconductor.org/biocLite.R")
biocLite()
```
下面是分析需要的packages

Differential expression analysis with limma
```
library(Biobase)
library(GEOquery)
library(limma)
```
#####加载需要的数据，我是把数据下载下来***.soft.gz再加载（Download GSE***.soft.gz file, put it in the current directory, and load it）
```
gse24551 <- getGEO(filename = 'GSE24551_family.soft.gz')
```
#####以后的分析都是基于[GEOquery package](https://bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html) 可以好好参考下
查看数据集GSE中的具体信息
```
head(Meta(gse))
# names of all the GSM objects contained in the GSE
names(GSMList(gse))
# and get the first GSM object on the list
GSMList(gse)[[1]]
# and the names of the GPLs represented
names(GPLList(gse))
```
Converting GSE to an ExpressionSet
由于这个GSE里面有两个GPL所以不能直接转换（好像是这样）
先看下这两个GPL
```
gsmplatforms<-lapply(GSMList(gse24551),function(x){Meta(x)$platform_id})
head(gsmplatforms)
```
有两个 "GPL5175"和'GPL11028'
然后我们先选择'GPL11028'
```
gsmlist=Filter(function(gsm){Meta(gsm)$platform_id=='GPL11028'},GSMList(gse24551))
length(gsmlist)
```
个数应该有160个病人的样本
紧接着我就想用GPL11028这个版本来注释基因提取表达等等（这里有点乱）
首先下载该GPL查看下基本信息
```
#Get the annotation GPL id (see Annotation: GPL10558)
gpl_name <- "GPL11028"
gpl <- getGEO(gpl_name, destdir=".")
Meta(gpl)$title
# Inspect the table of the gpl annotation object
colnames(Table(gpl))
# Get the gene symbol and entrez ids to be used for annotations
Table(gpl)[1:10, c(1, 6, 9, 12)]
dim(Table(gpl))
Table(gpl)$gene_symbol
```
提炼出有基因名字的探针号编号
```
# Get the gene expression data for all the probes with a gene symbol
geneProbes <- which(!is.na(Table(gpl)$gene_symbol))
probeids <- as.character(Table(gpl)$ID[geneProbes])
```
下面是提取"GPL11028"中的160样品中探针表达值，并且按照probeids中的顺序排列好（好像是这样）
```
# make the data matrix from the VALUE columns from each GSM
# being careful to match the order of the probesets in the platform
# with those in the GSMs
data.matrix <- do.call('cbind',lapply(gsmlist,function(x)
                                      {tab <- Table(x)
                                      mymatch <- match(probeids,tab$ID_REF)
                                      return(tab$VALUE[mymatch])
                                    }))
```
这这时候得到的探针数值都是未经转换的，不过我觉得做生存曲线只要比较大小，转不转换大小比较不变
下面就是把上面的matrix转变成[ExpressionSet](http://www.bio-info-trainee.com/1510.html)格式，这个格式可以看讲解
```
require(Biobase)
# go through the necessary steps to make a compliant ExpressionSet
rownames(data.matrix) <- probeids
colnames(data.matrix) <- names(gsmlist)
pdata <- data.frame(samples=names(gsmlist))
rownames(pdata) <- names(gsmlist)
pheno <- as(pdata,"AnnotatedDataFrame")
eset2 <-new('ExpressionSet',exprs=data.matrix,phenoData=pheno)
eset2
```
下面就是构建我们注释好的终极表达矩阵
```
probes <- intersect(probeids, rownames(exprs(eset2)))
length(probes)
geneMatrix <- exprs(eset2)[probes, ]
inds <- which(Table(gpl)$ID %in% probes)
# Check you get the same probes
head(probes)
head(as.character(Table(gpl)$ID[inds]))
# Create the expression matrix with gene ids
geneMatTable <- cbind(geneMatrix, Table(gpl)[inds, c(1, 3, 9, 12)])
head(geneMatTable)
```
这个大矩阵有点太大了电脑收不了，所以提取你想要的基因看看
```
########gene name extract expression info######
gene.name <- "在这里输入你想看的基因名字"
gene.exp<-geneMatTable[grepl(gene.name, geneMatTable$gene_symbol),]
write.table(gene.exp, file = paste("GSE24551_GPL11028_",gene.name,"_DataMatrix.txt"), append = FALSE, quote = F, sep = '\t',
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
```
然后生成了一个txt文件可以用EXCEL打开就成
上面的步骤得到了样品基因的表达值可以通过表达值进行排序了，不过他这个芯片检测的是基因的外显子，所以一个基因有多个探针表示，我就简单粗暴的把同一基因的探针值加和当做基因表达值排序了
现在还少的数据，就是病人生存期数据，下面就是搞这个
```
#loading packeges
library(stringr)
pheno.matrix <- do.call('rbind',lapply(GSMList(gse24551),function(x)
{info <- Meta(x)$characteristics_ch1[4:5]
#####Extract numbers from a character string########
info<-as.numeric(str_extract(info, "\\-*\\d+\\.*\\d*"))
return(info)
}))
#add column names and ouput Table
colnames(pheno.matrix)<-c("survival_years","survival_event")
write.table(pheno.matrix, file = "GSE24551_GPL11028_mataData.txt", append = FALSE, quote = F, sep = '\t',
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
```
Excel打开就是"survival_years"和"survival_event"
结合上面表达值排序在EXCEL里就可以做了，然后放到prism6出图
