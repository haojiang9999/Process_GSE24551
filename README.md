# Process_GSE24551
For Xiang
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
