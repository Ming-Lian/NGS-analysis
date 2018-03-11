<a name="content">目录</a>

[DESeq2使用指南](#title)
- [导入数据](#input-data)
	- [用tximport导入定量数据 (salmon)](#importing-quantification-data-using-tximport-package)
	- [Count matrix input](#count-matrix)





<h1 name="title">DESeq2使用指南</h1>

<a name="input-data"><h3>导入数据 [<sup>目录</sup>](#content)</h3></a>

<a name="importing-quantification-data-using-tximport-package"><h4>用tximport导入定量数据 (salmon) [<sup>目录</sup>](#content)</h4></a>

将由`salmon quant`产生的`quant.sf	`文件导入，然后构建gene-level DESeqDataSet object

先根据样本信息构建以下格式的数据框变量samples

```
# samples

##           pop center       run condition
## ERR188297 TSI  UNIGE ERR188297         A
## ERR188088 TSI  UNIGE ERR188088         A
## ERR188329 TSI  UNIGE ERR188329         A
## ERR188288 TSI  UNIGE ERR188288         B
## ERR188021 TSI  UNIGE ERR188021         B
## ERR188356 TSI  UNIGE ERR188356         B
```

run和condition两列是必须的

然后根据samples$run来获取对应的数据集的存储路径，保存于变量files中，并且读入 transcripts-genes对应关系的table

```
files <- file.path(dir,"salmon", samples$run, "quant.sf")
names(files) <- samples$run
tx2gene <- read.csv(file.path(dir, "tx2gene.csv"))
```

导入定量数据

```
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
```

最后，根据定量数据txi和样本信息samples构建DESeqDataSet

```
library("DESeq2")
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)
```

<a name="count-matrix"><h4>Count matrix input [<sup>目录</sup>](#content)</h4></a>

简单来说就是构建一个表达谱矩阵，然后用**DESeqDataSetFromMatrix**函数构造DESeqDataSet

详细内容请点 [这里](https://github.com/Ming-Lian/NGS-analysis/blob/master/RNA-seq.md#deseq2)




参考资料:

(1) [Bioconductor tutorial: Analyzing RNA-seq data with DESeq2](http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
