<p name="content">目录</p>

[Analysis pipeline for RNA-seq](#title)
- [测序数据下载](#download)
- [比对与定量](#map-quant)
	- [Salmon流程](#salmon)
		- [创建索引](#salmon-index)
		- [定量](#salmon-quant)
	- [subread流程](#subread)
		- [创建索引](#subread-index)
		- [比对](#subread-map)
		- [定量](#subread-quant)
	- [hisat2-stringtie流程](#hisat2-stringtie)
		- [hisat2创建索引](#hisat2-index)
		- [hisat2比对](#hisat2-map)
		- [stringtie转录本拼接](#stringtie-assm)
		- [stringtie定量](#stringtie-quant)
- [差异表达分析](#diff-exp)
	- [DESeq2](#deseq2)
	- [Ballgown](#ballgown)

<h1 name="title">Analysis pipeline for RNA-seq</h1>

<p align="center"> <img src=/picture/hisat2-stringtie-ballgown.png align="middle"> </p>

<a name="download"><h3 >测序数据下载 [<sup>目录</sup>](#content)</h3></a>

---

参见： https://github.com/Ming-Lian/Memo/blob/master/ChIP-seq-pipeline.md#get-data

<h3 name="map-quant">比对与定量</h3>

---

<a name="salmon"><h3 ><li>Salmon流程 [<sup>目录</sup>](#content)</li></h3></a>

<br size=1 />

不需要比对，直接对转录水平进行定量

<h4 name="salmon-index">创建索引</h4>

```
$ salmon index -t Arabidopsis_thaliana.TAIR10.28.cdna.all.fa.gz -i athal_index_salmon
```

<h4 name="salmon-quant">定量</h4>

```
#! /bin/bash
index=~/rna_seq_practice_2/data/ref/athal_index_salmon ## 指定索引文件夹
quant=~/rna_seq_practice_2/work/quant_salmon # 指定定量文件输出文件夹
for fn in ERR1698{194..209};
do
	sample=`basename ${fn}` # basename命令用于去掉路径信息，返回纯粹的文件名，如果指定的文件的扩展名，则将扩展名也一并去掉。
	echo "Processin sample ${sampe}"
	salmon quant -i $index -l A \
		-1 ${sample}_1.fastq.gz \
		-2 ${sample}_2.fastq.gz \
		-p 5 -o $quant/${sample}_quant
done
# quant_salmon.sh
```

<a name="subread"><h3><li>subread流程 [<sup>目录</sup>](#content)</li></h3></a>

<h4 name="subread-index">创建索引</h4>

```
$ gunzip Arabidopsis_thaliana.TAIR10.28.dna.genome.fa.gz
$ subread-buildindex -o athal_index_subread   Arabidopsis_thaliana.TAIR10.28.dna.genome.fa
```

<h4 name="subread-map">比对</h4>

```
#! /bin/bash
index=~/rna_seq_practice_2/data/ref/athal_index_subread
map=~/rna_seq_practice_2/work/map
for fn in ERR1698{194..209};
do
	sample=`basename ${fn}`
	echo "Processin sample ${sampe}" 
	subjunc -i $index \
		-r ${sample}_1.fastq.gz \
		-R ${sample}_2.fastq.gz \
		-T 5 -o $map/${sample}_subjunc.bam
# 比对的sam自动转为bam，但是并不按照参考基因组坐标排序
done
# map_subjunc.sh
```

<h4 name="subread-quant">定量</h4>

```
featureCounts=~/anaconda2/bin
# gff3=~/rna_seq_practice_2/data/ref/Arabidopsis_thaliana.TAIR10.28.gff3.gz
gtf=~/rna_seq_practice_2/data/ref/Arabidopsis_thaliana.TAIR10.28.gtf
count=~/rna_seq_practice_2/work/quant_subjunc
nohup $featureCounts/featureCounts  -T 5 -p -t exon -g gene_name -a $gtf -o  $count/counts.txt   *.bam &
nohup $featureCounts/featureCounts  -T 5 -p -t exon -g gene_id -a $gtf -o  $count/counts_id.txt   *.bam &
```

<a name="hisat2-stringtie"><h3><li>hisat2-stringtie流程 [<sup>目录</sup>](#content)</li></h3></a>

<h4 name="hisat2-index">hisat2创建索引</h4>

```
# build reference index
##  using the python scripts included in the HISAT2 package, extract splice-site and exon information from the gene
annotation fle
$ extract_splice_sites.py chrX_data/genes/chrX.gtf >chrX.ss
$ extract_exons.py chrX_data/genes/chrX.gtf >chrX.exon
##  build a HISAT2 index
$ hisat2-build --ss chrX.ss --exon chrX.exon chrX_data/genome/chrX.fa chrX_tran
```

<h4 name="hisat2-map">hisat2比对</h4>

```
$ hisat2 -p 10 --dta -x chrX_tran -1 reads1_1.fastq -2 reads1_2.fastq | samtools sort -@ 8 -O bam -o reads1.sort.bam 1>map.log 2>&1
```
`Usage: hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]`
> - -p Number of threads to use
> - --dta reports alignments tailored for transcript assemblers
> - -x Hisat2 index
> - -1 The 1st input fastq file of paired-end reads
> - -2 The 2nd input fastq file of paired-end reads
> - -S File for SAM output (default: stdout)


<h4 name="stringtie-assm">stringtie转录本拼接</h4>

```
$ stringtie -p 16 -G Ref/hg19/grch37_tran/Homo_sapiens.GRCh37.75.gtf -o Asm/read1.gtf -l prefix Map/read1.bam 1>Asm/read1_strg_assm.log 2>&1
# 样本间转录本合并
$ stringtie --merge -p 16 -G Ref/hg19/grch37_tran/Homo_sapiens.GRCh37.75.gtf -o Asm/merge.gtf Asm/mergelist.txt 1>Asm/strg_merge.log 2>&1
```
`Transcript merge usage mode: stringtie --merge [Options] { gtf_list | strg1.gtf ...}`

> - -p number of threads (CPUs) to use
> - -G reference annotation to include in the merging (GTF/GFF3)
> - -o output path/file name for the assembled transcripts GTF
> - -l name prefix for output transcripts (default: STRG)

<h4 name="stringtie-quant">stringtie定量</h4>

1. 以ballgown格式输出

```
$ stringtie -e -B -p 16 -G Asm/merge.gtf -o quant/read1/read1.gtf \
	Map/read1.bam 1>quant/read1/read1_strg_quant.log 2>&1
```
> - -e only estimate the abundance of given reference transcripts (requires -G)
> - -B enable output of Ballgown table files which will be created in the same directory as the output GTF (requires -G, -o recommended)

2. 以read count进行定量，作为DESeq2或edgeR的输入

```
$ stringtie -e -p 16 -G Asm/merge.gtf -o quant/read1/read1.gtf \
	Map/read1.bam 1>quant/read1/read1_strg_quant.log 2>&1
$ python prepDE.py -i sample_lst.txt
```
注意：
- stringtie的用法与上面相同，除了少了一个参数`-B`
- `prepDE.py`脚本需要到stringtie官网下载：http://ccb.jhu.edu/software/stringtie/dl/prepDE.py

prepDE.py参数
> - -i the parent directory of the sample sub-directories or a textfile listing the paths to GTF files [default: ballgown]

samplelist textfile 格式如下：
```
ERR188021 <PATH_TO_ERR188021.gtf>
ERR188023 <PATH_TO_ERR188023.gtf>
ERR188024 <PATH_TO_ERR188024.gtf>
ERR188025 <PATH_TO_ERR188025.gtf>
ERR188027 <PATH_TO_ERR188027.gtf>
ERR188028 <PATH_TO_ERR188028.gtf>
ERR188030 <PATH_TO_ERR188030.gtf>
ERR188033 <PATH_TO_ERR188033.gtf>
ERR188034 <PATH_TO_ERR188034.gtf>
ERR188037 <PATH_TO_ERR188037.gtf>
```

<a name="diff-exp"><h3 >差异表达分析 [<sup>目录</sup>](#content)</h3></a>

---

<h4 name="deseq2"><li>使用DESeq2进行差异分析</li></h4>

DESeq2要求输入的表达矩阵是**read counts**

构建 DESeqDataSet 对象

```
# 构建表达矩阵
count_table<-read.delim("count_multicov.txt",header=F)
count_matrix<-as.matrix(count_table[,c(-1,-2,-3,-4)])
rownames(count_matrix)<-count_table$V4
# 构建分组矩阵
group_list<-factor(c("control","control","control","treat","treat","treat"))
colData <- data.frame(row.names=colnames(count_matrix), group_list=group_list))
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
	colData = colData,
	design = ~ group_list)
```
表达矩阵的格式为：

| ` ` | SRR1039508 | SRR1039509 | SRR1039512 | SRR1039513 | SRR1039516  | SRR1039517 | SRR1039520 | SRR1039521 |
|:---------:|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|
| ENSG00000000003 |	679 | 448 |	873 | 408 |	1138 |	1047 |	770 |	572 |
| ENSG00000000005 | 0 |	0 |	0 |	0 |	0 |	0 |	0 |	0 |
| ENSG00000000419 |	467 |	515 |	621 |	365 |	587 |	799 |	417  |	508 |
| ENSG00000000457 |	260 |	211 |	263 |	164 |	245 |	331 |	233  |	229 |
| ENSG00000000460 |	60 | |	55 |	40 |	35 |	78 |	63 |	76 |	60 |
| ENSG00000000938 |	0 |	0 |	2 |	0 |	1 |	0 |	0 |	0 |
| ENSG00000000971 |	3251 |	3679 |	6177 |	4252 |	6721 |	11027 |	5176 |	7995 |
| ENSG00000001036 |	1433 |	1062 |	1733 |	881 |	1424 |	1439 |	1359  |	1109 |
| ENSG00000001084 |	519 |	380 |	595 |	493 |	820 |	714 |	696  |	704 |
| ENSG00000001167 |	394 |	236 |	464 |	175 |	658 |	584 |	360  |	269 |
| ENSG00000001460 |	172 |	168 |	264 |	118 |	241 |	210 |	155  |	177 |
| ENSG00000001461 |	2112 |	1867 |	5137 |	2657 |	2735 |	2751 |	2467 |	2905
| ENSG00000001497 |	524 |	488 |	638 |	357 |	676 |	806 |	493 |	475 |
| ENSG00000001561 |	71 |	51 |	211 |	156 |	23 |	38 |	134 |	172 |

分组矩阵的格式为：

| ` ` | dex | SampleName | cell |
|:---:|:---:|:---:|:---:|
|SRR1039508|untrt|GSM1275862|N61311|
|SRR1039509|trt|GSM1275863|N61311|
|SRR1039512|untrt|GSM1275866|N052611|
|SRR1039513|trt|GSM1275867|N052611|
|SRR1039516|untrt|GSM1275870|N080611|
|SRR1039517|trt|GSM1275871|N080611|
|SRR1039520|untrt|GSM1275874|N061011|
|SRR1039521|trt|GSM1275875|N061011|

差异表达分析

```
## 数据过滤
# dds <- dds[ rowSums(counts(dds)) > 1 ,]
# dim(dds)
dds <- DESeq(dds)
plotDispEsts(dds, main="Dispersion plot")
rld <- rlogTransformation(dds)
exprMatrix_rlog=assay(rld)
res <- results(dds, contrast=c("condition",'Day1','Day0'))
resOrdered <- res[order(res$padj),] 
res_Day1_Day0=as.data.frame(resOrdered)
```


<h4 name="ballgown"><li>使用Ballgown进行差异分析</li></h4>

紧接着stringtie的定量结果进行分析

stringtie的定量结果提供多种表达值的表示方法，有read counts, RPKM/FPKM, TPM 。其中read counts是原始reads计算，RPKM/FPKM 和 TPM 都是基因表达值的归一化后的，因为本身的某些缺点，主流科学家强烈要求它就被TPM取代了。

下面对 RPKM/FPKM 和 TPM 的说明摘抄自健明大大的简书：https://www.jianshu.com/p/e9d5d7206124 ，说得通俗易懂

> **TPM是什么？**
>
> 我不喜欢看公式，直接说事情，我有一个基因A，它在这个样本的转录组数据中被测序而且mapping到基因组了 5000个的reads，而这个基因A长度是10K，我们总测序文库是50M，所以这个基因A的RPKM值是 5000除以10，再除以50，为10. 就是把基因的reads数量根据基因长度和样本测序文库来normalization 。
>
> 那么它的TPM值是多少呢？ 这个时候这些信息已经不够了，需要知道该样本其它基因的RPKM值是多少，加上该样本有3个基因，另外两个基因的RPKM值是5和35，那么我们的基因A的RPKM值为10需要换算成TPM值就是 1,000,000 \*10/(5+10+35)=200,000， 看起来是不是有点大呀！其实主要是因为我们假设的基因太少了，一般个体里面都有两万多个基因的，总和会大大的增加，这样TPM值跟RPKM值差别不会这么恐怖的。



- 载入R包
```
require(ballgown)
require(dplyr)
require(genefilter)
setwd("/share/disk5/lianm")
```
- 载入stringtie输出的表达数据，并设置表型信息（即分组信息）
```
bg_Z1Z4<-ballgown(dataDir="SingleCell_process/Cleandata/Expression",samplePattern=“Z[14]T",meas="FPKM")
pData(bg_Z1Z4)<-data.frame(id=sampleNames(bg_Z1Z4),group=c(rep(1,num_group1),rep(0,num_group2)))
```
- 过滤低丰度的基因
```
bg_Z1Z4_filt<-subset(bg_Z1Z4,"rowVars(texpr(bg_Z1Z4))>1",genomesubset=T)
```
- 差异表达基因分析
```
result_genes<-stattest(bg_Z1Z4_filt,feature="gene",covariate="group",getFC=T)
result_genes<-data.frame(geneNames=geneNames(bg_Z1Z4_filt)[match(result_genes$id,geneIDs(bg_Z1Z4_filt))],geneIDs=geneIDs(bg_Z1Z4_filt)[match(result_genes$id,geneIDs(bg_Z1Z4_filt))],result_genes)
result_genes_sort<-arrange(result_genes,pval)
write.csv(result_genes_sort,file=paste("SingleCell_process/Cleandata/Expression/",name_group1,"_VS_",name_group2,"_geneDiff_results.csv",sep=""),row.names=F)
```
> 分组设置对差异表达分析的影响：
> - FC = group_1/group_0，所以分组标签互换后FC会变为原来的倒数

- 差异转录本分析
```
result_trans<-stattest(bg_Z1Z4_filt,feature="transcript",covariate="group",getFC=T)
result_trans<-data.frame(geneNames=geneNames(bg_Z1Z4_filt),geneIDs=geneIDs(bg_Z1Z4_filt),result_trans)
result_trans_sort<-arrange(result_transs,pval)
write.csv(result_trans_sort,file=paste("SingleCell_process/Cleandata/Expression/",name_group1,"_VS_",name_group2,"_transDiff_results.csv",sep=""),row.names=F)
```

---

参考资料：

(1) [【PANDA姐的转录组入门】拟南芥实战项目-定量、差异表达、功能分析](https://mp.weixin.qq.com/s/74xmKy30ecQ2ciCqOq8MdQ)

(2) [一步法差异分析](https://github.com/jmzeng1314/my-R/tree/master/DEG_scripts)

(3) Pertea M, Kim D, Pertea G M, et al. Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown.[J]. Nature Protocols, 2016, 11(9):1650.
