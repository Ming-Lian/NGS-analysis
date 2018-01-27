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

<h1 name="title">Analysis pipeline for RNA-seq</h1>

<h3 name="download">测序数据下载 [<sup>目录</sup>](#content)</h3>

参见： https://github.com/Ming-Lian/Memo/blob/master/ChIP-seq-pipeline.md#get-data

<h3 name="salmon">Salmon流程 [<sup>目录</sup>](#content)</h3>

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

<h3 name="subread">subread流程 [<sup>目录</sup>](#content)</h3>

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

<h3 name="hisat2-stringtie">hisat2-stringtie流程 [<sup>目录</sup>](#content)</h3>

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
> - -G    reference annotation to include in the merging (GTF/GFF3)
> - -o output path/file name for the assembled transcripts GTF
> - -l name prefix for output transcripts (default: STRG)

<h4 name="stringtie-quant">定量并以ballgown格式输出</h4>

```
$ stringtie -e -B -p 16 -G Asm/merge.gtf -o quant/read1/read1.gtf \
	Map/read1.bam 1>quant/read1/read1_strg_quant.log 2>&1
```
> - -e only estimate the abundance of given reference transcripts (requires -G)
> -  -B enable output of Ballgown table files which will be created in the same directory as the output GTF (requires -G, -o recommended)

<h3 name="diff-exp">差异表达分析 [<sup>目录</sup>](#content)</h3>

<h4 name="deseq2">使用DESeq2进行差异分析</h4>

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



---

参考资料：

(1) [【PANDA姐的转录组入门】拟南芥实战项目-定量、差异表达、功能分析](https://mp.weixin.qq.com/s/74xmKy30ecQ2ciCqOq8MdQ)

(2) [一步法差异分析](https://github.com/jmzeng1314/my-R/tree/master/DEG_scripts)
