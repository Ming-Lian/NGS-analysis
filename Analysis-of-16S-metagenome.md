<a name="content">目录</a>

[16S微生物组](#title)
- [1. 数据预处理：QIMME](#data-preprocess)
- [2. OTU聚类](#otu-cluster)
	- [2.1. UPARSE/USEARCH](#otu-cluster-using-uparse)




<h1 name="title">16S微生物组</h1>

<a name="data-preprocess"><h2>1. 数据预处理：QIMME [<sup>目录</sup>](#content)</h2></a>

首先，认识一下16S的测序数据：

> <p align="center"><img src=./picture/16S-metagenome-data-preprocess-1.png width=800 /></p>
> 
> 红色的是barcode，蓝色的是16S通用引物
> 
> 根据每个样品的barcode和通用引物，我们可以构建QIIME1所需的mapping_file，这是QIIME1运行过程中非常重要的一个文件，包括了样本ID，barcode，两端通用引物，样品分组，以及其他描述信息等。
> 
> <p align="center"><img src=./picture/16S-metagenome-data-preprocess-2.png width=800 /></p>
> 
> 该文件在样本拆分、样本筛选、分组和统计分析中都会用到。barcode在样本拆分过程中是必须的，如果利用其他方法进行样本拆分和质量控制，可以将barcode和引物部分设为空值。
> 
> mapping_file文件的格式需要用验证通过：
> 
> ```
> $ validate_mapping_file.py -m mapping_file.txt
> ```


1. 如果是双端序列，需要将其拼接成一条序列

	样本拆分前，需要**拼接双端序列**。利用QIIME1的`join_paired_ends.py`拼接双端序列，需要安装`fastq-join`

	中间相互重叠的序列的碱基质量值会以较高的一个为准进行校正
	
	```
	$ conda install fastq-join
	$ join_paired_ends.py \
		-j 30 \  # -j 最小overlap长度，可根据自己的测序数据设置
		-f SS_R1.fastq \
		-r SS_R2.fastq \
		-o SS_merged 
	```

2. **样本拆分**：根据不同样本的barcode将reads分开；

	(1) 然后在拆分样本前还需要**提取barcode序列**：

	```
	$ extract_barcodes.py \
		-f SS_merged/fastqjoin.join.fastq \
	    -m mapping_file.txt \
	    -c barcode_paired_stitched \
		--bc1_len 8 \
		--bc2_len 8 \
	    --rev_comp_bc2 \
		-o SS_barcode
	```

	`extract_barcodes.py`可以处理多种类型的barcode，其使用方法（-c）包括了：

	> - `barcode_single_end`：处理单个文件的左侧的barcode
	> - `barcode_paired_end`：处理一对文件的barcode，输出时--fastq1（-f）barcode 在前，--fastq2（-r）barcode紧接着输出
	> - `barcode_paired_stitched`：处理单个文件的两端的barcode，输出barcode序列时，左侧的barcode在前面，右侧尾部的barcode紧接着输出

	(2) 获得barcode序列文件后，就可以**对样本进行拆分**了
	
	```
	$ split_libraries_fastq.py \
		-i SS_barcode/reads.fastq \
		-o SS_fna \
		-m mapping_file.txt \
		-b SS_barcode/barcodes.fastq \
	    -q 19 \
		--max_barcode_errors 0 \
		--barcode_type 8
	```

	样本拆分结果文件：

	```
	SS_fna
	|--- histograms.txt # 统计了不同长度的序列的数量
	|--- seqs.fna # 之后QIIME1做OTU聚类所需要的文件，该文件每个序列的ID都是以Sample_n开头，n即该样品的第n个序列
	```

3. **进行质量控制**，去掉低质量序列，去除barcode和通用引物；

	由于Illumina测序reads通常3'端质量值下降明显，通常需要将两端序列拼接后再进行质量过滤，否则可能出现大量序列由于3'端序列被切掉而无法拼接到一起的现象

	这一步可以使用 cutadapt 去除引物

	```
	$ cutadapt \
		-b CCTACGGGAGGCAGCAG \
		-b GACTACHVGGGTWTCTAAT \
	    -e 0.2 \
		-o SS_trimmed SS_fna/seqs.fna
	```

<a name="otu-cluster"><h2>2. OTU聚类 [<sup>目录</sup>](#content)</h2></a>

<a name="otu-cluster-using-uparse"><h3>2.1. UPARSE/USEARCH [<sup>目录</sup>](#content)</h3></a>

上次我们介绍了如何利用QIIM1的脚本拼接双端序列、拆分样本和质控，最后得到干净的16S序列，这些序列也是UPARSE进行OTU聚类所需要的输入文件。但是还需要将fasta文件的格式为`>sample.seqid`，这样后续得到的OTU表会自动统计不同样本的OTU序列数量

```
>sample1.0
TCCGCAATGGACGAAAGTCTG…
>sample1.1
GGACGAAAGTCTGACGGGTGC…
>sample2.0
………
```

**(1) 将序列修剪到相同的长度**

例如250bp，同时要保证序列左侧对齐，这样来自同一模板的序列可以完全匹配，后续才可以准确的计算unique序列的丰度

```
$ usearch -fastx_truncate seq.fa -trunclen 250 -fastaout seq_250.fa
```

**(2) 合并所有样本中来自同一模板的相同的序列得到unique序列，并按照数量从多到少排序**

在该过程中，丰度低于2的unique序列默认会被丢弃，这些序列很可能来自测序过程中的错误序列，会产生不可信的OTU，对后续分析造成影响

之所以需要按照数量排序，是因为丰度越高的unique序列越有可能是真实的生物序列，OTU聚类按照从丰度高到低进行搜索，以97%的一致性为标准，合并所有一致性高于该阈值的unique序列，从而得到OTU序列

```
$ usearch -fastx_uniques seq_250.fa -fastaout seq_250_unique.fa -sizeout -relabel Uniq
```

得到合并后unique序列，并且按照降序排序的文件如下：

<p align="center"><img src=./picture/16S-metagenome-OUT-cluster-1.png width=800 /></p>

**(3) 对这些unique序列进行OTU聚类**

OTU聚类大多以97%的一致性作为阈值，之所以使用这个阈值是因为这是大多数人可以接受的一个阈值，Edgar推荐的全长16S阈值为99%，V4区的阈值为100%，但如果只做属水平的分析，一般OTU阈值设为97%就可以了

```
$ usearch -cluster_otus seq_250_unique.fa -otus seq_otus.fa -relabel allOTU
```

聚类的过程如下，首先丰度最高的序列作为OTU代表序列，和其他序列进行比对，若相似性小于3%，则算作该OTU的相同序列，大于3%的序列算作新的OTU序列，同一条序列两部分来自不同OTU序列的算作嵌合体序列并且被丢掉。

<p align="center"><img src=./picture/16S-metagenome-OUT-cluster-2.jpg width=600 /></p>

通过OTU聚类得到OTU序列文件如下：

<p align="center"><img src=./picture/16S-metagenome-OUT-cluster-3.png width=800 /></p>

**(4) 统计样本OTU丰度**

通过比对原始序列文件（seq.fa）和我们得到的OTU序列库（seq_otus.fa）对每个样本的OTU数据进行统计，得到OTU表

```
$ usearch \
	-usearch_global seq.fa \
	-db seq_otus.fa \
	-id 0.97 \
	-otutabout seq_otu_table.txt \
	-biomout seq_otu_table.json
```

其中`-otutabout`和`-biomout`分别可以输出tsv和json格式的OTU表，各种OTU表格格式可以通过QIIME1的`biom convert`命令进行转换，例如将txt格式转换为hdf5格式：

```
$ biom convert -i seq_otu_table.txt -o seq_otu_table.hdf5 --table-type="OTU table" --to-hdf5
```

**(5) OTU丰度标准化**

由于每个样本的测序数量不一样，因此得到的OTU表在进行统计分析前还需要进一步标准化,这里可以用QIIME1的脚本`normalize_table.py`，或者USEARCH v11提供的命令：

```
$ usearch -otutab_rare seq_otu_table.txt -sample_size 10000 -output seq_otu_table_10k.txt 
```

我们就得到了所有样本的OTU信息，我们后续可以对该OTU表进行统计分析、进化距离分析，并且利用序列比对将OTU序列和已知的菌种16S序列关联起来，通常可以进行属水平的菌种分析。

**(6) 计算OTU多样性**

USEARCH提供了计算OTU多样性的方法，例如计算alpha多样性：

```
$ usearch -alpha_div seq_otu_table_10k.txt -output alpha.txt 
```

USEARCH提供了chao1、shannon、simpson、richness等十几种alpha多样性指数，如果只计算某几种多样性可以用-metrics参数指定：

```
$ usearch -alpha_div seq_otu_table_10k.txt -output alpha.txt  -metrics chao1,simpson
```

计算beta多样性：

```
$ usearch -beta_div seq_otu_table_10k.txt -filename_prefix ./beta
```

同样的，USEARCH提供了bray curtis, euclidean, jaccard, manhatten等多种beta多样性的计算，同样可以使用-metrics参数来指定






---

参考资料：

(1) [【宇宙实验媛】16S微生物组（一）| 数据预处理](https://mp.weixin.qq.com/s?__biz=MzU1NDkzOTk2MQ==&mid=2247483967&idx=1&sn=9e779f4e2b588ee2b81488d3fb3f7a8e&scene=21#wechat_redirect)

(2) [【宇宙实验媛】16S微生物组（二）| OTU聚类](https://mp.weixin.qq.com/s/xAifuchwB97Jv0QVhpQDwQ)

(3) [USEARCH官网](http://www.drive5.com/usearch/)
