<a name="content">目录</a>

[16S微生物组](#title)
- [数据预处理](#data-preprocess)




<h1 name="title">16S微生物组</h1>

<a name="data-preprocess"><h2>数据预处理 [<sup>目录</sup>](#content)</h2></a>

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







---

参考资料：

(1) [【宇宙实验媛】16S微生物组（一）| 数据预处理](https://mp.weixin.qq.com/s?__biz=MzU1NDkzOTk2MQ==&mid=2247483967&idx=1&sn=9e779f4e2b588ee2b81488d3fb3f7a8e&scene=21#wechat_redirect)
