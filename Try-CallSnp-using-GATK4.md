<a name="content">目录</a>

[跑GATK4流程](#title)
- [1. 数据说明](#description-of-testdata)
- [2. 质控与数据预处理](#qc-preprocessing)
- [3. Mapping](#mapping)
- [4. 安装GATK4与数据准备](#install-gatk4)
- [5. 排序及标记重复](#sort-and-mark-duplicate)
- [6. 质量值校正](#recallbrate-base-quality-scores)
- [7. SNP、 INDEL位点识别](#snp-indel-identify)
	- [7.1. VCF格式](#vcf-format)
- [8. 变异位点过滤](#select-filt-variants)
- [9. 变异位点注释](#variants-annotation)

<h1 name="title">跑GATK4流程</h1>

该测试在`/share/disk5/lianm/GATK4_practice`文件夹下进行

```
|-GATK4_practice
	|----------fastq	# 保存fastq文件，包括原始fastq文件和预处理后得到的fastq文件
	|----------fastqc	# FastQC质控的输出结果
	|----------Ref	# 参考序列(hg19:chr17)以及索引（BWA索引）
	|----------sam	# 保存数据处理过程中产生的sam/bam文件
```

<a name="description-of-testdata"><h2>1. 数据说明 [<sup>目录</sup>](#content)</h2></p>

从一个人的 Pair-end 重测序数据集中，抽出chr17染色体7400000-7800000的序列，得到测试用的fastq文件

- T_1.fastq
- T_2.fastq

参考基因组选择hg19

<a name="qc-preprocessing"><h2>2. 质控与数据预处理 [<sup>目录</sup>](#content)</h2></p>

```
$ fastqc -f fastq T_1.fastq T_2.fastq -o fastqc/

$ cutadpat -q 20,20 -m 20 -o T_1.QC.fastq -p T_2.QC.fastq T_1.fastq T_2.fastq
```

<a name="mapping"><h2>3. Mapping [<sup>目录</sup>](#content)</h2></p>

```
# 先建索引
$ nohup wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz >download.log &
$ tar zxvf chromFa.tar.gz
$ cp chr17.fa /share/disk5/lianm/GATK4_practice/Ref
$ bwa index -a bwtsw chr17.fa

# 比对
$ bwa mem Ref/chr17.fa fastq/T_1.QC.fastq fastq/T_2.QC.fastq -R \
	'@RG\tID:T\tLB:T\tSM:T\tPL:ILLUMINA' > sam/T.chr17.sam
```

<a name="install-gatk4"><h2>4. 安装GATK4与数据准备 [<sup>目录</sup>](#content)</h2></p>

到GATK官网获取下载链接

```
$ nohup wget -c https://github.com/broadinstitute/gatk/releases/download/4.0.8.1/gatk-4.0.8.1.zip >gatk4_download.log &

$ unzip gatk-4.0.8.1.zip

# gatk不需要安装，只需要将其路径加入到PATH中即可
$ echo "export PATH=/share/disk5/lianm/biosoft/gatk-4.0.8.1:\$PATH" >>~/.bash_profile
$ . ~/.bash_profile
```

gatk各个工具的使用方法为

```
$ gatk [--java-options "-Xmx4G"] ToolName [GATK args]
```

若想查看某个工具的帮助文档：

```
$ gatk ToolName -h/--help
```

数据准备：

1、用Samtools为参考序列创建一个索引，这是为方便GATK能够快速地获取fasta上的任何序列做准备

```
$ samtools faidx Ref/chr17.fa	# 这一步运行很快，执行结束后会在chr17.fa所在文件夹下生成chr17.fai文件
```

2、生成`.dict`文件

```
$ gatk CreateSequenceDictionary -R Ref/chr17.fa -O Ref/chr17.dict
```

dic文件的内容：

```
@HD     VN:1.5
@SQ     SN:chr17        LN:81195210     M5:351f64d4f4f9ddd45b35336ad97aa6de     UR:file:/share/disk5/lianm/GATK4_practice/Ref/chr17.fa
```

<a name="sort-and-mark-duplicate"><h2>5. 排序及标记重复 [<sup>目录</sup>](#content)</h2></p>

```
# 排序
$ gatk SortSam -I sam/T.chr17.sam -O sam/T.chr17.sort.bam -R Ref/chr17.fa \
	-SO coordinate --CREATE_INDEX

# 标记重复，这一步会比较费时
$ gatk MarkDuplicates -I sam/T.chr17.sort.bam -O sam/T.chr17.markdup.bam \
	-M sam/T.chr17.metrics --CREATE_INDEX
```

参数说明：

```
# SortSam的参数说明：

--SORT_ORDER,-SO:SortOrder    Sort order of output file.   Required. Possible values: {

                              "queryname" (Sorts according to the readname. This will place read-pairs and other derived
                              reads (secondary and supplementary) adjacent to each other. Note that the readnames are
                              compared lexicographically, even though they may include numbers. In paired reads, Read1
                              sorts before Read2.)

                              "coordinate" (Sorts primarily according to the SEQ and POS fields of the record. The
                              sequence will sorted according to the order in the sequence dictionary, taken from from
                              the header of the file. Within each reference sequence, the reads are sorted by the
                              position. Unmapped reads whose mates are mapped will be placed near their mates. Unmapped
                              read-pairs are placed after all the mapped reads and their mates.)
                              duplicate (Sorts the reads so that duplicates reads are adjacent. Required that the
                              mate-cigar (MC) tag is present. The resulting will be sorted by library, unclipped 5-prime
                              position, orientation, and mate's unclipped 5-prime position.)
                              }

--CREATE_INDEX:Boolean        Whether to create a BAM index when writing a coordinate-sorted BAM file.  Default value:
                              false. Possible values: {true, false}

# MarkDuplicates的参数说明

--METRICS_FILE,-M:File        File to write duplication metrics to  Required.
```

<a name="recallbrate-base-quality-scores"><h2>6. 质量值校正 [<sup>目录</sup>](#content)</h2></p>

**建立较正模型**

```
$ gatk BaseRecalibrator -R Ref/chr17.fa -I sam/T.chr17.markdup.bam -O \
sam/T.chr17.recal.table --known-sites ../Ref/VCF/dbsnp_138.hg19.vcf --known-sites \
../Ref/VCF/Mills_and_1000G_gold_standard.indels.hg19.vcf --known-sites \
../Ref/VCF/1000G_phase1.indels.hg19.vcf
```

参数说明：

```
--known-sites:FeatureInput    One or more databases of known polymorphic sites used to exclude regions around known
                              polymorphisms from analysis.  This argument must be specified at least once. Required.
```


**质量校正**

```
$ gatk ApplyBQSR -R Ref/chr17.fa -I sam/T.chr17.markdup.bam -bqsr \
sam/T.chr17.recal.table -O sam/T.chr17.recal.bam
```

参数说明：

```
--bqsr-recal-file,-bqsr:File  Input recalibration table for BQSR  Required.
```

<a name="snp-indel-identify"><h2>7. SNP、 INDEL位点识别 [<sup>目录</sup>](#content)</h2></p>

- **生成gvcf文件**

	```
	$ gatk HaplotypeCaller -R Ref/chr17.fa -I sam/T.chr17.recal.bam -ERC GVCF \
	--dbsnp ../Ref/VCF/dbsnp_138.hg19.vcf -O calling/T.chr17.raw.snps.indels.vcf -L \
	chr17:7400000-7800000
	```
	
	参数说明：
	
	```
	--dbsnp,-D:FeatureInput		dbSNP file  Default value: null.
	
	--emit-ref-confidence,-ERC:ReferenceConfidenceMode
					Mode for emitting reference confidence scores  Default value: NONE.
					Possible values:
					{NONE, BP_RESOLUTION, GVCF}
	
	--intervals,-L:String		One or more genomic intervals over which to operate  This argument may be specified 0 or
					more times. Default value: null.
	```
	
	-L 规定识别突变位点的区域，如-L chr17:7400000-7800000 只识别17号染色体7400000-7800000 区域的突变位点。 全外显子组分析请用捕获区域bed文件。

- **通过gvcf检测变异**

```
$ gatk GenotypeGVCFs -R Ref/chr17.fa --dbsnp ../Ref/VCF/dbsnp_138.hg19.vcf \
--variant calling/T.chr17.raw.snps.indels.vcf -O calling/T.chr17.raw.snps.indels.genotype.vcf
```

参数说明：

```
--variant,-V:String           A VCF file containing variants  Required.
```

<a name="vcf-format"><h3>7.1. VCF格式 [<sup>目录</sup>](#content)</h3></p>

<p align="center"><img src=./picture/RunGATK4-VCF-format-1.png width=900 /></p>

<p align="center"><img src=./picture/RunGATK4-VCF-format-2.png width=900 /></p>

<p align="center"><img src=./picture/RunGATK4-VCF-format-3.png width=900 /></p>

> - **QUAL**： Phred格式(Phred_scaled)的质量值，表示在该位点存在variant的可能性；该值越高，则variant的可能性越大；
> 
> ```
> 计算方法： Phred值 = -10 * log (1-p) p为variant存在的概率; 通过计算公式可以看出值为10的表示错误概率为0.1，该位点为variant的概率为90%
> ```
> 
> - FILTER：过滤信息； GATK能使用其它的方法来进行过滤，过滤结果中通过则该值为” PASS”;若variant不可靠，则该项不为” PASS”或” .”
> - INFO ： variant的详细信息
> - FORMAT 和 T： sample的基因型的信息。代表这该名称的样品，是由BAM文件中的@RG下的 SM 标签决定的

<p align="center"><strong>FORMAT 和 T 部分</strong></p>

<p align="center"><img src=./picture/RunGATK4-VCF-format-4.png width=600 /></p>

<p align="center"><img src=./picture/RunGATK4-VCF-format-4-Genotype.png width=600 /></p>

- **GT** ：基因型（Genotype）
	- 两个数字中间用’ /'分 开，这两个数字表示双倍体的sample的基因型。
	- 0 表示样品中有ref的allele；
	- 1 表示样品中variant的allele；
	- 2表示有第二个variant的allele。
	- 0/0 表示sample中该位点为纯合的，和ref一致；
	- 0/1 表示sample中该位点为杂合的，有ref和variant两个基因型；
	- 1/1 表示sample中该位点为纯合的，和variant一致。

- **AD** ： （Allele Depth）每一种allele的reads覆盖度
	- 在diploid中则是用逗号分割的两个值，前者对应ref基因型，后者对应variant基因型；
- **DP** ： （Depth）该位点的覆盖度
- **GQ** ： 基因型的质量值（Genotype Quality）
	- Phred格式(Phred_scaled)的质量值，表示在该位点该基因型存在的可能性；该值越高，则Genotype的可能性越 大；
	- 计算方法： Phred值 = -10 * log (1-p) p为基因型存在的概率。
- **PL** ： 三种基因型的质量值（provieds the likelihoods of the given genotypes）
	- 这三种指定的基因型为(0/0,0/1,1/1)，这三种基因型的概率总和为1。
	- 该值越大，表明为该种基因型的可能性越小。 Phred值 = -10 * log (p) p为基因型存在的概率。

<p align="center"><strong>INFO部分</strong></p>

<p align="center"><img src=./picture/RunGATK4-VCF-format-5.png width=600 /></p>

<p align="center"><img src=./picture/RunGATK4-VCF-format-5-INFO.png width=600 /></p>

> **AC(Allele Count)** 表示该Allele的数目；
> 
> **AF(Allele Frequency)** 表示Allele的频率；
> 
> **AN(Allele Number)** 表示Allele的总数目；
> 
> **DP**： reads覆盖度。是一些reads被过滤掉后的覆盖度；
> 
> **FS**：使用Fisher’s精确检验来检测strand bias而得到的phred格式的p值。该值越小越好；
> 
> **HaplotypeScore**： Consistency of the site with at most two segregating haplotypes 单倍体得分，值越小越好；
> 
> **MQ**： RMS Mapping Quality 所有样本中比对质量的均方根。值越大越好；
> 
> **MQRankSum**： Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities，对突变位点和参考序列位点比对质量值进行秩和检验的Z值，该值越大越好；
> 
> **QD**： Variant Confidence/Quality by Depth 位点覆盖深度的可信度，值越大越好；
> 
> **ReadPosRankSum**： Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias 突变位点与序列末端的距离，距离末端越近，则越可能为假阳性位点，值越大越好；

<a name="select-filt-variants"><h2>8. 变异位点过滤 [<sup>目录</sup>](#content)</h2></p>

```
假阳性
	|--- 测序错误
	|--- 比对错误
假阴性
	|--- 覆盖度不足
	|--- 捕获不到
```

**Steps**

> - Extract the SNPs from the call set
> - Apply the filter to the SNP call set
> - Extract the Indels from the call set
> - Apply the filter to the Indel call set
> - Combine SNP and indel call set
> - Get passed call set

- **提取SNP位点**
	
	```
	$ gatk SelectVariants -R Ref/chr17.fa -V calling/T.chr17.raw.snps.indels.genotype.vcf \
	--select-type-to-include SNP -O filter/T.chr17.raw.snps.genotype.vcf
	```
	
	参数说明：
	
	```
	--select-type-to-include,-select-type:Type
	                              Select only a certain type of variants from the input file  This argument may be specified
	                              0 or more times. Default value: null. Possible values: {NO_VARIATION, SNP, MNP, INDEL,
	                              SYMBOLIC, MIXED}
	```

- **提取INDEL位点**

	```
	$ gatk SelectVariants -R Ref/chr17.fa -V calling/T.chr17.raw.snps.indels.genotype.vcf \
	--select-type-to-include INDEL -O filter/T.chr17.raw.indels.genotype.vcf
	```

- **SNP位点过滤**

	```
	$ gatk VariantFiltration -R Ref/chr17.fa -V filter/T.chr17.raw.snps.genotype.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || \
	ReadPosRankSum < -8.0" --filter-name "SNP_FILTER" -O filter/T.chr17.filter.snps.genotype.vcf
	```
	
	参数说明：
	
	```
	--filter-expression,-filter:String
	                              One or more expression used with INFO fields to filter  This argument may be specified 0
	                              or more times. Default value: null.
	
	--filter-name:String          Names to use for the list of filters  This argument may be specified 0 or more times.
	                              Default value: null.
	```

- **INDEL位点过滤**

	```
	$ gatk VariantFiltration -R Ref/chr17.fa -V filter/T.chr17.raw.indels.genotype.vcf --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < \
	-20.0" --filter-name "INDEL_FILTER" -O filter/T.chr17.filter.indels.genotype.vcf
	```

- **合并过滤后SNP、 INDEL文件**

	```
	$ gatk MergeVcfs -I filter/T.chr17.filter.snps.genotype.vcf -I \
	filter/T.chr17.filter.indels.genotype.vcf -O filter/T.chr17.filter.snps.indels.genotype.vcf
	```

- **提取PASS突变位点**

	```
	$ gatk SelectVariants -R Ref/chr17.fa -V filter/T.chr17.filter.snps.indels.genotype.vcf -O
	T.chr17.pass.snps.indels.genotype.vcf -select "vc.isNotFiltered()"
	```
	
	参数说明：
	
	```
	--selectExpressions,-select:String
	                              One or more criteria to use when selecting the data  This argument may be specified 0 or
	                              more times. Default value: null.
	```

<a name="variants-annotation"><h2>9. 变异位点注释 [<sup>目录</sup>](#content)</h2></p>


