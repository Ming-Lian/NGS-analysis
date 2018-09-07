<a name="content">目录</a>

[跑GATK4流程](#title)
- [1. 数据说明](#description-of-testdata)
- [2. 质控与数据预处理](#qc-preprocessing)
- [3. Mapping](#mapping)
- [4. 安装GATK4与数据准备](#install-gatk4)
- [5. 排序及标记重复](#sort-and-mark-duplicate)
- [6. 质量值校正](#recallbrate-base-quality-scores)
- [7. SNP、 INDEL位点识别](#snp-indel-identify)


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
	--dbsnp,-D:FeatureInput       dbSNP file  Default value: null.
	
	--emit-ref-confidence,-ERC:ReferenceConfidenceMode
								Mode for emitting reference confidence scores  Default value: NONE.
								Possible values:
								{NONE, BP_RESOLUTION, GVCF}
	
	--intervals,-L:String         One or more genomic intervals over which to operate  This argument may be specified 0 or
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
