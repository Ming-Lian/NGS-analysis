<a name="content">目录</a>

[Analysis Pipeline for Shotgun Metageomics](#title)
- [Bining: CONCOCT](#bining)
	- [Dependencies](#dependencies)
	- [Assembling Metagenomic Reads](#assemble)
	- [Cutting up contigs](#cutting-up-contigs)
	- [Map, Remove Duplicate and Quant Coverage](#map-remove-dup-quant-coverage)
	- [Run concoct](#run-concoct)







<h1 name="title">Analysis Pipeline for Shotgun Metageomics</h1>

<a name="bining"><h2>Bining: CONCOCT [<sup>目录</sup>](#content)</h2></a>

以目前主流的 Bining 工具 CONCOCT 为例

CONCOCT的算法原理，请点 [这里](Algorithms-in-Bioinformatics.md#bining-concoct)

可以使用开发者提供的测试数据，[下载地址](https://github.com/BinPro/CONCOCT-test-data/releases)

<a name="dependencies"><h3>Dependencies [<sup>目录</sup>](#content)</h3></a>

该软件的依赖：

Fundamental dependencies

```
python v2.7.*
gcc
gsl
```

Python packages

```
cython>=0.19.2
numpy>=1.7.1
scipy>=0.12.0
pandas>=0.11.0
biopython>=1.62b
scikit-learn>=0.13.1
```

Optional dependencies

```
For assembly, use your favorite, here is one
	* Vevet
		In velvet installation directory Makefile, set ‘MAXKMERLENGTH=128’, if this value is smaller in the default installation.
To create the input table (containing average coverage per sample and contig)
	* BEDTools version >= 2.15.0 (only genomeCoverageBed)
	* Picard tools version >= 1.110
	* samtools version >= 0.1.18
	* bowtie2 version >= 2.1.0
	* GNU parallel version >= 20130422
	* Python packages: pysam>=0.6
For validation of clustering using single-copy core genes
	* Prodigal >= 2.60
	* Python packages: bcbio-gff>=0.4
	* R packages: gplots, reshape, ggplot2, ellipse, getopt and grid
	* BLAST >= 2.2.28+
```

<a name="assemble"><h3>Assembling Metagenomic Reads [<sup>目录</sup>](#content)</h3></a>

```
# 将多个样本的测序数据fastq文件，按照双端分别进行合并
$ cat $CONCOCT_TEST/reads/Sample*_R1.fa > All_R1.fa
$ cat $CONCOCT_TEST/reads/Sample*_R2.fa > All_R2.fa

# 拼接
$ velveth velveth_k71 71 -fasta -shortPaired -separate All_R1.fa All_R2.fa
$ velvetg velveth_k71 -ins_length 400 -exp_cov auto -cov_cutoff auto
```

velveth:

> takes in a number of sequence files, produces a hashtable, then outputs two files in an output directory (creating it if necessary), Sequences and Roadmaps, which are necessary to velvetg. 
> 
> 语法：
> 
> ```
> ./velveth output_directory hash_length [[-file_format][-read_type] filename]
> ```

velvetg:

> Velvetg is the core of Velvet where the de Bruijn graph is built then manipulated

<a name="cutting-up-contigs"><h3>Cutting up contigs [<sup>目录</sup>](#content)</h3></a>

将大片段的contigs (>=20kb)，切成一个个10kb的小片段，当切到尾部只剩不到20kb时，停止切割，以防切得过碎

```
python $CONCOCT/scripts/cut_up_fasta.py -c 10000 -o 0 -m contigs/velvet_71.fa > contigs/velvet_71_c10K.fa
```
<a name="map-remove-dup-quant-coverage"><h3>Map, Remove Duplicate and Quant Coverage [<sup>目录</sup>](#content)</h3></a>

1. 使用 Bowtie2 执行 mapping 操作

2. 用 MarkDuplicates（Picard中的一个工具） 去除 PCR duplicates

3. 用 BEDTools genomeCoverageBed 基于 mapping 得到的 bam 文件计算每个contigs的coverage

**（1） Map, Remove Duplicate**

其中1、2步操作可以由CONCOCT中提供的脚本`map-bowtie2-markduplicates.sh`完成

先要自行建好这些contigs的bowtie2索引

```
# index for contigs
$ bowtie2-build contigs/velvet_71_c10K.fa contigs/velvet_71_c10K.fa
```
用`map-bowtie2-markduplicates.sh`脚本完成 `mapping` -> `remove duplicate`

```
for f in $CONCOCT_TEST/reads/*_R1.fa; do
    mkdir -p map/$(basename $f);
    cd map/$(basename $f);
    bash $CONCOCT/scripts/map-bowtie2-markduplicates.sh -ct 1 -p '-f' $f $(echo $f | sed s/R1/R2/) pair $CONCOCT_EXAMPLE/contigs/velvet_71_c10K.fa asm bowtie2;
    cd ../..;
done
```

> - `-c` option to compute coverage histogram with genomeCoverageBed
> - `-t` option is number of threads
> - `-p` option is the extra parameters given to bowtie2. In this case -f
> - `-k` 保留中间文件

随后的5个参数：

> - pair1, the fasta/fastq file with the #1 mates
> - pair2, the fasta/fastq file with the #2 mates
> - pair_name, a name for the pair used to prefix output files
> - assembly, a fasta file of the assembly to map the pairs to
> - assembly_name, a name for the assembly, used to postfix outputfiles
> - outputfolder, the output files will end up in this folder

如果要自己逐步执行第1、2两步，则可以通过以下方式实现：

```
# Index reference, Burrows-Wheeler Transform
$ bowtie2-build SampleA.fasta SampleA.fasta

# Align Paired end, sort and index
bowtie2 \
	-p 32 \
	-x SampleA.fasta \
	-1 $Data/SampleA.1.fastq \
	-2 $Data/SampleA.2.fastq | \
	samtools sort -@ 18 -O BAM -o SampleA.sort.bam
samtools index SampleA.sort.bam

# Mark duplicates and index
java -Xms32g -Xmx32g -XX:ParallelGCThreads=15 -XX:MaxPermSize=2g -XX:+CMSClassUnloadingEnabled \
    -jar picard.jar MarkDuplicates \
    I=./SampleA.sort.bam \
    O=./SampleA.sort.md.bam \
    M=./SampleA.smd.metrics \
    VALIDATION_STRINGENCY=LENIENT \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
    REMOVE_DUPLICATES=TRUE # 该参数默认为false，即在输出中不过滤duplicate，但是会对这些记录的flag进行修改标记
samtools index ./SampleA.sort.md.bam
```

**（2）Quant Coverage**

第3步，计算每个contigs的coverage，用`gen_input_table.py`脚本

```
# usage: gen_input_table.py [-h] [--samplenames SAMPLENAMES] [--isbedfiles] fastafile bamfiles [bamfiles ...]
# --samplenames 写有样品名的文件，每个文件名一行
# --isbedfiles  如果在上一步map时运行了genomeCoverageBed，则可以加上此参数后直接用 *smds.coverage文件。如果没运行genomeCoverageBed，则不加此参数，依旧使用bam文件。

$ python $CONCOCT/scripts/gen_input_table.py --isbedfiles \
	--samplenames <(for s in Sample*; do echo $s | cut -d'_' -f1; done) \
	../contigs/velvet_71_c10K.fa */bowtie2/asm_pair-smds.coverage \
	> concoct_inputtable.tsv
```

注：

> 这个脚本可以接受两种类型的输入
> - （1）对bamfiles执行`genomeCoverageBed (bedtools genomecov`得到的`*smds.coverage`文件，此时要使用`--isbedfiles`参数，这样脚本只执行下面提到的第2步操作——计算每条contig的平均depth（又称为这条contig的abundance）；
> - （2）原始的bamfiles，则脚本要执行下面提到的两步操作；

也可以自己写命令逐步实现，这样有利于加深对工具的理解

1. 计算每条contig的depth分布（histograms）

	<p align="center"><img src=./picture/Metagenome-Tools-CONCOCT-genomecov-1.png width=800 /></p>

	```
	$ bedtools genomecov -ibam ./SampleA.smds.bam > ./SampleA.smds.coverage
	```
	
	`bedtools genomecov`默认计算histograms，如输出为`chr1   0  980  1000`，则说明在contig chr1上depth=0的碱基数为980bp，该contig长度为1000bp
	
	例如：

	> ```
	> $ cat A.bed
	> chr1  10  20
	> chr1  20  30
	> chr2  0   500
	> 
	> $ cat my.genome
	> chr1  1000
	> chr2  500
	> 
	> $ bedtools genomecov -i A.bed -g my.genome
	> chr1   0  980  1000  0.98
	> chr1   1  20   1000  0.02
	> chr2   1  500  500   1
	> genome 0  980  1500  0.653333
	> genome 1  520  1500  0.346667
	> ```
	> 
	> 输出格式为：
	> 
	> - chromosome
	> - depth of coverage from features in input file
	> - number of bases on chromosome (or genome) with depth equal to column 2
	> - size of chromosome (or entire genome) in base pairs
	> - size of chromosome (or entire genome) in base pairs

2. 计算每条contig的平均depth

	<p align="center"><img src=./picture/Metagenome-Tools-CONCOCT-genomecov-2.png width=500 /></p>

	有两种计算方法：

	<p align="center"><img src=./picture/Metagenome-Tools-CONCOCT-genomecov-3.png width=400 /></p>

	或

	<p align="center"><img src=./picture/Metagenome-Tools-CONCOCT-genomecov-4.png width=400 /></p>

	第二种计算方法本质上就是加权平均

	```
	awk 'BEGIN {pc=""} 
	{
	    c=$1;
	    if (c == pc) {
	        cov=cov+$2*$5;
	    } else {
	      print pc,cov;
	      cov=$2*$5;
	    pc=c}
	} END {print pc,cov}' SampleA.smds.coverage | tail -n +2 > SampleA.smds.coverage.percontig
	```

**（3）Generate linkage table**

接着要构建 linkage per sample between contigs，**目前不是很理解它这一步的目的**

```
# usage: bam_to_linkage.py [-h] [--samplenames SAMPLENAMES] [--regionlength REGIONLENGTH] [--fullsearch] [-m MAX_N_CORES] [--readlength READLENGTH] [--mincontiglength MINCONTIGLENGTH] fastafile bamfiles [bamfiles ...]
# --samplenames 写有样品名的文件，每个文件名一行
# --regionlength contig序列中用于linkage的两端长度 [默认 500]
# --fullsearch 在全部contig中搜索用于linkage
# -m 最大线程数，每个ban文件对应一个线程
# --readlength untrimmed reads长度 [默认 100]
# --mincontiglength 识别的最小contig长度 [默认 0]

cd $CONCOCT_EXAMPLE/map
python bam_to_linkage.py -m 8 --regionlength 500 --fullsearch --samplenames sample.txt $DATA/SampleA.fasta ./SampleA.smds.bam > SampleA_concoct_linkage.tsv
mv SampleA_concoct_linkage.tsv ../concoct-input

# 输出文件格式
# 共2+6*i列 (i样品数)，依次为contig1、contig2、nr_links_inward_n、nr_links_outward_n、nr_links_inline_n、nr_links_inward_or_outward_n、read_count_contig1_n、read_count_contig2_n
# where n represents sample name. 
# Links只输出一次，如 contig1contig2 输出，则 contig2contig1 不输出

# contig1: Contig linking with contig2
# contig2: Contig linking with contig1
# nr_links_inward: Number of pairs confirming an inward orientation of the contigs -><-
# nr_links_outward: Number of pairs confirming an outward orientation of the contigs <--> 
# nr_links_inline: Number of pairs confirming an outward orientation of the contigs ->->
# nr_links_inward_or_outward: Number of pairs confirming an inward or outward orientation of the contigs. This can be the case if the contig is very short and the search region on both tips of a contig overlaps or the --fullsearch parameter is used and one of the reads in the pair is outside
# read_count_contig1/2: Number of reads on contig1 or contig2. With --fullsearch read count over the entire contig is used, otherwise only the number of reads in the tips are counted.
```


<a name="run-concoct"><h3>Run concoct [<sup>目录</sup>](#content)</h3></a>



参考资料：

(1) [CONCOCT’s documentation](http://concoct.readthedocs.io/en/latest/index.html)

(2) [Manual for Velvet](https://www.ebi.ac.uk/~zerbino/velvet/Manual.pdf)

(3) [BEDtools官网](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html)

(4) [【Yue Zheng博客】宏基因组binning-CONCOCT](http://www.zhengyue90.com/?p=182)
