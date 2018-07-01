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

其中1、2步操作可以由CONCOCT中提供的脚本`map-bowtie2-markduplicates.sh`完成

先要自行建好这些contigs的bowtie2索引

```
# index for contigs
$ bowtie2-build contigs/velvet_71_c10K.fa contigs/velvet_71_c10K.fa
```
用`map-bowtie2-markduplicates.sh`脚本完成 `mapping` -> `remove duplicate` -> `quantify coverage`

```
for f in $CONCOCT_TEST/reads/*_R1.fa; do
    mkdir -p map/$(basename $f);
    cd map/$(basename $f);
    bash $CONCOCT/scripts/map-bowtie2-markduplicates.sh -ct 1 -p '-f' $f $(echo $f | sed s/R1/R2/) pair $CONCOCT_EXAMPLE/contigs/velvet_71_c10K.fa asm bowtie2;
    cd ../..;
done
```

> -c option to compute coverage histogram with genomeCoverageBed
> -t option is number of threads
> -p option is the extra parameters given to bowtie2. In this case -f

随后的5个参数：

> - pair1, the fasta/fastq file with the #1 mates
> - pair2, the fasta/fastq file with the #2 mates
> - pair_name, a name for the pair used to prefix output files
> - assembly, a fasta file of the assembly to map the pairs to
> - assembly_name, a name for the assembly, used to postfix outputfiles
> - outputfolder, the output files will end up in this folder

第3步，计算每个contigs的coverage，用`gen_input_table.py`脚本

```
python $CONCOCT/scripts/gen_input_table.py --isbedfiles \
	--samplenames <(for s in Sample*; do echo $s | cut -d'_' -f1; done) \
	../contigs/velvet_71_c10K.fa */bowtie2/asm_pair-smds.coverage \
	> concoct_inputtable.tsv
```

接着要构建 linkage per sample between contigs，**目前不是很理解它这一步的目的**

```
python $CONCOCT/scripts/bam_to_linkage.py -m 8 \
	--regionlength 500 --fullsearch \
	--samplenames <(for s in Sample*; do echo $s | cut -d'_' -f1; done) \
	../contigs/velvet_71_c10K.fa Sample*/bowtie2/asm_pair-smds.bam \
	> concoct_linkage.tsv
```

<a name="run-concoct"><h3>Run concoct [<sup>目录</sup>](#content)</h3></a>



参考资料：

(1) [CONCOCT’s documentation](http://concoct.readthedocs.io/en/latest/index.html)

(2) [Manual for Velvet](https://www.ebi.ac.uk/~zerbino/velvet/Manual.pdf)
