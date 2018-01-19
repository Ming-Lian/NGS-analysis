<a name="content">目录</a>

[Analysis pipeline for MeRIP-seq](#title)
- [比对参考基因组](#map)
	- [Tophat](#tophat)
	- [HISAT2](#hisat2)
- [Peak calling](#peak)
	- [MACS2](#macs2)
	- [PeakRanger](#peakranger)
- [Peaks注释](#peak-anno)
	- [CEAS](#ceas)
- [Motif识别](#motif)
	- [HOMER](#homer)
	- [MEME](#meme) 
- [Differential binding analysis](#diff-bind)
	- [Merge peaks](#merge-peaks)
	- [Preparing ChIP-seq count table](#cout-table)
	- [Differential binding by DESeq2](#deseq2)





<a name="title"><h1>Analysis pipeline for MeRIP-seq</h1></a>

![](/m6A_seq.jpg "MeRIP-seq")

<p align="center"><strong>Schematic diagram of the MeRIP-seq protocol</strong></p>

由于m6A-seq数据分析的原理与过程和ChIP-seq十分相似，所以这里略过前面的质控，简单说明比对和peak calling步骤，具体内容可以参考[**ChIP-seq分析流程**](https://github.com/Ming-Lian/Memo/blob/master/ChIP-seq-pipeline.md)

<a name="map"><h3>比对参考基因组 [<sup>目录</sup>](#content)</h3></a>

---
在 ChIP-seq 中一般用 BWA 或者 Bowtie 进行完全比对就可以了，但是在 MeRIP-seq 中，由于分析的 RNA ，那么就存在**可变剪切**，对于存在可变剪切的 mapping 用 **Tophat** 或者 Tophat 的升级工具 **HISAT2** 更合适

<a name="tophat"><h4><u>Tophat</u></h4></a>

```
# build reference index
## <1> build genome index
$ bowtie2-build  hg19.fa hg19
## <2> build transcriptome index
$ tophat -p 8 -G hg19.gtf --transcriptome-index=Ref/hg19/hg19_trans/know hg19

# mapping
$ tophat -p 8 --transcriptome-index=Ref/hg19/hg19_trans/know -o outdir hg19 reads1_1.fastq reads1_2.fastq
# Only uniquely mapped reads with mapping quality score ≥20 were kept for the subsequent analysis for each sample
$ samtools view -q 20 -O bam -o outdir/accepted_hits.highQual.bam outdir/accepted_hits.bam
```
Tophat参数
> - -p Number of threads to use
> - -G Supply TopHat with a set of gene model annotations and/or known transcripts, as a GTF 2.2 or GFF3 formatted file. 
> - --transcriptome-index TopHat should be first run with the -G/--GTF option together with the --transcriptome-index option pointing to a directory and a name prefix which will indicate where the transcriptome data files will be stored. Then subsequent TopHat runs using the same --transcriptome-index option value will directly use the transcriptome data created in the first run (no -G option needed after the first run). 
> - -o Sets the name of the directory in which TopHat will write all of its output

samtools view 参数
> - -q only include reads with mapping quality >= INT [0]

<a name="hisat2"><h4><u>HISAT2</u></h4></a>

```
# build reference index
##  using the python scripts included in the HISAT2 package, extract splice-site and exon information from the gene
annotation fle
$ extract_splice_sites.py chrX_data/genes/chrX.gtf >chrX.ss
$ extract_exons.py chrX_data/genes/chrX.gtf >chrX.exon
##  build a HISAT2 index
$ hisat2-build --ss chrX.ss --exon chrX.exon chrX_data/genome/chrX.fa chrX_tran

# mapping
$ hisat2 -p 10 --dta -x chrX_tran -1 reads1_1.fastq -2 reads1_2.fastq | samtools sort -@ 8 -O bam -o reads1.sort.bam 1>map.log 2>&1
```
`Usage: hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]`
> - -p Number of threads to use
> - --dta reports alignments tailored for transcript assemblers
> - -x Hisat2 index
> - -1 The 1st input fastq file of paired-end reads
> - -2 The 2nd input fastq file of paired-end reads
> - -S File for SAM output (default: stdout)

<a name="peak"><h3>Peak calling [<sup>目录</sup>](#content)</h3></a>

---

<a name="macs2"><h4><u>MACS2</u></h4></a>

参考ChIP-seq分析流程中的[peak calling](https://github.com/Ming-Lian/Memo/blob/master/ChIP-seq-pipeline.md#peak-calling)过程

<a name="peakranger"><h4><u>PeakRanger</u></h4></a>

```
peakranger ccat --format bam SRR1042594.sorted.bam SRR1042593.sorted.bam  \
Xu_MUT_rep1_ccat_report --report --gene_annot_file hg19refGene.txt -q 0.05 -t 4 
```

<a name="peak-anno"><h3>Peaks注释 [<sup>目录</sup>](#content)</h3></a>

---

<a name="ceas"><h4><u>CEAS</u></h4></a>

哈佛刘小乐实验室出品的软件，可以跟MACS软件call到的peaks文件无缝连接，实现peaks的注释以及可视化分析

CEAS需要三种输入文件：
> - Gene annotation table file (sqlite3)
> > 可以到CEAS官网上下载：http://liulab.dfci.harvard.edu/CEAS/src/hg18.refGene.gz ，也可以自己构建：到UCSC上下载，然后用`build_genomeBG`脚本转换成split3格式
> - BED file with ChIP regions (TXT)
> > 需要包含chromosomes, start, and end locations，这样文件可以由 peak-caller （如MACS2）得到
> - WIG file with ChiP enrichment signal (TXT)
> >  如何得到wig文件可以参考[samtools操作指南：以WIG文件输出测序深度](https://github.com/Ming-Lian/NGS-analysis/blob/master/samtools%E6%93%8D%E4%BD%9C%E6%8C%87%E5%8D%97.md#output-wig)

CEAS的使用方法很简单：
```
ceas --name=H3K36me3_ceas --pf-res=20 --gn-group-names='Top 10%,Bottom 10%'  \
-g hg19.refGene -b  ../paper_results/GSM1278641_Xu_MUT_rep1_BAF155_MUT.peaks.bed \
-w ../rawData/SRR1042593.wig
```
> - --name Experiment name. This will be used to name the output files.
> - --pf-res Wig profiling resolution, DEFAULT: 50bp
> - --gn-group-names The names of the gene groups in --gn-groups. The gene group names are separated by commas. (eg, --gn-group-names='top 10%,bottom 10%').
> - -g Gene annotation table
> - -b BED file of ChIP regions
> - -w WIG file for either wig profiling or genome background annotation.

<a name="motif"><h3>Motif识别</h3></a>
---

<a name="homer"><h4><u>HOMER</u></h4></a>

安装旧版本的HOMER比较复杂，因为旧版依赖于调用其他几个工具：
> - blat
> - Ghostscript
> - weblogo
> **Does NOT work with version 3.0!!!!**

新版HOMER安装很简单，主要是通过`configureHomer.pl`脚本来安装和管理HOMER
```
cd ~/biosoft
mkdir homer &&  cd homer
wget http://homer.salk.edu/homer/configureHomer.pl

# Installing the basic HOMER software
perl configureHomer.pl -install

# Download the hg19 version of the human genome
perl configureHomer.pl -install hg19
```
安装好后可以进行 Motif Identification
```
# 提取对应的列给HOMER作为输入文件
# change 
#		chr1	1454086	1454256	MACS_peak_1	59.88 
#to   
#		MACS_peak_1	chr1	1454086	1454256	+
$ awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' macs_peaks.bed >homer_peaks.bed

# MeRIP-seq 中 motif 的长度为6个 nt
$ findMotifsGenome.pl homer_peaks.bed hg19 motifDir -size 200 -len 8,10,12

# 自己指定background sequences，用bedtools shuffle构造随机的suffling peaks
$ bedtools shuffle -i peaks.bed -g <GENOME> >peaks_shuffle.bed
# 用参数"-bg"指定background sequences
$ findMotifsGenome.pl homer_peaks.bed hg19 motifDir -bg peaks_shuffle.bed -size 200 -len 8,10,12
```
`Usage: findMotifsGenome.pl <pos file> <genome> <output directory> [additional options]`
> 注意：<genome> 参数只需要写出genome的序号，不需要写出具体路径

最后得到的文件夹里面有一个详细的网页版报告

![](http://homer.ucsd.edu/homer/ngs/peakMotifs.output.png)

<a name="meme"><h4><u>MEME</u></h4></a>

下载安装MEME
```
cd ~/biosoft
mkdir MEMEsuite &&  cd MEMEsuite
## http://meme-suite.org/doc/download.html
wget  http://meme-suite.org/meme-software/4.11.2/meme_4.11.2_1.tar.gz
tar zxvf meme_4.11.2_1.tar.gz 
cd meme_4.11.2/
./configure --prefix=$HOME/my-bin/meme --with-url="http://meme-suite.org"
make 
make install 
```

```
# 先提取peaks区域所对应的序列
bedtools getfasta -fi input.fasta -bed input.bed -fo output.fasta
# Motif identification
meme output.fasta -dna -mod oops -pal
```
`Usage: meme <dataset> [optional arguments]`
> - < dataset > File containing sequences in FASTA format
> - -dna Sequences use DNA alphabet
> - -mod Distribution of motifs,3 options: oops | zoops | anr
> - -pal Force palindromes (requires -dna)

<a name="diff-bind"><h3>Differential binding analysis [<sup>目录</sup>](#content)</h3></a>

---

<a name="merge-peaks"><h4><u>Merge peaks</u></h4></a>

当ChIP-seq数据中有多分组，多样本以及多个重复时，需要进行样本间peaks的merge

```
bedtools intersect -a Mcf7H3k27acUcdAlnRep1_peaks.filtered.bed -b Mcf7H3k27acUcdAlnRep2_peaks.filtered.bed -wa | cut -f1-3 | sort | uniq > Mcf7Rep1_peaks.bed
bedtools intersect -a Mcf7H3k27acUcdAlnRep1_peaks.filtered.bed -b Mcf7H3k27acUcdAlnRep2_peaks.filtered.bed -wb | cut -f1-3 | sort | uniq > Mcf7Rep2_peaks.bed
bedtools intersect -a Panc1H3k27acUcdAlnRep1_peaks.filtered.bed -b Panc1H3k27acUcdAlnRep2_peaks.filtered.bed -wa | cut -f1-3 | sort | uniq > Panc1Rep1_peaks.bed
bedtools intersect -a Panc1H3k27acUcdAlnRep1_peaks.filtered.bed -b Panc1H3k27acUcdAlnRep2_peaks.filtered.bed -wb | cut -f1-3 | sort | uniq > Panc1Rep2_peaks.bed

rm *filtered*

cat *bed | sort -k1,1 -k2,2n | bedtools merge > merge.bed 
```

<a name="count-table"><h4><u>Preparing ChIP-seq count table</u></h4></a>

用**bedtools**
```
# Make a bed file adding peak id as the fourth colum
$ awk '{$3=$3"\t""peak_"NR}1' OFS="\t" merge.bed > bed_for_multicov.bed
# 输入的bam文件要提前做好index，可同时提供多个bam文件
$ bedtools multicov -bams input1.bam input2.bam ... -bed bed_for_multicov.bed > counts_multicov.txt
```
> - NR 表示awk开始执行程序后所读取的数据行数
> - OFS Out of Field Separator，输出字段分隔符

用**featureCounts** (subread工具包中的组件）
```
# Make a saf(simplified annotation format) file for featureCount in the subread package,

shown below:
GeneID	Chr	Start	End	Strand
497097	chr1	3204563	3207049	-
497097	chr1	3411783	3411982	-
497097	chr1	3660633	3661579	-
...

$ awk -F "\t" '{$1="peak_"NR FS$1;$4=$4FS"."}1' merge.bed > subread.saf
$ featureCounts -T 4 -a subread.saf -F SAF -o counts_subread.txt ../../data/*bam
```
`Usage: featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] ...`
> - -a Name of an annotation file
> - -F Specify format of the provided annotation file. Acceptable formats include 'GTF' (or compatible GFF format) and 'SAF'. 'GTF' by default
> - -o Name of the output file including read counts
> - -T Number of the threads

<a name="deseq2"><h4><u>Differential binding by DESeq2</u></h4></a>




参考资料：

(1) Zhang C, Chen Y, Sun B, et al. m(6)A modulates haematopoietic stem and progenitor cell specification[J]. Nature, 2017, 549(7671):273.

(2) Pertea M, Kim D, Pertea G, et al. Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie, and Ballgown[J]. Nature Protocols, 2016, 11(9):1650. 

(3) [ChIP-seq-pipeline](https://github.com/Ming-Lian/Memo/blob/master/ChIP-seq-pipeline.md)

(4) [ChIPseq pipeline on jmzeng1314's github](https://github.com/jmzeng1314/NGS-pipeline/tree/master/CHIPseq)

(5) [ChIPseq pipeline on crazyhottommy's github](https://github.com/crazyhottommy/ChIP-seq-analysis)
