<a name="content">目录</a>

[基因组结构变异 (structural variation) 检测](#title)
- [1. 基于NGS的SVs检测原理](#svs-detection-principle)
	- [1.1. Read Pair方法](#detect-method-read-pair)
	- [1.2. Split Read方法](#detect-method-split-read)
	- [1.3. Read Depth方法](#detect-method-read-depth)
	- [1.4. 基于 de novo assembly](#detect-method-denovo-assembly)
- [2. 专注CNV的检测](#focus-on-cnv-discovery)
	- [2.1. 基于芯片的检测方法](#detect-by-array)
	- [2.2. CNV-seq](#cnv-seq)
	- [2.3. Somatic CNV 检测](#scna)
	- [2.4. 实操一：CNV检测](#inaction-cnv-calling)
		- [2.4.1. GATK4 somatic CNV discovery](#gatk4)
		- [2.4.2. CNVnator](#cnvnator)
	- [2.5. 实操二：CNV区域注释基因](#inaction-cnv-annotation)
- [3. 专注SV的检测](#focus-on-sv-discovery)
	- [3.1. lumpy](#lumpy)
	- [3.2. delly](#delly)
- [补充部分](#extend)
	- [对于multiple mapping情况的处理](#extend-deal-with-multimapping)
	- [reference genome中多个相似拷贝对CNV calling的影响](#similar-copies-confound-cnv-calling)
	- [mean-shift算法](#mean-shift-algorithm)






<h1 name="title">基因组结构变异 (structural variation) 检测</h1>

<p align="center"><img src=./picture/StructralVariation-SVs-types.png width=800 /></p>

检测SVs的意义：

> - SVs对基因组的影响比起SNP更大，一旦发生往往会给生命体带来重大影响，比如导致出生缺陷、癌症等；
> - 有研究发现，基因组上的SVs比起SNP而言，更能代表人类群体的多样性特征；
> - 稀有且相同的一些结构性变异往往和疾病（包括癌症）的发生相互关联甚至还是其直接的致病诱因。

<a name="svs-detection-principle"><h2>1. 基于NGS SVs检测原理 [<sup>目录</sup>](#content)</h2></a>

目前已有的检测SVs的策略：

- Read Pair，一般称为Pair-End Mapping，简称RP或者PEM；
- Split Read，简称SR；
- Read Depth，简称RD，也有人将其称为RC——Read Count的意思，它与Read Depth是同一回事，顾名思义都是利用read覆盖情况来检测变异的方法；
- 序列从头组装（de novo Assembly， 简称AS）的方法。

<p align="center"><img src=./picture/StructralVariation-SVs-detection-principle-outline.png width=800 /></p>

<p align="center"><img src=./picture/StructralVariation-SVs-detection-principle-outline-2.jpg width=800 /></p>

<a name="detect-method-read-pair"><h3>1.1. Read Pair方法 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/StructralVariation-SVs-detection-principle-ReadPair-1.png width=800 /></p>

(1) 这个**插入片段长度的分布**是RP方法进行变异检测的一个关键信息

<p align="center"><img src=./picture/StructralVariation-SVs-detection-principle-ReadPair-2.jpg width=800 /></p>

如果插入片段长度有异常，即明显偏离分布中心的片段长度，则它可能是隐含的变异信号！

> 举个例子，如果我们发现它这个计算出来的插入片段长度与正态分布的中心相比大了200bp（假设这个200bp已经大于3个标准差了），那么就意味着参考基因组比read1和read2所在的片段要长200bp，通过类似这样的方式，我们就可以发现read1和read2所在的序列片段相比与参考基因组而言发生了200bp的删除（Deletion）。

(2) 通过比**对read1和read2之间的序列位置关系**，还能够发现更多非线性的序列变异

比如，序列倒置（Inversion），因为，按照PE的测序原理，read1和read2与参考基因组相比对，**正好是一正一负**，要么是read1比上正链，read2比上负链，要么是反过来，而且read1和read2都应处于同一个染色体上，如果不是这种现象，那么就很可能是序列的非线性结构性变异所致

RP方法存在的缺陷：

> - RP对于大长度Deletion（通常是大于1kbp）比较敏感，准确性也高，而50bp-200bp的这个范围内的变异是它的一个检测暗区。
> 
> 对于Deletion的检测，由于要求插入片段长度的变化要具有统计意义上的显著性，所以它所能检测到的片段长度就会受插入片段长度的标准差（SD）所影响。简单来说，就是越大的序列删除越偏离正常的长度中心，才越容易被检测到
> 
> - 它所能检测的Insertion序列，长度无法超过插入片段的长度
> 
> 如果这个Insertion序列很长——举个极端的例子——整个插入片段都是Insertion序列，那么你会发现read1和read2根本就不会比对上基因组，它在基因组上一点信号都没有，你甚至都不知道有这个序列的存在。另外，Insertion的检测精度也同时受限于插入片段长度的标准差。

PEM方法共有两种策略来鉴定SVs/CNVs：

- **clustering approach**：使用**预先定义的距离**（predefined distance）来检测异常（discordant）的PE reads，那这里就产生了一个问题：**如何得到合适的预定义距离？**
- **model-based approach**：对全基因范围的PE read的插入长度进行**概率分布建模**，然后基于统计检验得到异常的PE reads；




<p align="center">使用PEM方法检测CNV的工具列表</p>

|	Tool	|	URL	|	Language	|	Input	|	Comments	|
|:---|:---|:---|:---|:---|
|	BreakDancer	|	http://breakdancer.sourceforge.net/	|	Perl, C++	|	Alignment files	|	Predicting nsertions, deletions, inversions, inter- and intra-chromosomal translocations	|
|	PEMer	|	http://sv.gersteinlab.org/pemer/	|	Perl, Python	|	FASTA	|	Using simulation-based error models to call SVs	|
|	VariationHunter	|	http://compbio.cs.sfu.ca/strvar.htm	|	C	|	DIVETa	|	Detecting insertions, deletions and inversions	|
|	commonLAW	|	http://compbio.cs.sfu.ca/strvar.htm	|	C++	|	Alignment files	|	Aligning multiple samples simultaneously to gain accurate SVs using maximum parsimony model	|
|	GASV	|	http://code.google.com/p/gasv/	|	Java	|	BAM	|	A geometric approach for classification and comparison of structural variants	|
|	Spanner	|	N/A	|	N/A	|	N/A	|	Using PEM to detect tandem duplications	|

<a name="detect-method-split-read"><h3>1.2. Split Read方法 [<sup>目录</sup>](#content)</h3></a>

SR算法的核心也是对非正常PE比对数据的利用

对比上面提到的RP方法：

> RP中的非正常比对，通常是read1和read2在距离或者位置关系上存在着不正常的情形，而它的**一对PE read都是能够“无伤”地进行比对的**；
> 
> 但SR一般是指这两条PE的read，**有一条能够正常比对上参考基因组，但是另一条却不行**的情形。
> 
> 这个时候，比对软件（比如BWA）会尝试把这条没能够正常比上基因组的read在插入片段长度的波动范围内，使用更加宽松的Smith-Waterman局部比对方法，尝试搜索这条read最终可能比对得上的位置。如果这条read有一部分能够比上，那么BWA会对其进行软切除——soft-clip（CIGAR序列中包含S的那些read，图9中也有这类比对情况），标记能够比上的子序列（但不能比上的序列还是留在原read中，这也是软切除的含义）。

<p align="center"><img src=./picture/StructralVariation-SVs-detection-principle-SplitRead.png width=800 /></p>

SR的一个优势在于，它所检测到的SVs断点能**精确到单个碱基**，但是也和大多数的RP方法一样，无法解决复杂结构性变异的情形。

<p align="center">使用Split-Read方法检测CNV的工具列表</p>

|	Tool	|	URL	|	Language	|	Input	|	Comments	|
|:---|:---|:---|:---|:---|
|	AGE	|	http://sv.gersteinlab.org/age	|	C++	|	FASTA	|	A dynamic-programming algorithm using optimal alignments with gap excision to detect breakpoints	|
|	Pindel	|	http://www.ebi.ac.uk/~kye/pindel/	|	C++	|	BAM /FASTQ	|	Using a pattern growth approach to identify breakpoints of various SVs	|
|	SLOPE	|	http://www-genepi.med.utah.edu/suppl/SLOPE	|	C++	|	SAM/FASTQ/MAQb	|	Locating SVs from targeted sequencing data	|
|	SRiC	|	N/A	|	N/A	|	BLAT output	|	CalibratingSV calling using realistic error models	|


<a name="detect-method-read-depth"><h3>1.3. Read Depth方法 [<sup>目录</sup>](#content)</h3></a>

其实CNV实质上是序列Deletion或Duplication，是可以归类于Deletion和Insertion这个大的分类的，只是由于它的发生有着其独特的特点，而且往往还比较长，所以也就习惯了独立区分。

Read Depth方法一般用于CNVs的检测，其检测原理为：

> 基因组区域的Read Depth与其拷贝数相关
> 
> 全基因组测序（WGS）得到的覆盖深度呈现出来的是一个泊松分布——因为基因组上任意一个位点被测到的几率都是很低的——是一个小概率事件，在很大量的测序read条件下，其覆盖就会呈现一个泊松分布，如下图。
> 
> <p align="center"><img src=./picture/StructralVariation-SVs-detection-principle-ReadDepth-poission.png width=800 /></p>
> 
> **拷贝数增加会使得该区域的Read Depth高于期望值，而拷贝数缺失使得该区域的Read Depth低于高于期望值**

目前有两种利用Read depth信息检测CNV的策略：

- 一种是，通过检测样本在参考基因组上read的深度分布情况来发现CNV，这类适用于单样本，也是用的比较多的一个方法；

- 另一种则是通过识别并比较两个样本在基因组上存在丢失和重复倍增的区域，以此来获得彼此相对的CNV，适用于case-control模型的样本，或者肿瘤样本Somatic CNV的发现，这有点像CGH芯片的方法，因此又被称作CNV-seq，具体信息见 [2.1. CNV-seq](#cnv-seq)。

CNVnator使用的是第一种策略，同时也广泛地被用于检测大的CNV，当然还有很多冷门的软件，这里就不再列举了；CNV-seq使用的则是第二种策略。

<p align="center">使用Read-Depth方法检测CNV的工具列表</p>

|	Tool	|	URL	|	Language	|	Input	|	Comments	|
|:---|:---|:---|:---|:---|
|	SegSeqa	|	http://www.broad.mit.edu/cancer/pub/solexa_copy_numbers/	|	Matlab	|	Aligned read positions	|	Detecting CNV breakpoints using massively parallel sequence data	|
|	CNV-seqa	|	http://tiger.dbs.nus.edu.sg/cnv-seq/	|	Perl, R	|	Aligned read positions	|	Identifying CNVs using the difference of observed copy number ratios	|
|	RDXplorerb	|	http://rdxplorer.sourceforge.net/	|	Python, Shell	|	BAM	|	Detecting CNVs through event-wise testing algorithm on normalized read depth of coverage	|
|	BIC-seqa	|	http://compbio.med.harvard.edu/Supplements/PNAS11.html	|	Perl, R	|	BAM	|	Using the Bayesian information criterion to detect CNVs based on uniquely mapped reads	|
|	CNAsega	|	http://www.compbio.group.cam.ac.uk/software/cnaseg	|	R	|	BAM	|	Using flowcell-to-flowcell variability in cancer and control samples to reduce false positives	|
|	cn.MOPSb	|	http://www.bioinf.jku.at/software/cnmops/	|	R	|	BAM/read count matrices	|	Modelling of read depths across samples at each genomic position using mixture Poisson model	|
|	JointSLMb	|	http://nar.oxfordjournals.org/content/suppl/2011/02/16/gkr068.DC1/JointSLM_R_Package.zip	|	R	|	SAM/BAM	|	Population-based approach to detect common CNVs using read depth data	|
|	ReadDepth	|	http://code.google.com/p/readdepth/	|	R	|	BED files	|	Using breakpoints to increase the resolution of CNV detection from low-coverage reads	|
|	rSW-seqa	|	http://compbio.med.harvard.edu/Supplements/BMCBioinfo10-2.html	|	C	|	Aligned read positions	|	Identifying CNVs by comparing matched tumor and control sample	|
|	CNVnator	|	http://sv.gersteinlab.org/	|	C++	|	BAM	|	Using mean-shift approach and performing multiple-bandwidth partitioning and GC correction	|
|	CNVnorma	|	http://www.precancer.leeds.ac.uk/cnanorm	|	R	|	Aligned read positions	|	Identifying contamination level with normal cells	|
|	CMDSb	|	https://dsgweb.wustl.edu/qunyuan/software/cmds	|	C, R	|	Aligned read positions	|	Discovering CNVs from multiple samples	|
|	mrCaNaVar	|	http://mrcanavar.sourceforge.net/	|	C	|	SAM	|	A tool to detect large segmental duplications and insertions	|
|	CNVeM	|	N/A	|	N/A	|	N/A	|	Predicting CNV breakpoints in base-pair resolution	|
|	cnvHMM	|	http://genome.wustl.edu/software/cnvhmm	|	C	|	Consensus sequence from SAMtools	|	Using HMM to detect CNV	|


<a name="detect-method-denovo-assembly"><h3>1.4. 基于 de novo assembly [<sup>目录</sup>](#content)</h3></a>

其实从上面看下来，**SVs检测最大的难点实际上是read太短导致的**。就因为read太短，我们不能够在比对的时候横跨基因组重复区域；就因为read太短，很多大的Insertion序列根本就没能够看到信息；就因为read太短，比对才那么纠结，我们才需要用各种数学模型来 猜测这个变异到底应该是什么等等。

但从理论上来讲，三代测序和de novo assembly 的方法应该要算是基因组结构性变异检测上最有效的方法，它们都能够检测所有类型的结构性变异。

<p align="center"><img src=./picture/StructralVariation-SVs-detection-principle-denovo-assembly-1.png width=800 /></p>



<a name="focus-on-cnv-discovery"><h2>2. 专注CNV的检测 [<sup>目录</sup>](#content)</h2></a>

<a name="detect-by-array"><h3>2.1. 基于芯片的检测方法 [<sup>目录</sup>](#content)</h3></a>

最早的CNV检测实验方法，是采用**染色体核型分析 (karyotyping)**和**FISH（荧光原味杂交）**

在2003年，出现了**arrayCGH**方法和**SNP芯片**方法，这两种全基因范围的高通量检测方法，它们的检测原理是：

> 依据拷贝数与杂交信号的强度正相关，对特定的待检测基因组区段设计特异性杂交探针，将配对的Test Genome和Reference Genome分别用不同颜色的荧光标记，然后与芯片进行杂交，比较两种荧光信号的强弱

<p align="center"><img src=./picture/StructralVariation-CNV-seq-aCGH.jpg width=500 /></p>

比如SNP6.0芯片就是其中的一款代表产品，TCGA里面主要是通过Affymetrix SNP6.0 array这款芯片来测拷贝数变异，值得注意的是，并不是只有TCGA利用了SNP6.这个芯片数据，著名的CCLE计划也对一千多细胞系处理了SNP6.0芯片，数据也是可以下载的。

对SNP6.0的拷贝数芯片来说，通常是用`PICNIC`等软件处理原始数据，就可以得到的segment记录文件，每个样本一个结果，下面是示例结果：

```
Chromosome  Start   End Num_Probes  Segment_Mean
1   61735   1510801 226 -0.0397
1   1627918 1672603 17  -0.92
1   1687587 16153497    8176    0.0077
1   16153536    16153925    5   -2.7441
1   16154201    16155010    4   -0.8711
1   16165661    72768498    34630   0.0048
1   72768916    72811148    46  -1.7394
1   72811904    95674710    14901   0.0026
1   95676511    95676518    2   -1.6636
```

表明了某条染色体的某个区域内，SNP6.0芯片设计了多少个探针，芯片结果的拷贝数值是多少(这个区域的拷贝数用Segment_Mean)。通常二倍体的Segment_Mean值为0，可以用-0.2和0.2来作为该区域是否缺失或者扩增，也有人选择0.4作为阈值。

但是这些基于芯片杂交的检测方法存在一些缺点：

> - 存在杂交信号噪声；
> - 只能对已知的变异进行检测，无法发现新型的变异和罕见变异；
> - 低分辨率；
> - 有限的基因组覆盖度，有些区域因为一些原因没有设计检测探针，因此也就覆盖不到这些区域；








<a name="cnv-seq"><h3>2.2. CNV-seq [<sup>目录</sup>](#content)</h3></a>

CNV-seq：通过低倍全基因组测序检测CNV技术，2009年由Xie等开发提出

**原理：**

> CNV-seq检测原理派生自aCGH（比较基因组杂交芯片），但不针对芯片数据，而是使用NGS测序数据进行检测。将等量的待测样本DNA和正常样本DNA分别进行建库进行NGS测序后，与reference sequence进行比对，通过比较待测样本和正常样本每个滑动窗口内reads数目的多少（待测样本reads数/正常对照样本reads数）来确定待测样本每个位点的拷贝数。

其实简单来说，就是前面提到的Read Depth方法

它可以检出全基因组水平的大片段的CNV

**检测效果：**

> 采用双盲实验设计，对染色体结构异常的样本进行检测，滑动窗口大小为20kb的情况下，CNV-seq和高密度SNP-array对已知致病CNVs都能达到100%的检出。
> 
> 与中等密度SNP-array相比，CNV-seq表现更优。
> 
> 基于**PCR-free**的WGS文库构建方法，**GC偏好性性**更小，数据波动更小，表现更稳定。因为没有捕获设计（与WES比较），reads在全基因组范围**覆盖均一性**和**随机性**相比捕获测序的方法更好

<p align="center"><img src=./picture/StructralVariation-CNV-seq-vs-aCGH.jpg width=800 /></p>

<p align="center"><img src=./picture/StructralVariation-CNV-seq-advantages.jpg width=800 /></p>

<p align="center"><img src=./picture/StructralVariation-CNV-seq-performance.png width=800 /></p>

<p align="center">CNV-seq与芯片平台灵敏度比较</p>

从上图可以看到CNV-seq的检测灵敏度明显要高于芯片平台，而且**以mate-pair方式建库的NGS检测方法的检出率更高**

<a name="scna"><h3>2.3. Somatic CNV 检测 [<sup>目录</sup>](#content)</h3></a>

SCNA: somatic copy number alterations，一般指的是癌症中tumor tissue中的拷贝数变异

SCNA检测存在的难度：

- 进行tumor tissue的bulk取样测序，组织中的细胞异质性高，信号比较弱
- 可能受到正常的germline细胞的污染
- 来自父本和母本的reads混在一起




<a name="tools"><h3>2.4. CNV检测工具 [<sup>目录</sup>](#content)</h3></a>

<a name="gatk4"><h4>2.4.1. GATK4 somatic CNV discovery [<sup>目录</sup>](#content)</h4></a>

<p align="center"><img src=./picture/StructralVariation-CNV-discovery-GATK4-workflow.png width=800 /></p>

注：

> 目前还不懂用GATK4检测CNV流程中每一步的原理，只是简单地将**生物技能树**中的内容复制过来，而且近期也没有对这个分析流程进行深究的打算，仅当做一个知识的补充——毕竟GATK4并不是目前用于call SNV的主流工具

**1\. 首先制作外显子坐标记录文件**

```
# bed to intervals_list
$ cat exon_probe.hg38.gene.bed|awk '{print "chr"$0}' >hg38.chr.bed
$ java -jar ~/biosoft/picardtools/2.9.2/picard.jar  BedToIntervalList \
	I=hg38.chr.bed \
	O=list.interval_list \
	SD=/home/jianmingzeng/biosoft/GATK/resources/bundle/hg38/Homo_sapiens_assembly38.dict

# Preprocess Intervals
$ gatk  PreprocessIntervals \
	-L list.interval_list \
	--sequence-dictionary /home/jianmingzeng/biosoft/GATK/resources/bundle/hg38/Homo_sapiens_assembly38.dict \
	--reference /home/jianmingzeng/biosoft/GATK/resources/bundle/hg38/Homo_sapiens_assembly38.fasta  \
	--padding 250 \
	--bin-length 0 \
	--interval-merging-rule OVERLAPPING_ONLY \
	--output targets.preprocessed.interval.list
```

**2\. 然后把bam文件转为外显子reads数**

```
for i in  /PATH/TO/*.bam
do
	j=$(basename "$i" _recal.bam)
	echo $j
	## step1 : CollectReadCounts
	gatk  --java-options "-Xmx10G -Djava.io.tmpdir=./" CollectReadCounts \
		-I $i \
		-L $tl \
		-R $GENOME \
		--format HDF5  \
		--interval-merging-rule OVERLAPPING_ONLY \
		--output $j.clean_counts.hdf5

	## step2 : CollectAllelicCounts，这一步非常耗时，而且占空间
	gatk  --java-options "-Xmx10G -Djava.io.tmpdir=./" CollectAllelicCounts \
		-I $i \
		-L $tl \
		-R $GENOME \
		-O $j.allelicCounts.tsv
done
```

**3\. 接着合并所有的normal样本的数据创建** `cnvponM.pon.hdf5`

```
$ gatk --java-options "-Xmx20g" CreateReadCountPanelOfNormals \
	--minimum-interval-median-percentile 5.0 \
	--output cnvponM.pon.hdf5 \
	--input counts/OSCC_01_N.clean_counts.hdf5 \
	--input counts/OSCC_04_N.clean_counts.hdf5 
```

值得注意的是这个cnvponM.pon.hdf5文件， h5py文件是存放两类对象的容器，数据集(dataset)和组(group)，dataset类似数组类的数据集合，和numpy的数组差不多。group是像文件夹一样的容器，它好比python中的字典，有键(key)和值(value)。group中可以存放dataset或者其他的group。”键”就是组成员的名称，”值”就是组成员对象本身(组或者数据集)。

**5\. 最后走真正的CNV流程**

```
for i in ./counts/*
do
	j=$(basename "$i" .clean_counts.hdf5)
	gatk  --java-options "-Xmx20g" DenoiseReadCounts \
		-I $i \
		--count-panel-of-normals q \
		--standardized-copy-ratios $j.clean.standardizedCR.tsv \
		--denoised-copy-ratios $j.clean.denoisedCR.tsv
done


for i in denoisedCR/*
do
	j=$(basename "$i" .clean.denoisedCR.tsv)
	## ModelSegments的时候有两个策略，是否利用CollectAllelicCounts的结果
	gatk   --java-options "-Xmx20g" ModelSegments \
		--denoised-copy-ratios $i \
		--output segments \
		--output-prefix $j

	gatk   --java-options "-Xmx20g" CallCopyRatioSegments \
-I segments/$j.cr.seg \
-O segments/$j.clean.called.seg
```

<a name="cnvnator"><h4>2.4.2. CNVnator [<sup>目录</sup>](#content)</h4></a>

CNVnator的CNV检测原理：

> 1. 切分bins，定量bins的RD信号值；
> 
> 	首先将整个基因组切分成连续而不重叠的bins，每个bins的大小相同。然后计算每个bin的RD信号强度，由于在之前的研究中发现GC含量会影响RD信号的定量
> 	
> 	<p align="center"><img src=./picture/StructralVariation-CNV-discovery-CNVnator-GCcorrection.png width=800 /></p>
> 	
> 	因此根据GC含量对RD信号进行了校正，公式如下：
> 	
> 	<p align="center"><img src=./picture/StructuralVariation-CNV-discovery-CNVnator-GCcorrection-formula.jpg width=800 /></p>
> 
> 2. 分割不同CN的片段 (segments)。从全基因组范围的bins的RD信号分布图中，找出不同的CN (Copy Number) 的区域，确定相邻CN区域的breakpoint
> 
> 	该过程用到了图像处理领域常用的一个算法：mean shift 算法，该算法在下文补充部分有进一步的说明，点 [这里](#mean-shift-algorithm) 查看
> 
> 	<p align="center"><img src=./picture/StructuralVariation-CNV-discovery-CNVnator-meanshift.jpg width=600 /></p>
> 
> 	对于图中的每一个点（每一个点表示一个bin）获得它的mean-shift向量——指向与它最相似的那些bins（要么朝左要么朝右）；
> 
> 	根据这些bins的mean-shift向量的朝向，确定相邻CN区域的breakpoint——当相邻两个bins的mean-shift向量的朝向为背靠背形式时，breakpoint位于这两个bins之间；


```
# 1.提取mapping信息
$ cnvnator -root Sample1.root -tree Sample1.sorted.bam -unique 

# 2.生成质量分布图HISTOGRAM
$ cnvnator -root Sample1.root -his 100  -d hg38/Homo_sapiens_assembly38.fasta 

# 3.生成统计结果
$ cnvnator -root Sample1.root -stat 100 

# 4.RD信息分割partipition
$ cnvnator -root Sample1.root -partition 100 

# 5.变异检出
$ cnvnator -root Sample1.root -call 100 > Sample1.cnvnator.vcf
```

<a name="focus-on-sv-discovery"><h2>3. 专注SV的检测 [<sup>目录</sup>](#content)</h2></a>

<a name="lumpy"><h3>3.1. lumpy [<sup>目录</sup>](#content)</h3></a>

```
$ lumpyexpress \
	-B Sample1.sorted.bam \
	-S Sample1.discordants.sorted.bam \
	-D Sample1.splitters.sorted.bam \
	-o Sample1.lumpu.sv.vcf
```

<a name="delly"><h3>3.2. delly [<sup>目录</sup>](#content)</h3></a>

```
# SV检测
$ delly call \
	-g hg38/Homo_sapiens_assembly38.fasta \
	-o Sample1.delly.sv.bcf \
	-n Sample1.sorted.bam

# 过滤结果
$ delly filter \
	-f germline \
	-p \
	-q 20 Sample1.delly.sv.bcf \
	-o Sample1.delly.sv.filter.bcf
```

<a name="extend"><h2>补充部分 [<sup>目录</sup>](#content)</h2></a>

<a name="extend-deal-with-multimapping"><h3>对于multiple mapping情况的处理 [<sup>目录</sup>](#content)</h3></a>


大多数的reads都能在基因组中找到它唯一的mapping位置，但是如果这条reads是来源于重复序列区域，那么它就会得到multiple mapping的位置，且每条mapping的结果的分值都相同，此时对于这些multiple mapping的reads如何处置呢？

目前有两种主要的处理办法：

- 一种处理方式是在mapping之前，将reference中的重复序列过滤（mask）；
- 另一种方法是，从这些得分相同的multiple mapping位置中随机挑选一个，最为这条reads最后的归属；

对于paired-end的测序数据，还可以利用额外的信息来辅助multiple mapping情况的处置：

> 利用PE reads mapping的插入片段距离和双端的mapping方向，将那些距离或方向存在异常的mapping位置丢弃，以此来排除最不可能的mappning情况

<a name="similar-copies-confound-cnv-calling"><h3>reference genome中多个相似拷贝对CNV calling的影响 [<sup>目录</sup>](#content)</h3></a>

当参考基因组中存在多个相同或相似的拷贝区域时，会对CNV calling带来一定的误导

例如，在参考基因组（单倍体）中存在两个完全相同的片段重复，分别记为A和B，而在检测的样本中
只存在A拷贝，由于人是二倍体，样本会存在两份A拷贝，那么来着样本的A区域的reads在mapping过程中会随机均匀地分布在参考基因组的A区域和B区域（分值相同的multiple mapping随机指定一个mapping位置），则A区域和B区域的Read Depth为实际的一半，因此这两个基因组区域都会被鉴定为Deletion

为了解决这种情况，CNVnator提出了一种解决方案：

> 对于那些multiple mapping的reads，会被赋予一个mapping quality，为0，这些reads被称为q0 reads
> 
> 统计每个called CNV segments的q0 reads的比例，然后统计这些比例的频数分布，得到下面的直方图：
> 
> <p align="center"><img src=./picture/StructuralVariation-extend-deal-with-similar-segements-confound.png width=800 /></p>
> 
> 若一个CNV片段的q0比例大于50%，则暗示着它有可能是在参考基因组中存在高度相似的多个基因组区域，此时的CNV calling的结果可能不准确，需要进行过滤


<a name="mean-shift-algorithm"><h3>mean-shift算法 [<sup>目录</sup>](#content)</h3></a>

假设我们有一堆点，和一个小的圆形窗口，我们要完成的任务就是将这个窗口移动到最大灰度密度处（也就是点最多的地方）。如下图所示：

<p align="center"><img src=./picture/StructralVariation-extend-mean-shift-algorithm.png width=800 /></p>

初始窗口是蓝色的C1，它的圆心为蓝色方框的C1_o，而窗口中所有点质心却是C1_r，很明显圆心和点的质心没有重合。所以移动圆心C1_o到质心C1_r，这样我们就得到了一个新的窗口。这时又可以找到新的窗口内所以点的质心，大多数情况下还是不重合，所以重复上面的操作直到：将新窗口的圆心和它所包含的点的质心重合，这样我们的窗口会落在像素值（和）最大的地方。如上图C2是窗口的最后位置，它包含的像素点最多。

每一次迭代中需要移动圆心到质心，这个移动的方向称为**Mean Shift向量**

那这个Mean Shift向量怎么计算呢？

> 对于给定d维空间 R<sup>d</sup> 中的n个样本点x<sub>i</sub>,i=1,2,…,n在xd点的Mean Shift向量的基本形式定义为： 
> 
> <p align="center"><img src=./picture/StructralVariation-extend-mean-shift-algorithm-2.png width=200 /></p>
> 
> 其中，S<sub>h</sub>是一个半径为h的高维球区域：S<sub>h</sub>(x)={y:(y−x)T(y−x)≤h<sup>2</sup>}，k表示在这n个样本点中有k个落入球S<sub>h</sub>中。
> 
> 直观上来看，这k个样本点在x处的偏移向量即为：对落入Sh区域中的k个样本点相对于点x的偏移向量求和然后取平均值；
> 
> 几何解释为：如果样本点xi服从一个概率密度函数为f(x)的分布，由于非零的概率密度函数的梯度指向概率密度增加最大的方向，因此从平均上来说，Sh区域内的样本点更多的落在沿着概率密度梯度的方向。因此，Mean Shift向量Mh(x)应该指向概率密度梯度的方向

---

参考资料：

(1) [【简书】一篇文章说清楚基因组结构性变异检测的方法](https://www.jianshu.com/p/4c8e109f0e6a)

(2) [【简书】TCGA的28篇教程之CNV那点事](https://www.jianshu.com/p/eadfc45f1f18?utm_campaign=haruki&utm_content=note&utm_medium=reader_share&utm_source=weixin)

(3) [【基因学苑】好东西应该再次分享——CNV数据分析](https://mp.weixin.qq.com/s/5xuqHTjYLDd-_1jzLXEeOw)

(4) Zhou B, Ho SS, Zhang X, et al.Whole-genome sequencing analysis of CNV using low-coverage and paired-end strategies is efficient and outperforms array-based CNV analysis Journal of Medical Genetics 2018;55:735-743.

(5)  Xie C , Tammi M T . CNV-seq, a new method to detect copy number variation using high-throughput sequencing[J]. BMC Bioinformatics, 2009, 10(1):80-0.

(6) Zhao M, Wang Q, Wang Q, Jia P, Zhao Z. Computational tools for copy number variation (CNV) detection using next-generation sequencing data: features and perspectives. BMC Bioinformatics. 2013;14 Suppl 11(Suppl 11):S1. 

(7) [GATK4 best practice for somatic CNV discovery](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11147)

(8) [【github】GATK4官方somatic CNV鉴定流程](cnv_somatic_panel_workflow.wdl
)

(9) [【生信技能树】GATK4的CNV流程-hg38](https://mp.weixin.qq.com/s/Lvfy7Y352WhLMuzawMvroA)

(10) [【基因学苑】一个人全基因组完整数据分析脚本](https://mp.weixin.qq.com/s/TBOhU_4d3iPQIBDRfv3e0A)

(11) [【CSDN博客】OpenCV学习笔记-MeanShift](https://blog.csdn.net/qq_36387683/article/details/80598587)

(12) [【CSDN博客】OpenCV之均值漂移(Mean Shift)算法](https://blog.csdn.net/qq_23968185/article/details/51804574)

(13) Abyzov A, Urban AE, Snyder M, Gerstein M. CNVnator: an approach to discover, genotype, and characterize typical and atypical CNVs from family and population genome sequencing. Genome Res. 2011;21(6):974-84. 

