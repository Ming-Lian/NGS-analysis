<a name="content">目录</a>

[基因组结构变异 (structural variation) 检测](#title)
- [1. 基于NGS的SVs检测原理](#svs-detection-principle)
	- [1.1. Read Pair方法](#detect-method-read-pair)
	- [1.2. Split Read方法](#detect-method-split-read)
	- [1.3. Read Depth方法](#detect-method-read-depth)
	- [1.4. 基于 de novo assembly](#detect-method-denovo-assembly)
- [2. 专注CNV的检测](#focus-on-cnv-discovery)
	- [2.1. CNV-seq](#cnv-seq)
	- [2.2. CNV检测工具](#tools)
		- [2.2.1. GATK4 somatic CNV discovery](#gatk4)
		- [2.2.2. CNVnator](#cnvnator)
- [3. 专注SV的检测](#focus-on-sv-discovery)
	- [3.1. lumpy](#lumpy)
	- [3.2. delly](#delly)







<h1 name="title">基因组结构变异 (structural variation) 检测</h1>

<p align="center"><img src=./picture/StructralVariation-SVs-types.png width=800 /></p>

检测SVs的意义：

> - SVs对基因组的影响比起SNP更大，一旦发生往往会给生命体带来重大影响，比如导致出生缺陷、癌症等；
> - 有研究发现，基因组上的SVs比起SNP而言，更能代表人类群体的多样性特征；
> - 稀有且相同的一些结构性变异往往和疾病（包括癌症）的发生相互关联甚至还是其直接的致病诱因。

<a name="svs-detection-principle"><h2>1. 基于NGSSVs检测原理 [<sup>目录</sup>](#content)</h2></a>

目前已有的检测SVs的策略：

- Read Pair，一般称为Pair-End Mapping，简称RP或者PEM；
- Split Read，简称SR；
- Read Depth，简称RD，也有人将其称为RC——Read Count的意思，它与Read Depth是同一回事，顾名思义都是利用read覆盖情况来检测变异的方法；
- 序列从头组装（de novo Assembly， 简称AS）的方法。

<p align="center"><img src=./picture/StructralVariation-SVs-detection-principle-outline.png width=800 /></p>

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

<a name="detect-method-read-depth"><h3>1.3. Read Depth方法 [<sup>目录</sup>](#content)</h3></a>

目前有两种利用Read depth信息检测CNV的策略：

- 一种是，通过检测样本在参考基因组上read的深度分布情况来发现CNV，这类适用于单样本，也是用的比较多的一个方法；

- 另一种则是通过识别并比较两个样本在基因组上存在丢失和重复倍增的区域，以此来获得彼此相对的CNV，适用于case-control模型的样本，或者肿瘤样本Somatic CNV的发现，这有点像CGH芯片的方法，因此又被称作CNV-seq，具体信息见 [2.1. CNV-seq](#cnv-seq)。

CNVnator使用的是第一种策略，同时也广泛地被用于检测大的CNV，当然还有很多冷门的软件，这里就不再列举了；CNV-seq使用的则是第二种策略。

<a name="detect-method-denovo-assembly"><h3>1.4. 基于 de novo assembly [<sup>目录</sup>](#content)</h3></a>

其实从上面看下来，**SVs检测最大的难点实际上是read太短导致的**。就因为read太短，我们不能够在比对的时候横跨基因组重复区域；就因为read太短，很多大的Insertion序列根本就没能够看到信息；就因为read太短，比对才那么纠结，我们才需要用各种数学模型来 猜测这个变异到底应该是什么等等。

但从理论上来讲，三代测序和de novo assembly 的方法应该要算是基因组结构性变异检测上最有效的方法，它们都能够检测所有类型的结构性变异。

<p align="center"><img src=./picture/StructralVariation-SVs-detection-principle-denovo-assembly-1.png width=800 /></p>

<a name="focus-on-cnv-discovery"><h2>2. 专注CNV的检测 [<sup>目录</sup>](#content)</h2></a>

<a name="cnv-seq"><h3>2.1. CNV-seq [<sup>目录</sup>](#content)</h3></a>

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

<a name="tools"><h3>2.2. CNV检测工具 [<sup>目录</sup>](#content)</h3></a>

<a name="gatk4"><h4>2.2.1. GATK4 somatic CNV discovery [<sup>目录</sup>](#content)</h4></a>

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

<a name="cnvnator"><h4>2.2.2. CNVnator [<sup>目录</sup>](#content)</h4></a>

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

---

参考资料：

(1) [【简书】一篇文章说清楚基因组结构性变异检测的方法](https://www.jianshu.com/p/4c8e109f0e6a)

(2) [【基因学苑】好东西应该再次分享——CNV数据分析](https://mp.weixin.qq.com/s/5xuqHTjYLDd-_1jzLXEeOw)

(3) Zhou B, Ho SS, Zhang X, et al.Whole-genome sequencing analysis of CNV using low-coverage and paired-end strategies is efficient and outperforms array-based CNV analysis Journal of Medical Genetics 2018;55:735-743.

(4) [GATK4 best practice for somatic CNV discovery](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11147)

(5) [【github】GATK4官方somatic CNV鉴定流程](cnv_somatic_panel_workflow.wdl
)

(6) [【生信技能树】GATK4的CNV流程-hg38](https://mp.weixin.qq.com/s/Lvfy7Y352WhLMuzawMvroA)

(7) [【基因学苑】一个人全基因组完整数据分析脚本](https://mp.weixin.qq.com/s/TBOhU_4d3iPQIBDRfv3e0A)

