<a name="content">目录</a>

[snp-calling](#title)
- [GATK4流程](#gatk4)
	- [1. 准备配套数据](#gatk4-prepare-necessary-datasets)
	- [2. BWA: Map to Reference](#gatk4-map-ref)
	- [3. 前期处理](#gatk4-post-alignment-processing)
		- [3.1. 去除PCR重复](#gatk4-remove-read-duplicates)
			- [3.1.1. duplicates的产生原因](#reason-of-duplicates)
			- [3.1.2. PCR bias的影响](#influence-of-pcr-bias)
			- [3.1.3. 探究samtools和picard去除read duplicates的方法](#principle-of-remove-duplicates)
			- [3.1.4. 操作：排序及标记重复](#operate-remove-read-duplicates)
		- [3.2. 质量值校正](#gatk4-recallbrate-base-quality-scores)
	- [4. SNP、 INDEL位点识别与过滤](#gatk4-snp-indel-identify-and-filter)
		- [4.1. SNP calling 策略的选择](#gatk4-choice-for-snp-calling-strategies)
		- [4.2. Germline SNPs + Indels](#gatk4-germline-snps-indels)
			- [4.2.1. SNP、 INDEL位点识别](#gatk4-germline-snps-indels-identify)
			- [4.2.2. SNP、 INDEL位点过滤](#gatk4-germline-snps-indels-filter)






<h1 name="title">snp-calling</h1>

<a name="gatk4"><h1 align="center">GATK4流程</h1></p>

<p align="center"><img src=./picture/GATK4-pipeline.png width=800 /></p>

<a name="gatk4-prepare-necessary-datasets"><h2>1. 准备配套数据 [<sup>目录</sup>](#content)</h2></p>

要明确你的参考基因组版本了！！！ b36/b37/hg18/hg19/hg38，记住**b37和hg19并不是完全一样**的，有些微区别哦！！！

**1、下载hg19**

这个下载地址非常多，常用的就是NCBI，ensembl和UCSC了，但是这里推荐用这个脚本下载（下载源为UCSC）：

```
# 一个个地下载hg19的染色体
for i in $(seq 1 22) X Y M;
do echo $i;
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr${i}.fa.gz;
done

gunzip *.gz

# 用cat按照染色体的顺序拼接起来，因为GATK后面的一些步骤对染色体顺序要求非常变态，如果下载整个hg19，很难保证染色体顺序是1-22，X,Y,M
for i in $(seq 1 22) X Y M;
do cat chr${i}.fa >> hg19.fasta;
done

rm -fr chr*.fasta
```

<a name="gatk4-map-ref"><h2>2. BWA: Map to Reference [<sup>目录</sup>](#content)</h2></p>

1. 建立参考序列索引

	```
	$ bwa index -a bwtsw ref.fa
	```

	参数`-a`用于指定建立索引的算法：
	
	> - bwtsw 适用于>10M
	> - is 适用于参考序列<2G (默认-a is)

	可以不指定`-a`参数，bwa index会根据基因组大小来自动选择合适的索引方法

2. 序列比对

	```
	$ bwa mem ref.fa sample_1.fq sample_2.fq -R '@RG\tID:sample\tLB:sample\tSM:sample\tPL:ILLUMINA' \
		2>sample_map.log | samtools sort -@ 20 -O bam -o sample.sorted.bam 1>sample_sort.log 2>&1
	```

	`-R` 选项为必须选项，用于定义头文件中的SAM/BAM文件中的read group和sample信息

	> ```
	> The file must have a proper bam header with read groups. Each read group must contain the platform (PL) and sample (SM) tags. 
	> For the platform value, we currently support 454, LS454, Illumina, Solid, ABI_Solid, and CG (all case-insensitive)
	> 
	> The GATK requires several read group fields to be present in input files and will fail with errors if this requirement is not satisfied
	> 
	> ID：输入reads集的ID号； LB： reads集的文库名； SM：样本名称； PL：测序平台
	> ```

	查看SAM/BAM文件中的read group和sample信息：
	
	```
	$ samtools view -H /path/to/my.bam | grep '^@RG'
	@RG ID:0    PL:solid    PU:Solid0044_20080829_1_Pilot1_Ceph_12414_B_lib_1_2Kb_MP_Pilot1_Ceph_12414_B_lib_1_2Kb_MP   LB:Lib1 PI:2750 DT:2008-08-28T20:00:00-0400 SM:NA12414  CN:bcm
	@RG ID:1    PL:solid    PU:0083_BCM_20080719_1_Pilot1_Ceph_12414_B_lib_1_2Kb_MP_Pilot1_Ceph_12414_B_lib_1_2Kb_MP    LB:Lib1 PI:2750 DT:2008-07-18T20:00:00-0400 SM:NA12414  CN:bcm
	@RG ID:2    PL:LS454    PU:R_2008_10_02_06_06_12_FLX01080312_retry  LB:HL#01_NA11881    PI:0    SM:NA11881  CN:454MSC
	@RG ID:3    PL:LS454    PU:R_2008_10_02_06_07_08_rig19_retry    LB:HL#01_NA11881    PI:0    SM:NA11881  CN:454MSC
	@RG ID:4    PL:LS454    PU:R_2008_10_02_17_50_32_FLX03080339_retry  LB:HL#01_NA11881    PI:0    SM:NA11881  CN:454MSC
	...
	```

	read group信息不仅会加在头信息部分，也会在比对结果的每条记录里添加一个 `RG:Z:*` 标签

	```
	$ samtools view /path/to/my.bam | grep '^@RG'
	EAS139_44:2:61:681:18781    35  1   1   0   51M =   9   59  TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA B<>;==?=?<==?=?=>>?>><=<?=?8<=?>?<:=?>?<==?=>:;<?:= RG:Z:4  MF:i:18 Aq:i:0  NM:i:0  UQ:i:0  H0:i:85 H1:i:31
	EAS139_44:7:84:1300:7601    35  1   1   0   51M =   12  62  TAACCCTAAGCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA G<>;==?=?&=>?=?<==?>?<>>?=?<==?>?<==?>?1==@>?;<=><; RG:Z:3  MF:i:18 Aq:i:0  NM:i:1  UQ:i:5  H0:i:0  H1:i:85
	EAS139_44:8:59:118:13881    35  1   1   0   51M =   2   52  TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA @<>;<=?=?==>?>?<==?=><=>?-?;=>?:><==?7?;<>?5?<<=>:; RG:Z:1  MF:i:18 Aq:i:0  NM:i:0  UQ:i:0  H0:i:85 H1:i:31
	EAS139_46:3:75:1326:2391    35  1   1   0   51M =   12  62  TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA @<>==>?>@???B>A>?>A?A>??A?@>?@A?@;??A>@7>?>>@:>=@;@ RG:Z:0  MF:i:18 Aq:i:0  NM:i:0  UQ:i:0  H0:i:85 H1:i:31
	...
	```

	若原始SAM/BAM文件没有read group和sample信息，可以通过AddOrReplaceReadGroups添加这部分信息

	```
	$ java -jar picard.jar AddOrReplaceReadGroups \
	      I=input.bam \
	      O=output.bam \
	      RGID=4 \
	      RGLB=lib1 \
	      RGPL=illumina \
	      RGPU=unit1 \
	      RGSM=20
	```

<a name="gatk4-post-alignment-processing"><h2>3. 前期处理 [<sup>目录</sup>](#content)</h2></p>

在进行本部分的操作之前先要做好以下两部的准备工作

1、创建GATK索引。用Samtools为参考序列创建一个索引，这是为方便GATK能够快速地获取fasta上的任何序列做准备

```
$ samtools faidx database/chr17.fa	# 该命令会在chr17.fa所在目录下创建一个chr17.fai索引文件
```

2、生成.dict文件

```
$ gatk CreateSequenceDictionary -R database/chr17.fa -O database/chr17.dict
```

GATK4的调用语法：

```
gatk [--java-options "-Xmx4G"] ToolName [GATK args]
```
<a name="gatk4-remove-read-duplicates"><h3>3.1. 去除PCR重复 [<sup>目录</sup>](#content)</h3></p>

<a name="reason-of-duplicates"><h4>3.1.1. duplicates的产生原因 [<sup>目录</sup>](#content)</h4></p>

<p align="center"><img src=./picture/GATK4-pipeline-remove-duplicates-reason-of-duplicates.jpg width=900/></p>

- **PCR duplicates（PCR重复）**

PCR扩增时，同一个DNA片段会产生多个相同的拷贝，第4步测序的时候，这些来源于同！一！个！拷贝的DNA片段会结合到Fellowcell的不同位置上，生成完全相同的测序cluster，然后被测序出来，这些相同的序列就是duplicate

- **Cluster duplicates**

生成测序cluster的时候，某一个cluster中的DNA序列可能搭到旁边的另一个cluster的生成位点上，又再重新长成一个相同的cluster，这也是序列duplicate的另一个来源，这个现象在Illumina HiSeq4000之后的Flowcell中会有这类Cluster duplicates

- **Optical duplicates（光学重复）**

某些cluster在测序的时候，捕获的荧光亮点由于光波的衍射，导致形状出现重影（如同近视散光一样），导致它可能会被当成两个荧光点来处理。这也会被读出为两条完全相同的reads

- **Sister duplicates**

它是文库分子的两条互补链同时都与Flowcell上的引物结合分别形成了各自的cluster被测序，最后产生的这对reads是完全反向互补的。比对到参考基因组时，也分别在正负链的相同位置上，在有些分析中也会被认为是一种duplicates。

<a name="influence-of-pcr-bias"><h4>3.1.2. PCR bias的影响 [<sup>目录</sup>](#content)</h4></p>

1. DNA在打断的那一步会发生一些损失， 主要表现是会引发一些碱基发生颠换变换（嘌呤-变嘧啶或者嘧啶变嘌呤） ， 带来假的变异。 PCR过程会扩大这个信号， 导致最后的检测结果中混入了假的结果；

2. PCR反应过程中也会带来新的碱基错误。 发生在前几轮的PCR扩增发生的错误会在后续的PCR过程中扩大， 同样带来假的变异；

3. 对于真实的变异， PCR反应可能会对包含某一个碱基的DNA模版扩增更加剧烈（这个现象称为PCR Bias） 。 因此， 如果反应体系是对含有reference allele的模板扩增偏向强烈， 那么变异碱基的信息会变小， 从而会导致假阴。

<p align="center"><img src=./picture/GATK4-pipeline-remove-duplicates-1.png width=900/></p>

<a name="principle-of-remove-duplicates"><h4>3.1.3. 探究samtools和picard去除read duplicates的方法 [<sup>目录</sup>](#content)</h4></p>

**1、samtools**

samtools 去除 duplicates 使用 **rmdup**

```
$ samtools rmdup [-sS] <input.srt.bam> <out.bam>
```

只需要开始`-s`的标签， 就可以对单端测序进行去除PCR重复。其实对单端测序去除PCR重复很简单的，因为比对flag情况只有0,4,16，只需要它们比对到染色体的起始终止坐标一致即可，flag很容易一致。

> Remove potential PCR duplicates: if multiple read pairs have identical external coordinates, only retain the pair with highest mapping quality. I

<p align="center"><img src=./picture/GATK4-pipeline-remove-duplicates-samtools-principle-1.png width=900/></p>

但是对于双端测序就有点复杂了~

> In the paired-end mode, this command ONLY works with FR orientation and requires ISIZE is correctly set. It does not work for unpaired reads (e.g. two ends mapped to different chromosomes or orphan reads).

<p align="center"><img src=./picture/GATK4-pipeline-remove-duplicates-samtools-principle-2.png width=900/></p>

很明显可以看出，去除PCR重复不仅仅需要它们比对到染色体的起始终止坐标一致，尤其是**flag**，在双端测序里面一大堆的flag情况，所以我们的94741坐标的5个reads，一个都没有去除！

这样的话，双端测序数据，用samtools rmdup效果就很差，所以很多人建议用picard工具的MarkDuplicates 功能

**2、picard**

picard对于单端或者双端测序数据并没有区分参数，可以用同一个命令！

对应单端测序，picard的处理结果与samtools rmdup没有差别，不过这个java软件的缺点就是**奇慢无比**




<a name="operate-remove-read-duplicates"><h4>3.1.4. 操作：排序及标记重复 [<sup>目录</sup>](#content)</h4></p>

<p align="center"><img src=./picture/GATK4-pipeline-remove-duplicates-2.png width=800/></p>

<p align="center"><img src=./picture/GATK4-pipeline-remove-duplicates-3.png width=800/></p>

1、排序（SortSam）

- 对sam文件进行排序并生成bam文件，将sam文件中同一染色体对应的条目按照坐标顺序从小到大进行排序
- GATK4的排序功能是通过`picard SortSam`工具实现的。虽然`samtools sort`工具也可以实现该功能，但是在GATK流程中还是推荐用picard实现，因为SortSam会在输出文件的头信息部分添加一个SO标签用于说明文件已经被成功排序，且**这个标签是必须的**，GATK需要检查这个标签以保证后续分析可以正常进行
- `https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_SortSam.php`

```
# 使用GATK命令
$ gatk SortSam -I mapping/T.chr17.sam -O preprocess/T.chr17.sort.bam -R database/chr17.fa -SO coordinate --CREATE_INDEX
# 使用picard命令
$ java -jar picard.jar SortSam \
      I=input.bam \
      O=sorted.bam \
      SORT_ORDER=coordinate
```

<p align="center"><img src=./picture/GATK4-pipeline-remove-duplicates-4.png width=900/></p>

如何检查是否成功排序？

```
$ samtools view -H /path/to/my.bam
@HD     VN:1.0  GO:none SO:coordinate
@SQ     SN:1    LN:247249719
@SQ     SN:2    LN:242951149
@SQ     SN:3    LN:199501827
@SQ     SN:4    LN:191273063
@SQ     SN:5    LN:180857866
@SQ     SN:6    LN:170899992
@SQ     SN:7    LN:158821424
@SQ     SN:8    LN:146274826
@SQ     SN:9    LN:140273252
@SQ     SN:10   LN:135374737
@SQ     SN:11   LN:134452384
@SQ     SN:12   LN:132349534
@SQ     SN:13   LN:114142980
@SQ     SN:14   LN:106368585
@SQ     SN:15   LN:100338915
@SQ     SN:16   LN:88827254
@SQ     SN:17   LN:78774742
@SQ     SN:18   LN:76117153
@SQ     SN:19   LN:63811651
@SQ     SN:20   LN:62435964
@SQ     SN:21   LN:46944323
@SQ     SN:22   LN:49691432
@SQ     SN:X    LN:154913754
@SQ     SN:Y    LN:57772954
@SQ     SN:MT   LN:16571
@SQ     SN:NT_113887    LN:3994
...
```

若随后的比对记录中的contig那一列的顺序与头文件的顺序一致，且在头信息中包含`SO:coordinate`这个标签，则说明，该文件是排序过的

2、标记重复（Markduplicates）

- 标记文库中的重复
- `https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_markduplicates_MarkDuplicates.php`


```
gatk MarkDuplicates -I preprocess/T.chr17.sort.bam -O preprocess/T.chr17.markdup.bam -M preprocess/T.chr17.metrics --CREATE_INDEX
```

<p align="center"><img src=./picture/GATK4-pipeline-remove-duplicates-5.png width=900/></p>

<a name="gatk4-recallbrate-base-quality-scores"><h3>3.2. 质量值校正 [<sup>目录</sup>](#content)</h3></p>

检测碱基质量分数中的系统错误，需要用到 GATK4 中的 BaseRecalibrator 工具

**碱基质量分数重校准（Base quality score recalibration，BQSR)**，就是利用机器学习的方式调整原始碱基的质量分数。它分为两个步骤:

> - 利用已有的snp数据库，建立相关性模型，产生重校准表( recalibration table)
> - 根据这个模型对原始碱基进行调整，只会调整非已知SNP区域。

注：如果不是人类基因组，并且也缺少相应的已知SNP数据库，可以通过严格SNP筛选过程（例如结合GATK和samtools）建立一个snp数据库。

1、建立较正模型

质量值校正，这一步需要用到variants的known-sites，所以需要先准备好已知的snp，indel的VCF文件：

```
# 下载known-site的VCF文件，到Ensembl上下载
$ wget -c -P Ref/mouse/mm10/vcf ftp://ftp.ensembl.org/pub/release-93/variation/vcf/mus_musculus/mus_musculus.vcf.gz >download.log &
$ cd Ref/mouse/mm10/vcf && gunzip mus_musculus.vcf.gz && mv mus_musculus.vcf dbsnp_150.mm10.vcf
# 建好vcf文件的索引，需要用到GATK工具集中的IndexFeatureFile，该命令会在指定的vcf文件的相同路径下生成一个以".idx"为后缀的文件
$ gatk IndexFeatureFile -F dbsnp_150.mm10.vcf

# 建立较正模型
$ gatk BaseRecalibrator -R Ref/mouse/mm10/bwa/mm10.fa -I PharmacogenomicsDB/mouse/SAM/ERR118300.enriched.markdup.bam -O \
PharmacogenomicsDB/mouse/SAM/ERR118300.recal.table --known-sites Ref/mouse/mm10/vcf/dbsnp_150.mm10.vcf
```


2、质量值校准

```
# 质量校正
$ gatk ApplyBQSR -R Ref/mouse/mm10/bwa/mm10.fa -I PharmacogenomicsDB/mouse/SAM/ERR118300.enriched.markdup.bam -bqsr \
PharmacogenomicsDB/mouse/SAM/ERR118300.recal.table -O PharmacogenomicsDB/mouse/SAM/ERR118300.recal.bam
```

<a name="gatk4-snp-indel-identify-and-filter"><h2>4. SNP、 INDEL位点识别 [<sup>目录</sup>](#content)</h2></p>

<a name="gatk4-choice-for-snp-calling-strategies"><h3>4.1. SNP calling 策略的选择 [<sup>目录</sup>](#content)</h3></p>

当你有多个samples，然后你call snp时候，你是应该将所有sample分开call完之后，再merge在一起。还是直接将所有samples，同时用作input然后call snp呢？

这里需要知道有哪些snp calling的策略：

- **single sample calling**：每一个sample的bam file都进行单独的snp calling，然后每个sample单独snp calling结果再合成一个总的snp calling的结果。

- **batch calling**： 一定数目群集的bamfiles 一起calling snps，然后再merge在一起

- **joint calling**： 所有samples的BAM files一起call 出一个包含所有samples 变异信息的output

一般来说，如果条件允许（computational power等），使用joint calling ，即将所有samples同时call是比较优的选择

原因：

**1、对于低频率的变异具有更高更好的检测sensitivity**

在joint calling中，由于所有samples中所有的位点都是同时call的，换句话说就是，所有位点的信息都是share的，因此可以如果某些samples中个别位点是低频率的，但是可能在其他samples中，该位点的频率比较高，因此可以准确的对低频位点有更加好的calling 效果。

<p align="center"><img src=./picture/GATK4-pipeline-snp-calling-strategies.jpg width=800 /></p>

如左图，在1到n的samples中，碱基G只出现在其中两个samples中。如果我们将这些samples单独的call snps，这个低频率G的位点将会被忽略。但是joint calling 却可以允许将这些碱基G出现的频率进行累加，将该低频率的突变也call出来。

在右图中，上面的sample是一个跟ref一样具有纯合位点，下面的sample在一些位点中有部分数据缺少的情况，例如rs429358.如果将这两个samples分开call snp，然后merge一起，这样会错误地将这两个samples，在这个位点上看作具有类似的变异频率，但是其实由于下面的samples在改位点的区域存在数据缺少，只能看作non-informative。

**2、更好的过滤掉假阳性的变异callling的结果**

现在使用过滤变异的方法例如VQSR等利用的统计模型，都基于一个比较大的samples size。joint calling 这种方法可以提供足够的数据，确保过滤这一步是统一应用于所有samples的。


<a name="gatk4-germline-snps-indels"><h3>4.2. Germline SNPs + Indels [<sup>目录</sup>](#content)</h3></p>

将一个或多个个体放在一起call snp，得到一个 joint callset，该snp calling的策略称为**joint calling**

<p align="center"><img src=./picture/GATK4-pipeline-Germline-SNPs-Indels.png width=800 /></p>

- 第一步，单独为每个样本生成后续分析所需的中间文件——gVCF文件。这一步中包含了对原始fastq数据的质控、比对、排序、标记重复序列、BQSR和HaplotypeCaller gVCF等过程。这些过程全部都适合在单样本维度下独立完成。值得注意的是，与单样本模式不同，该模式中每个样本的gVCF应该成为这类流程的标配，在后续的步骤中我们可以通过gVCF很方便地完成群体的Joint Calling；

- 第二步，依据第一步完成的gVCF对这个群体进行Joint Calling，从而得到这个群体的变异结果和每个人准确的基因型（Genotype），最后使用VQSR完成变异的质控。

<a name="gatk4-germline-snps-indels-identify"><h4>4.2.1. SNP、 INDEL位点识别 [<sup>目录</sup>](#content)</h4></p>

1、单独为每个样本生成后续分析所需的中间文件——gVCF文件

```
$ gatk HaplotypeCaller -R Ref/chr17.fa -I sam/T.chr17.recal.bam -ERC GVCF --dbsnp ../Ref/VCF/dbsnp_138.hg19.vcf \
	-O calling/T.chr17.raw.snps.indels.vcf -L chr17:7400000-7800000
```

HaplotypeCaller可以同时call snp 和indel，通过局部 de-novo 组装而非基于mapping的结果。一旦程序在某个区域发现了变异信号，它会忽略已有的mapping信息，在该区域执行reads重组装

如何对指定的区域call snp？

> 使用 `-L` 参数指定识别突变位点的区域
> 
> 1. 提供字符串指定单个interval，如`-L chr17:7400000-7800000` 只识别17号染色体7400000-7800000 区域的突变位点。或者使用 `-XL` 参数排除某个区域
> 
> 2. 提供一个保存interval lists的文件。GATK支持多种格式的interval lists文件
> 
> 	- Picard-style `.interval_list`
> 
> 		由 SAM-like header 和 intervals 信息组成，坐标系统为 1-based
> 
> 		```
> 		@HD     VN:1.0  SO:coordinate
> 		@SQ     SN:1    LN:249250621    AS:GRCh37       UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta   M5:1b22b98cdeb4a9304cb5d48026a85128     SP:Homo Sapiens
> 		@SQ     SN:2    LN:243199373    AS:GRCh37       UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta   M5:a0d9851da00400dec1098a9255ac712e     SP:Homo Sapiens
> 		1       30366   30503   +       target_1
> 		1       69089   70010   +       target_2
> 		1       367657  368599  +       target_3
> 		1       621094  622036  +       target_4
> 		1       861320  861395  +       target_5
> 		1       865533  865718  +       target_6
> 		```
> 
> 	- GATK-style `.list` or `.intervals`
> 
> 		这种格式很简单，intervals需要写成这种格式：`<chr>:<start>-<stop>`，坐标系统为 1-based
> 
> 	- BED files `.bed`
> 
> 		BED3格式：`<chr> <start> <stop>`，坐标系统为 0-based，GATK只接受 1-based 坐标系统，因此GATK会根据文件后缀 `.bed` 识别bed文件格式，然后会将 0-based 转换为 1-based
> 
> 	**注意：**
> 	> - intervals 必须按照 reference 的坐标进行排序
> 	> - 可以提供多个intervals集合，文件格式不必一致：`-L intervals1.interval_list -L intervals2.list -L intervals3.bed ...`

2、合并多个GVCF文件得到GenomicsDB，为joint genotyping做准备。若是单样本跳过该步骤

```
$ gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
      -V data/gvcfs/mother.g.vcf.gz \
      -V data/gvcfs/father.g.vcf.gz \
      -V data/gvcfs/son.g.vcf.gz \
      --genomicsdb-workspace-path my_database \
      -L 20
```

也可以使用CombineGVCFs实现该功能。CombineGVCFs继承自GATK4之前的版本，但是它的运行速度不如GenomicsDB

```
$ gatk CombineGVCFs \
   -R reference.fasta \
   --variant sample1.g.vcf.gz \
   --variant sample2.g.vcf.gz \
   -O cohort.g.vcf.gz
```

3、联合call snp，得到 joint-called SNP&indel

```
# 单样本
$ gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R Homo_sapiens_assembly38.fasta \
   -V input.g.vcf.gz \
   -O output.vcf.gz

# 多样本
$ gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R Homo_sapiens_assembly38.fasta \
   -V gendb://my_database \
   -O output.vcf.gz
```

<a name="gatk4-germline-snps-indels-filter"><h4>4.2.2. SNP、 INDEL位点过滤[<sup>目录</sup>](#content)</h4></p>

- 方法一：通过质量校正来过滤（Filter Variants by Variant (Quality Score) Recalibration）

	1、用机器学习的方法基于已知的变异位点对caller给出的原始 variant quality score 进行校正  (VQSR)，包含两步：（1）构建校正模型；（2）应用模型进行质量值校正

	由于评价SNP和Indel质量高低的标准不同，因此需要分SNP和Indel两种不同的模式，分别进行校正
	
	```
	# SNP mode
	## 构建校正模型
	$ gatk VariantRecalibrator \
	   -R Homo_sapiens_assembly38.fasta \
	   -V input.vcf.gz \
	   --resource hapmap,known=false,training=true,truth=true,prior=15.0:hapmap_3.3.hg38.sites.vcf.gz \
	   --resource omni,known=false,training=true,truth=false,prior=12.0:1000G_omni2.5.hg38.sites.vcf.gz \
	   --resource 1000G,known=false,training=true,truth=false,prior=10.0:1000G_phase1.snps.high_confidence.hg38.vcf.gz \
	   --resource dbsnp,known=true,training=false,truth=false,prior=2.0:Homo_sapiens_assembly38.dbsnp138.vcf.gz \
	   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
	   -mode SNP \
	   -O output.snps.recal \
	   --tranches-file output.snps.tranches \
	   --rscript-file output.snps.plots.R
	## 应用模型进行质量值校正
	$ gatk ApplyVQSR \
	   -R Homo_sapiens_assembly38.fasta \
	   -V input.vcf.gz \
	   -O output.snps.VQSR.vcf.gz \
	   --truth-sensitivity-filter-level 99.0 \
	   --tranches-file output.snps.tranches \
	   --recal-file output.snps.recal \
	   -mode SNP

	# Indel mode
	## 上一步SNP mode产生的输出作为这一步的输入
	## 构建校正模型
	$ gatk VariantRecalibrator \
	   -R Homo_sapiens_assembly38.fasta \
	   -V input.snps.VQSR.vcf.gz \
	   --resource hapmap,known=false,training=true,truth=true,prior=15.0:hapmap_3.3.hg38.sites.vcf.gz \
	   --resource omni,known=false,training=true,truth=false,prior=12.0:1000G_omni2.5.hg38.sites.vcf.gz \
	   --resource 1000G,known=false,training=true,truth=false,prior=10.0:1000G_phase1.snps.high_confidence.hg38.vcf.gz \
	   --resource dbsnp,known=true,training=false,truth=false,prior=2.0:Homo_sapiens_assembly38.dbsnp138.vcf.gz \
	   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
	   -mode INDEL \
	   -O output.indels.recal \
	   --tranches-file output.indels.tranches \
	   --rscript-file output.indels.plots.R
	## 应用模型进行质量值校正
	$ gatk ApplyVQSR \
	   -R Homo_sapiens_assembly38.fasta \
	   -V input.snps.VQSR.vcf.gz \
	   -O output.snps.indels.VQSR.vcf.gz \
	   --truth-sensitivity-filter-level 99.0 \
	   --tranches-file output.indels.tranches \
	   --recal-file output.indels.recal \
	   -mode INDEL
	```

	执行变异校正需要满足两个条件：

	> - 大量高质量的已知variants作为训练集，而这对于许多的物种是不满足的
	> - 数据集不能太小，因为它需要足够的数据集来识别 good vs. bad variants
	
	因此对于只有一两个样本、靶向测序、RNA-seq、非模式动物的数据，都不推荐进行变异质量值校正。若以上两个条件都不满足，则需要采取直接过滤 (hard-filtering) 的策略

- 方法二：直接过滤 (hard-filtering)

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

---

参考资料：

(1) [生信菜鸟团：GATK使用注意事项](http://www.bio-info-trainee.com/838.html)

(2) 小天师兄《全外显子组测序分析》

(3) [RNA-Seq是否可以替代WES完成外显子的变异检测?二代测序的四种Read重复是如何产生的?](https://mp.weixin.qq.com/s/RfEt-O-R2njje5Asu3WpwQ)

(4) [生信菜鸟团：仔细探究samtools的rmdup是如何行使去除PCR重复reads功能的](http://www.bio-info-trainee.com/2003.html)

(5) [生信技能树论坛：GATK之BaseRecalibrator](http://www.biotrainee.com/thread-1402-1-1.html)

(6) [GATK Forum: Collected FAQs about input files for sequence read data (BAM/CRAM)](https://gatkforums.broadinstitute.org/gatk/discussion/1317/collected-faqs-about-bam-files)

(7) [GATK Forum: What input files does the GATK accept / require?](https://software.broadinstitute.org/gatk/documentation/article.php?id=1204)

(8) [GATK Forum: Collected FAQs about interval lists](https://software.broadinstitute.org/gatk/documentation/article.php?id=1319)

(9) [GATK Forum：GATK Best Practices](https://software.broadinstitute.org/gatk/best-practices/)

(10) [公众号碱基矿工：GATK4全基因组数据分析最佳实践 ，我以这篇文章为标志，终结当前WGS系列数据分析的流程主体问题 | 完全代码](https://mp.weixin.qq.com/s/Sa019WuSg8fRQgkWAIG4pQ)

(11) [生信菜鸟团：生信笔记：call snp是应该一起call还是分开call？](https://mp.weixin.qq.com/s/XVkFuU2zWLY5r6grDwzkcA)
