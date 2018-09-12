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

	`-R` 选项为必须选项，用于定义头文件

	> ```
	> The file must have a proper bam header with read groups. Each read group must contain the platform (PL) and sample (SM) tags. 
	> For the platform value, we currently support 454, LS454, Illumina, Solid, ABI_Solid, and CG (all case-insensitive)
	> 
	> The GATK requires several read group fields to be present in input files and will fail with errors if this requirement is not satisfied
	> 
	> ID：输入reads集的ID号； LB： reads集的文库名； SM：样本名称； PL：测序平台
	> ```


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



---

参考资料：

(1) [生信菜鸟团：GATK使用注意事项](http://www.bio-info-trainee.com/838.html)

(2) 小天师兄《全外显子组测序分析》

(3) [RNA-Seq是否可以替代WES完成外显子的变异检测?二代测序的四种Read重复是如何产生的?](https://mp.weixin.qq.com/s/RfEt-O-R2njje5Asu3WpwQ)

(4) [生信菜鸟团：仔细探究samtools的rmdup是如何行使去除PCR重复reads功能的](http://www.bio-info-trainee.com/2003.html)

(5) [生信技能树论坛：GATK之BaseRecalibrator](http://www.biotrainee.com/thread-1402-1-1.html)
