<a name="content">目录</a>

[基因型填充](#title)
- [问题描述](#description)
- [技术来源的基因型缺失](#genotype-missing-drived-from-tech)
- [缺失的判断：缺失率](#missing-rate)
- [基因型缺失的影响](#effect)
- [基因型填充的原理](#principle-of-imputation)
- [实现工具](#tools)
	- [IMPUTE2](#tools-impute2)
		- [6.1.1. 两种应用场景](#tools-impute2-2-scenario)
		- [6.1.2. Best Practices](#tools-impute2-best-practices)
		- [6.1.3. 1000 Genomes Imputation Cookbook](#tools-impute2-1kgp-cookbook)
			- [6.1.3.1. Before Imputation](#tools-impute2-1kgp-cookbook-before-imputation)
			- [6.1.3.2. Pre-Phasing](#tools-impute2-1kgp-cookbook-pre-phasing)
				- [6.1.3.2.1. using IMPUTE2](#tools-impute2-1kgp-cookbook-pre-phasing-using-impute2)
				- [6.1.3.2.2. using SHAPEIT (recommended)](#tools-impute2-1kgp-cookbook-pre-phasing-using-shapeit)
			- [6.1.3.3. Imputation](#tools-impute2-1kgp-cookbook-imputation)



<h1 name="title">基因型填充</h1>

<a name="description"><h2>问题描述 [<sup>目录</sup>](#content)</h2></a>

基因型缺失：样本中没有被测序数据覆盖到的区域，基因型就属于未知的，我们将之称为缺失位点

<p align="center"><img src=./picture/Genotype-Imputation-description.png width=800 /></p>

基因型数据的缺失又分为**遗传性缺失**和**检测性缺失**：

> - 遗传性缺失：个体遗传信息的变异（例如，这个位点DNA片段真实缺失）导致的基因型缺失
> 
> - 检测性缺失：由于检测技术的局限、错误等导致的信息丢失。各类基因型检测技术都会产生检测性的基因型缺失。

<a name="genotype-missing-drived-from-tech"><h2>技术来源的基因型缺失 [<sup>目录</sup>](#content)</h2></a>

1. 全基因组重测序技术

	全基因组重测序理论上应该覆盖整个基因组，因此**未覆盖的区域**都可以被定义为缺失。那么群体研究中的低深度测序（一般平均深度低于10X），不可避免会产生大量随机缺失。

2. 简化基因组测序

	简化基因组测序是通过酶切，并富集限制性内切酶周边的片段并进行测序的策略。针对简化基因组，我们称的缺失一般指的是没有被检测到的酶切片段相关的位点。简化基因组的缺失，主要与**酶切效率**有关。酶切效率越高，缺失率越低。

3. 外显子测序以及目标区域捕获测序

	同简化基因组测序类似，基于探针杂交的DNA捕获以及测序技术，同样会产生大量的缺失。这种缺失主要是由于**探针杂交捕获的效率**所致。

4. SNP芯片

	SNP芯片利用芯片杂交后的荧光信号，来判断某个位点的基因型。SNP芯片同样也会产生大量缺失。但在实际的研究中，SNP 芯片主要面临的问题是**芯片型号不同**，甚至来源不同的厂商，那么芯片中包含的SNP位点也不同。当来源不同的数据一起分析的时候，将面临数据不一致的问题。简单说来，就是你有的我没有，我有的你没有。如下图，Affymetrix和illuminate两大SNP 芯片厂商生产的人类芯片就使用的是不同的SNP集，当放在一起分析的时候将面临SNP不一致的问题。

	<p align="center"><img src=./picture/Genotype-Imputation-missing-resource-snpArray.jpg width=800 /></p>

注意！

>** 基因型缺失是一个相对性的概念**。以上缺失的概念都是针对同种技术的比较。不同的技术比较，也可以定义为缺失。例如，同样一份样本，我们使用全部以上4种技术检测。如果以全基因组高深度测序（>30X）为参照标准，后续的3种技术都有大量位点没有检测到，处于基因型缺失的状态。

<a name="missing-rate"><h2>缺失的判断：缺失率 [<sup>目录</sup>](#content)</h2></a>

分为**样本水平的缺失率**和**位点水平的缺失率**

例如下图，0、1、2 分别代表三种检测到的基因型，图中缺失位点使用“？”表示。那么样本1的缺失率=20%（总体10个位点，有两个位点缺失），而位点2的缺失率=60%（总体5个位点，有3个位点缺失）

<p align="center"><img src=./picture/Genotype-Imputation-missing-rate.png width=800 /></p>

<a name="effect"><h2>基因型缺失的影响 [<sup>目录</sup>](#content)</h2></a>

基因型缺失最直接的影响就是这个位置的**信息缺失**，从而影响下游分析（包括遗传图谱构建，QTL定位，选择压力分析，GWAS分析等）的信息完整性和准确性。


例如，（b）中红色的点是（a）中缺失的位点。而与性状关联的SNP位点，恰恰位于虚线所在的区域内。这些显著位点在（a）中是缺失的，所以（a）没有检测到关联信号，从丢失了非常关键的信息

<p align="center"><img src=./picture/Genotype-Imputation-missing-effect.jpg width=500 /></p>

基因型缺失对GWAS分析、选择压力分析影响都比较大

<a name="principle-of-imputation"><h2>基因型填充的原理 [<sup>目录</sup>](#content)</h2></a>

原理：

> 基于**家系样本的遗传特性**。具有已知亲缘关系的个体之间具有共享的单体型（haplotype），这些由有限个遗传标记所构成的单体型随祖先一起遗传，反映连锁不平衡。
> 
> 在具有相同单体型的家系中，**遗传标记少的样本可以参照遗传标记多的样本进行基因型填充**。对于没有亲缘关系的样本，以上理论也基本适用，主要的差别在于**无血缘关系的样本之间共享的单体型比家系样本之间的要短很多**。对无亲缘关系样本进行基因型填充需要一个高密度遗传标记构成的单体型图谱作为参照。
> 
> 通过对比待填充样本和参考模板，找到两者之间共有的单体型，然后就可以将匹配上的参考模板中的位点复制到目标数据集中。

常见imputation的基本逻辑包括两步：

> - 从目标位点/区域非缺失的位点中，总结这个区域的基因型规律，并分类。其实就是分析各个区域的单体型组成；
> 
> - 根据某样本缺失位点的上下其他非缺失位点，判断这个区域属于哪种单倍型。然后根据所属单倍型的基因型补充该样本的缺失位点；

<p align="center"><img src=./picture/Genotype-Imputation-example.jpg width=800 /></p>

根据缺失样本有限的基因型信息（仅有3个位点），就可以判断这个样本与参考单倍型集中的哪种单倍型最为相似（图中分别对应紫色、绿色、黄色三种单倍型）。然后，将对应的最相似的单倍型赋予给该样本，从而让该样本获得完整的基因型，图b

<a name="tools"><h2>实现工具 [<sup>目录</sup>](#content)</h2></a>

(1) **计算密集型**，比如IMPUTE、 IMPUTE2、MACH、 和fastPHASE/BIMBAM

这种类型的方法在填充的过程中**充分考虑到全部可以观察到的基因型信息**，使得对缺失值的估算更加精确；但以上大部分软件都是针对人类的开发的。人类种群的遗传特性是个体杂合率较高、近交率低、系谱关系来源随机。很多植物，尤其作物的遗传特性则和人类相反。

(2) **计算高效型**，比如PLINK、TUNA、WHAP和BEAGLE

此种算法**仅仅关注与特定位点相邻的一小部分标记**的基因型，因此在计算上更加快捷

<a name="tools-impute2"><h3>IMPUTE2 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/Genotype-Imputation-tools-IMPUTE2-1.jpg width=800 /></p>

<a name="tools-impute2-2-scenario"><h3>两种应用场景 [<sup>目录</sup>](#content)</h3></a>

Impute2的基因填充 (genotype imputation) 分为两种应用情景：

1. ONE REFERENCE PANEL

	<p align="center"><img src=./picture/Genotype-Imputation-tools-IMPUTE2-2.jpg width=800 /></p>

	```
	./impute2 \
	 -m ./Example/example.chr22.map \
	 -h ./Example/example.chr22.1kG.haps \
	 -l ./Example/example.chr22.1kG.legend \
	 -g ./Example/example.chr22.study.gens \
	 -strand_g ./Example/example.chr22.study.strand \
	 -int 20.4e6 20.5e6 \
	 -Ne 20000 \
	 -o ./Example/example.chr22.one.phased.impute2
	```

	参数说明：

	> - `-m <file>`: 目标区域重组率图谱文件(Fine-scale recombination map for the region to be analyzed)，记录的是基因组中各个位点的重组率和彼此间物理距离的关系
	> 
	> 
	> 这个文件应该包含三列：
	> 
	> ```
	> (1) physical position: in base pairs
	> (2) recombination rate: between current position and next position in map (in cM/Mb)
	> (3) genetic map position: in cM
	> 
	> 例如：
	> position COMBINED_rate(cM/Mb) Genetic_Map(cM)
	> 35326 0.251801 0.000000
	> 35411 0.482009 0.000021
	> 40483 0.598191 0.002466
	> ```
	> 
	> ---
	> 
	> - `-h <file 1> <file 2>`: 已知的单体型信息文件，每行表示一个SNP位点，每列表示一个单体型 (one row per SNP and one column per haplotype)
	> 
	> 
	> 所有的allele必须表示成0或1的形式
	> 
	> 一旦用`-h`参数指定一个单体型文件，就需要用`-l`参数指定一个对应的Legend文件
	> 
	> Impute2允许同时指定两个单体型文件： `-h <file 1> <file 2>`
	> 
	> ---
	> 
	> - `-l <file 1> <file 2>`：与单体型文件对应的Legend文件，保存的是对每个SNP位点的描述信息
	> 
	> 
	> 这个文件包含四列：
	> 
	> ```
	> rsID, physical position (in base pairs), allele 0, and allele 1
	> 
	> 最后两列的 allele 0 和 allele 1 是对碱基组成的说明
	> ```
	> 
	> ---
	> 
	> - `-g <file>`: 包含目标研究群体的genotypes的文件，即[Genotype File Format](http://www.stats.ox.ac.uk/%7Emarchini/software/gwas/file_format.html#Genotype_File_Format)，对它进行后续的基因型填充 (impute) 和分型 (phase)
	> 
	> 
	> 该文件每行表示一个SNP，前五列分别为：
	> 
	> ```
	> (1) SNP ID：这一列一般表示为染色体号
	> (2) RS ID of the SNP
	> (3) base-pair position of the SNP
	> (4) the allele coded A
	> (5) the allele coded B
	> ```
	> 
	> 紧接着的3列是群体中的一个个体的三种可能的基因型：AA，AB或BB
	> 
	> 再接着3列是第二个个体的，以此类推，示例文件如下：
	> 
	> ```
	> SNP1 rs1 1000 A C 1 0 0 1 0 0
	> SNP2 rs2 2000 G T 1 0 0 0 1 0
	> SNP3 rs3 3000 C T 1 0 0 0 1 0
	> SNP4 rs4 4000 C T 0 1 0 0 1 0
	> SNP5 rs5 5000 A G 0 1 0 0 0 1
	> ```
	> 
	> ---
	> 
	> - `-strand_g <file>`: 指定SNP所在的链的方向
	> 
	> 该文件每行表示一个SNP，包含两列（列之间用单个空格隔开）：（1）SNP所在的碱基位置；（2）链的方向，`+`或`-`；
	> 
	> - `-int <lower> <upper>`: 用于基因型推断的基因组间隔的长度，可以以长格式表示，如 ` -int 5420000 10420000`，也可以以指数形式表示，如 `-int 5.42e6 10.42e6`
	> 
	> - `-Ne <int>`: 这个参数的说明看不懂，把原文贴在下面：
	> 
	> ```
	> "Effective size" of the population (commonly denoted as Ne in the population genetics literature) from which your dataset was
	>  sampled. This parameter scales the recombination rates that IMPUTE2 uses to guide its model of linkage disequilibrium patterns.
	> When most imputation runs were conducted with reference panels from HapMap Phase 2, we suggested values of 11418 for imputation 
	> from HapMap CEU, 17469 for YRI, and 14269 for CHB+JPT. 
	> ```
	> - `-o <file>`: 输出文件名，文件格式与`-g`参数指定的文件相同，即都是Genotype File Format

2. TWO REFERENCE PANELS

	<p align="center"><img src=./picture/Genotype-Imputation-tools-IMPUTE2-3.jpg width=800 /></p>

	在这种应用情景中，用到了两个refrence panel，分别记作 panel 0 和 panel 1

	例如，panel 0 可以是1000 Genomes Project的haplotype，包含了基因组中几乎全部的常见SNPs；panel 1 可以是HapMap Phase 3的haplotype，仅包含了基因组中的部分的常见SNPs；panel 3 是用商用SNPs芯片得到的一系列的case和control的样本的genotype

<a name="tools-impute2-best-practices"><h3>6.1.2. Best Practices [<sup>目录</sup>](#content)</h3></a>

1. **基因型填充前 (pre-imputation) 进行genotypes质控**

	过滤低质量的变异位点和样本

2. **保证分析中使用的基因组坐标系统一致**

	NCBI build number (e.g., "b36" or "b37") 对应于 UCSC version (e.g., "hg18" or "hg19")

3. **选择reference panel**

	之前的GWAS研究中，研究人员一般都是选择与对应人群遗传距离最相近的reference panel，而Impute2推荐**使用worldwide reference panel**，程序能够从中选出最合适的haplotype用于基因型填充

	这样做的好处是：

	> - 不需要费时费力去挑选haplotypes来构造reference panel；
	> 
	>
	> ```
	> Good results can be obtained in any study population by tuning a single software parameter (-k_hap) with a simple rule of thumb
	>```
	>
	> - 该策略适用于各种人群的研究；
	>
	> ```
	> Our group and others have used this approach to successfully impute populations ranging from homogeneous isolates to recent and
	> complex admixtures
	> ```
	> 
	> - 填充效果往往优于研究人员自己挑选构造的小reference panel
	> 
	> ```
	> This is because individuals from "diverged" populations may still share genomic segments of recent common ancestry, and IMPUTE2
	>  can use this haplotype sharing to improve accuracy. At the same time, the software can ignore haplotypes that are not helpful.
	> ```
	> 
	> - 对于大的reference panel，Impute2也能进行高效地处理，不需要担心会带来的计算负担

<a name="tools-impute2-1kgp-cookbook"><h3>6.1.3. 1000 Genomes Imputation Cookbook [<sup>目录</sup>](#content)</h3></a>

<a name="tools-impute2-1kgp-cookbook-before-imputation"><h4>6.1.3.1. Before Imputation [<sup>目录</sup>](#content)</h4></a>

1. **对Genotype数据进行质控**

	包括样本水平的质控和marker水平的质控
	
	> - 样本水平：
	> 	- call rate
	> 	- heterozygosity
	> 	- relatedness between genotyped individuals
	> 	- correspondence between sex chromosome genotypes and reported gender
	> - marker水平：
	> 	- call rates
	> 	- deviations from Hardy-Weinberg Equilibrium
	> 	- excluding low frequency SNPs，for older genotyping platforms,
	
	质控代码参考自：`http://sites.google.com/site/mikeweale/software/gwascode`

2. **将Genotype数据转换为Build 37**

	目前的1000 Genome Project的数据使用的是NCBI genome build 37 (hg19)的坐标系统，因此在基因型填充之前需要保证你的Genotype文件也是hg19的坐标系统，且位点是落在正链上

	若坐标系统不一致，可以使用LiftOver进行坐标转换，但是转换过程中可能有少量的SNP转换失败

3. **将Genotype文件转换为IMPUTE格式**

	在格式转换之前，需要先按照坐标进行排序

	`GTOOL`可以将`PLINK PED`转换为`IMPUTE`格式



<a name="tools-impute2-1kgp-cookbook-pre-phasing"><h4>6.1.3.2. Pre-Phasing [<sup>目录</sup>](#content)</h4></a>

对于大规模的reference panels，基因型填充建议分两步进行：

> - pre-phasing：推断每个样本的单体型
> - imputation：对分型得到的单体型 (phased haplotypes) 中缺失的allele进行基因型填充

`IMPUTE2` 或 `SHAPEIT` 都可以执行pre-phasing操作，Drs. Bryan Howie 和 Jonathan Marchini推荐使用`SHAPEIT`进行pre-phasing，因为该工具采用的phasing方法更准确

pre-phasing采用**滑动窗口法** (Sliding Window Analyses) 进行：

> 将一个染色体分割成若干Mb的块 (blocks)，对这个块中的genotypes进行phasing
> 
> 用这种phasing方法可能遇到的**两种棘手的情况**：
> 
> - 若每个块的大小一致，可能不同的块之间SNP的分布密度存在很大差异；
> - 若某个块正好跨过着丝点 （centromere）,而着丝点附近区域的SNP密度极低，往往意味着是 large gap in 1000 Genome SNPs；
> 
> 若SNP密度太低是很难进行phasing的，此时可以采取的**解决策略**是：
> 
> - 若某一个块的SNP密度太低（例如少于200个），可以将它与邻居的块合并成一个更大的块一起phasing；
> - 避免构造跨着丝点的块；

<a name="tools-impute2-1kgp-cookbook-pre-phasing-using-impute2"><h4>6.1.3.2.1. using IMPUTE2 [<sup>目录</sup>](#content)</h4></a>


<a name="tools-impute2-1kgp-cookbook-pre-phasing-using-shapeit"><h4>6.1.3.2.2. using SHAPEIT (recommended) [<sup>目录</sup>](#content)</h4></a>



---

参考资料：

(1) [【简书】群体遗传学习笔记-基因型缺失数据的填充](https://www.jianshu.com/p/dafd1e6e4a98)

(2) [Impute2官方文档](http://mathgen.stats.ox.ac.uk/impute/impute2_overview.html)

(3) [Genotype File Format](http://www.stats.ox.ac.uk/%7Emarchini/software/gwas/file_format.html#Genotype_File_Format)

(4) [IMPUTE2: 1000 Genomes Imputation Cookbook](https://genome.sph.umich.edu/wiki/IMPUTE2:_1000_Genomes_Imputation_Cookbook)

(5) Weale M (2010) Quality Control for Genome-Wide Association Studies. Methods Mol. Biol. 628:341–372
