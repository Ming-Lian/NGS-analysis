<a name="content">目录</a>

[Algorithms in Bioinformatics](#title)
- [1. Sequence Alignment](#alignment)
	- [1.1. DIAMOND](#alignment-diamond)
	- [1.2. BLAT](#blat)
	- [1.3. BLASR](#blasr)
	- [1.4. BLAST](#blast)
		- [1.4.1. TCR/BCR克隆鉴定](#tcr-bcr-identification)
	- [1.5. 基于后缀树的快速序列比对](#fast-alignment-based-on-suffix-tree)
- [2. Motif Finding](#motif)
	- [2.1. MEME：EM算法](#motif-em)
- [3. Bining for Metageonome](#bining)
	- [3.1. CONCOCT](#bining-concoct)
- [4. Genome Assembly](#genome-assembly)
	- [4.1. 构建de Bruijn graph](#assembly-construct-dbg)
	- [4.2. 简化DBG得到assembly graph：filigree edges](#assembly-simplify-graph)
- [5. variants calling](#variants-calling)
	- [5.1. 背景知识](#introduction-to-variants-calling)
		- [5.1.1. base calling过程的错误及校正](#base-calling-error-and-correction)
		- [5.1.2. 比对过程的错误及校正](#mapping-error-and-correction)
	- [5.2. snp calling](#snp-calling)
		- [5.2.1. samtools/bcftools](#snp-calling-using-samtools-bcftools)

<h1 name="title">Algorithms in Bioinformatics</h1>

<a name="alignment"><h2>1. Sequence Alignment [<sup>目录</sup>](#content)</h2></a>

<a name="alignment-diamond"><h3>1.1. DIAMOND [<sup>目录</sup>](#content)</h3></a>

DIAMOND is a new high-throughput program for aligning DNA reads or protein sequences against a protein reference database such as NR, at up to 20,000 times the speed of BLAST, with high sensitivity. 

> In fast mode, about 20,000 times faster than BLASTX, while reporting about 80-90% of all matches that BLASTX finds, with an e-value of at most 1e-5
> 
> In sensitive mode, DIAMOND ist about 2,500 times faster than BLASTX, finding more than 94% of all matches

**算法原理探究：**

常规的序列比对工具，以 BLASTX 为代表，其采用的是 seed-and-extend 的策略

> In this two-phase approach, users search first for matches of seeds (short stretches of the query sequence) in the reference database, and this is followed by an 'extend' phase that aims to compute a full alignment

在进行比对之前，需要对整个参考序列建好种子索引，即 k-mer seeds index，然后在对用户提交的查询序列集进行线性扫描。

而 DIAMOND 对此进行了改进：

> **1\. Double Index**
> 
> 对参考序列与查询序列同时建立 index，得到它们的 seeds 和 locations 的列表。对这两个列表进行字典排序和关联。
> 
> ```
> Double indexing takes advantage of the cache hierarchy by increasing data locality, 
> thus reducing the demands on main memory bandwidth.
> ```
> 
> **2\. Spaced Seeds**
> 
> BLASTX 使用步移法获得连续的 k-mer seeds，这使得为了获得较高的灵敏度 (sensitivity) 要求 seeds 足够短 (length 3–6 amino acids)，这增加了搜索的时间复杂度。为了既能保证高灵敏度又能使用长 seeds，DIAMOND 使用了 spaced seeds，即在部分位置使用更长的 seeds。
> 
> 这些位点的数量与分布被分别叫作 spaced seeds 的权重 (weight) 和 形状 (shape)。
> 
> ```
> To achieve high sensitivity, DIAMOND uses a set of four carefully chosen shapes11 of length 15–24 and 
> weight 12 by default. The most sensitive version of DIAMOND uses 16 shapes of weight 9.
> ```

<a name="blat"><h3>1.2. BLAT [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/Algorithms-Bioinf-alignment-BLAT.png width=600/></p>

> Example showing the creation of non-overlapping k-mers from the target database and overlapping k-mers from the query sequence, for k=3. Coordinates of the database sequences are used to clump the matches into larger alignments (full process not shown).

<a name="blasr"><h3>1.3. BLASR [<sup>目录</sup>](#content)</h3></a>

BLASR是第一个针对PacBio序列的比对工具，2012年发表在《BMC Bioinformatics》期刊上，由PacBio研究团队开发，并且一直在更新，目前Google引用次数为433（截止2018.08.14）。

其主要思想如下图所示：

<p align="center"><img src=./picture/Algorithms-Bioinf-alignment-BLASR-1.png width=800 /></p>

包含三步：

> - （A）首先通过BWT-FM压缩或后缀数组（suffix array）索引技术对基因组进行转换，寻找与待比对序列相似性比较高的候选区域或区间；
> - （B）然后对这些候选区间进行稀疏动态比对（sparse dynamic programming）得到初步比对结果；
> - （C）最后运用动态规划算法进行详细的序列比对，得到最终的比对结果。

由于BLASR运用BWT-FM索引技术，因而可以大大提高搜索速度，并且降低内存消耗

下图是BLASR与传统二代比对工具的运行比较结果。可以看出不论运行速度方面还是消耗内存方面，BLASR均优于BWA-SW算法

<p align="center"><img src=./picture/Algorithms-Bioinf-alignment-BLASR-2.png width=800 /></p>

<a name="blast"><h3>1.4. BLAST [<sup>目录</sup>](#content)</h3></a>

<a name="tcr-bcr-identification"><h4>1.4.1. TCR/BCR克隆鉴定 [<sup>目录</sup>](#content)</h4></a>

采用两轮比对的策略，下图以鉴定V基因亚型为例：

<p align="center"><img src=./picture/Algorithms-Bioinf-alignment-BLAST-TCR-1.png width=800 /></p>

两轮比对的合理性分析：

> 第一轮比对是采用**全局比对**策略来鉴定 Non-CDR3 region；
>
> 第二轮比对是采用比对起始点固定（为第一轮鉴定出的CDR3起始点）的**局部比对**来鉴定CDR3 region；
>
> 为什么鉴定 Non-CDR3 region使用的是全局比对 (global alignment)，而鉴定CDR3 region使用的是局部比对 (local alignment)？
>
> 这是由TCR/BCR在进行V(D)J基因重组方式决定的，重组是从V、(D)、(J)基因座中各随机选择一个进行重排，得到 $V_i-(D_j)-J_k$形式的克隆重组形式，同时，会在两个基因的连接处发生随机的Indel，即在$V_i$的3'端、$D_j$的5'和3'端、$J_k$的5'端发生Indel，而在非CDR3区域，即$V_i$的5'端和$J_k$的3'端不发生改变
>
> 所以Non-CDR3 region理论上与reference高度一致，使用全局比对即可，而CDR3 region变异比较强，要使用更为敏感的局部比对

那么，在对CDR3区域进行局部比对时，设置多少的mismatch合理呢？

根据现有数据的mismatch的统计来设置阈值

<p align="center"><img src=./picture/Algorithms-Bioinf-alignment-BLAST-TCR-2.png width=800 /></p>

因此阈值设定为：

> - TRBV/J, TRAV/J：0；
>
> - IGHV/J：2；
>
>- IGKV/J, IGLV/J： 7；

根据CDR3区域的完整性（考虑碱基长度是否为3的倍数，即考虑是否存在移码突变，以及终止密码子导致的无义突变），鉴定出来的CDR3区域会被分为两大类：

> - in frame
>
> - out of frame

<a name="fast-alignment-based-on-suffix-tree"><h3>1.5. 基于后缀树的快速序列比对 [<sup>目录</sup>](#content)</h3></a>

一个长度为n的字符串S，它的后缀树定义为一棵满足如下条件的树：

> 1. 从根到树叶的路径与S的后缀一一对应。即每条路径惟一代表了S的一个后缀；
>
> 2. 每条边都代表一个非空的字符串；
>
> 3. 所有内部节点（根节点除外）都有至少两个子节点。

由于并非所有的字符串都存在这样的树，因此S通常使用一个终止符号进行填充（通常使用$）。 

以字符串S=banana为例，来说明如何构建后缀树：

（1）首先列出S的所有后缀子字符串：

- banana$
- anana$
- nana$
- na$
- a$
- $

（2）根据上面列出的后缀子字符串，构建suffix trie（后缀字典树）




















<a name="motif"><h2>2. Motif Finding [<sup>目录</sup>](#content)</h2></a>

<a name="motif-em"><h3>2.1. MEME：EM算法 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/MeRIP-seq-meme-principle-1.png width=800/></p>

<p align="center"><img src=./picture/MeRIP-seq-meme-principle-2.png width=800/></p>

<p align="center"><img src=./picture/MeRIP-seq-meme-principle-3.png width=800/></p>

<a name="bining"><h2>3. Bining for Metageonome [<sup>目录</sup>](#content)</h2></a>

<a name="bining-concoct"><h3>3.1. CONCOCT [<sup>目录</sup>](#content)</h3></a>

结合序列组成特征 (sequence composition) 和跨样本覆盖度特征 (coverage across multiple samples) 进行bining

在进行bining之前需要将所有样本的reads进行混拼 (coassembly) 得到contigs

- Sequence composition features

以 k-mer 长度等于5为例

将相互之间成反向互补关系的 5-mers pairs 记做一种，则总共有512种 5-mers，对每一条contig计算其各自 5-mers 的组成频率从而构造出一个长度为v=512的向量 Z<sub>i</sub>：

<p align="center">Z<sub>i</sub> = (Z<sub>i,1</sub>, ..., Z<sub>i,v</sub>)</p>

为了保证每个 5-mers 频率的计数非零（为后面的对数转换做准备），进行伪计数处理，然后用该序列的 5-mers的总数进行标准化，得到新的向量 Z<sub>i</sub><sup>'</sup>：

<p align="center"><img src=./picture/Algorithms-Bioinf-bining-composition-formula.png height=80 /></p>

- Contig coverage features

用段序列比对软件，将各个样本（总共有M个样本）的reads比对到contigs上，计算每条contigs在每个样本中的 coverage (Mapped reads * read length / contig length)，得到表示 congtig i 的 coverage 的向量 Y<sub>i</sub>：

<p align="center">Y<sub>i</sub> = (Y<sub>i,1</sub>, ..., Y<sub>i,M</sub>)</p>

两轮标准化处理：

> 伪计数处理：额外添加一条比对到该 contig 上的 read；再用该**样本内**所有 contigs（contigs总数为N）的 coverage 进行标准化，得到新的向量 Y<sub>i</sub><sup>'</sup>：
> 
> <p align="center"><img src=./picture/Algorithms-Bioinf-bining-coverage-formula-1.png height=80 /></p>
> 
> 然后再在**contig内部**进行标准化，得到新的向量 Y<sub>i</sub><sup>''</sup>：
> 
> <p align="center"><img src=./picture/Algorithms-Bioinf-bining-coverage-formula-2.png height=80 /></p>

- Combine two features

将表示某一个 contig i 的序列组成特征的向量 Z<sub>i</sub><sup>'</sup> 和 coverage 特征向量 Y<sub>i</sub><sup>''</sup>合并成组成一个新的特征向量 X<sub>i</sub>（向量长度为E=V+M），同时进行对数转换：

<p align="center">X<sub>i</sub> = { log(Z<sub>i</sub><sup>'</sup>) , log(Y<sub>i</sub><sup>''</sup>) }

- PCA，保留能解释至少90%的方差的主成分（共保留前D个主成分，D < E）

矩阵维数变化：N x E => N x D

- cluster contigs into bins

使用高斯混合模型 (Gaussian mixture model)

<a name="genome-assembly"><h2>4. Genome Assembly [<sup>目录</sup>](#content)</h2></a>

<a name="assembly-construct-dbg"><h3>4.1. 构建de Bruijn graph [<sup>目录</sup>](#content)</h3></a>


<a name="assembly-simplify-graph"><h3>4.2. 简化DBG得到assembly graph：filigree edges [<sup>目录</sup>](#content)</h3></a>

由原始reads推断出的de Bruijn graph 不能就直接用于拼接，需要进行剪枝

传统的拼接工具通过设定一个全局的reads覆盖度阈值来将那些只有较少reads支持的分支当做是测序错误，将这些分支剪去。这种图简化方法在一般的单一物种基因组拼接任务中是行得通的，但是应用在metagenome的拼接中则明显不合适：

metagenome中一种菌往往以strain mixtures形式出现，且不同菌株之间丰度差异很大，当我们想忽视strain mixtures中的菌株间的差异，构造出一个代表该菌的consensus genome时，则我们将图中代表稀有菌株的edge剪去，而保留菌种共有的edge

<p align="center"><img src=./picture/Algorithms-Bioinf-assembly-simplify-graph-1.png height=600 /></p>

采用相邻edge的coverage ratio来剪枝

<p align="center"><img src=./picture/Algorithms-Bioinf-assembly-simplify-graph-2.png ></p>

if $ratio\times cov(e_i) < cov(v)$，去除该edge

<a name="variants-calling"><h2>5. variants calling [<sup>目录</sup>](#content)</h2></a>

<a name="introduction-to-variants-calling"><h3>5.1. 背景知识 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/Algorithms-Bioinf-variants-calling-1.png ></p>

变异检测的常规步骤：

> - 将一个或多个样本的reads比对到参考基因组；
> - SNP calling：检出变异位点；
> - genotype calling：鉴定出个体的每个变异位点的基因型 (genotype)；

影响变异检测准确性的因素：

> - base-calling 的错误率；
> - 比对 (alignment) 的错误率；
> - 低覆盖度的测序 (low-coverage sequencing, <5× per site per individual, on average)，这使得：对于二倍体个体的两条同源染色体的某个位点，有很大概率只采样到其中一条染色体；

变异检测的准确率会影响到下游的分析，包括：

> - 鉴定罕见变异 (rare mutations)
> - 评估 allele 频率
> - 相关性分析 (association mapping)

一个提高准确率的策略是进行**靶向深度测序**（sequence target regions deeply (at >20× coverage)）

但是，随着大样本检测需求的增加，中覆盖度 (5-20X) 或 低覆盖度 (<5X) 是一个更加经济的选择。而且对于人群中低频变异的检测，大样本是一个重要前提，低覆盖度对低频变异检测的影响不大；对于相关性分析 (GWAS)，依赖于大样本和每个样本鉴定的准确性，但是往往这两者是不可兼得的，那么相对于单样本鉴定的准确性来说，大样本更加重要

许多用于提高变异检测准确性和评估准确率的方法被陆续提出，这些方法大都采用了概率统计模型

> 对一个基本概率统计量 "genotype likelihoods" 进行建模
>
> genotype likelihoods：包含了对 base calling，alignment 和 assembly 步骤中的错误率的综合评估，它利用了 allele frequencies、LD（连锁不平衡）模式等先验信息 (prior information)
>
> 最后给出的分析结果包括：
>
> - SNP call 和对应的不确定性评估
> - genotype call 和对应的不确定性评估
>
> 它们都有具体的统计学意义

<a name="base-calling-error-and-correction"><h4>5.1.1. base calling过程的错误及校正 [<sup>目录</sup>](#content)</h4></a>

- **不同测序仪的错误类型及产生的原因**

	（1）454测序仪——indel错误

	<p align="center"><img src=./picture/Algorithms-Bioinf-variants-calling-2.png ></p>

	454采用的是焦磷酸测序法，每一轮测序加入一种 ddNTP 作为反应的底物，即每一轮反应连接上的碱基组成是已知的，如果这种 ddNTP 可以在当前位点与模板链互补则成功连接延伸，释放出焦磷酸而发出荧光，从而测出当前位点的碱基组成，否则不发荧光

	但是如果在当前位点和它之后的若干个位点的碱基组成相同，即是同聚物形式，则当前这轮的测序反应是这几个连续位点的连续反应，由于这个连续反应之间的时间间隔几乎可以忽略不仅，则相当于一次性释放出多个焦磷酸，产生比单个反应强得多的荧光

	由于前面已经提到，这轮反应检测的碱基组成是已知的，所以454测序的 base calling 过程要做的是区别这种同聚反应的同聚物的长度，理论上长度越长，荧光强度越高，但是实际上，同聚物长度与荧光强度之间并不存在稳定的正比关系，相同长度的同聚物产生的荧光强度具有较大的方差，这使得 base calling 过程容易产生较高比例的 indel 错误

	（2）Illumina测序仪

	<p align="center"><img src=./picture/Algorithms-Bioinf-variants-calling-3.png ></p>

	Illumina测序仪由于测序原理与454不同，它采用的是边合成边测序的原理，每成功加上一个碱基之后，由于刚加上的 ddNTP 的3号位被连接上了荧光基团，阻止了下一个 ddNTP 的连接，只有在被检测完荧光并且将占位的荧光基团切除后才能进行后面的合成反应，所以每次只能最多只有一个碱基的合成，因此不存在454中所谓的 indel 错误

	Illumina 测序仪的 base calling 的错误主要来自于根据检测到的荧光推断出当前连接的碱基的组成这一步，称为 miscall，在当时（2010年附近），Illumina 测序仪的 miscall error rate 在 1% 附近

	导致 miscall 的主要原因是同一个cluster（在测序芯片即flowcell上，一个cluster一般来说是来源于同一个ssDNA片段的扩增产物）中的不同的DNA片段合成过程的**不同步**

	> 在每一轮反应中，大多数DNA片段是同步的，只有少数片段本来应该进行连接延伸，但是由于一些原因，比如部分反应空间ddNTP的浓度偏低，导致反应成功率下降，使得它的反应没有成功进行而相对于这个cluster的其他片段滞后了一些，或者，在曝光检测之后的荧光切除过程中有少部分的DNA片段上的荧光基团没有被成功切除，导致了后续反应的滞后，一般来说，第二种情况发生的可能性会大一些
	>
	> 因为每一轮不同步滞后的那些DNA片段毕竟占少数，所以，对这个cluster整体的荧光信息的影响比较小，所以在测序反应刚开始的时候，base calling 的准确率是比较高的；但是随着反应的持续进行，不同步部分不断累积越来越多，那么干扰信息也就越来越强，也就越来越难判断当前碱基组成到底是哪个，这也就解释了为什么反应越到后面测序质量越差

- **Phred 质量值校正**

	碱基质量的表示方式：Phred quality score (Q score)

	$$Q_{Phred}=-10\log_{10} P(error)$$

	则可以根据Q值算出它的实际错误率：

	$$P(error)=10^{-Q/10}$$

	若 Phred score 等于20，意味着测序错误率为 1%

	Phred碱基质量值是由base-calling算法评估出来的，但是它们可能并不能准确地反映真实的base-calling错误率

	> 针对不同测序平台提出 的 base calling 算法，包括：
	>
	> - 454 —— Pyrobayes
	>
	> - SOLiD —— Rsolid
	>
	> - Illumina —— Ibis、BayesCall
	>
	> 这些 base calling 算法基本都是测序仪厂家开发出来的，将原始测序错误率降低了 ~5-30%

	（1）SOAPsnp的校正策略

<a name="mapping-error-and-correction"><h4>5.1.2. 比对过程的错误及校正 [<sup>目录</sup>](#content)</h4></a>

（1）比对过程中可能存在的问题

比对过程需要区分测序错误和实际的碱基差异，同时给出合适的比对质量值 (well-calibrated mapping quality values)，因为后续的 variant calling 需要利用到它来计算后验概率 (posterior probabilities)

在进行比对的过程中，选择合适的mismatch位点数是一个比较重要的问题，它需要在准确性和read depth之间做权衡 (trade-off)

> 不同物种它适用的可容忍的mismatch数是不一样的
>
> 例如，果蝇不同个体之间的基因组变异程度相对于我们人来说是比较大的，如果用适用于人基因组mapping的可容忍mismatch数来进行果蝇基因组的mapping，明显是过于严格的，这样会导致大量的reads无法成功比对到参考基因组上，特别是那些多态性位点分布密度比较高的区域，reads基本上很难比对上，这些区域的变异的检出率也就比较低
>
> 再比如，将果蝇的参数用在人上面，那么比对的标准相对设低了，这就容易带来大量的错误比对结果

对于基因组中的高变区域，reads的mapping很困难，一个解决方案是使用更长的reads或者双端测序，但是对于diversity极高的区域，比如 MHC (major histocompatibility complex ) 区域，这些方法还是显得力不从心，这个时候就采用 de novo 拼接的策略了

结合长reads和de novo 拼接的方法，基因组中高变区域鉴定的大多数情况都可以得到有效地解决



























<a name="snp-calling"><h3>5.2. snp calling [<sup>目录</sup>](#content)</h3></a>

<a name="snp-calling-using-samtools-bcftools"><h4>5.2.1. samtools/bcftools [<sup>目录</sup>](#content)</h4></a>

（1）李恒对千人基因组计划的思考：

- 为什么选择低覆盖多样本的策略而非高覆盖度少样本的策略？

	如果采用的是高覆盖的的测序策略，那么基本上可以利用reads之间的相互校验来得到更为准确的信息

	但是如果采用的是低覆盖多样本的策略，我们可以减少采样的波动/不稳定性，可以鉴定出存在多个样本中的variants，且能得到许多在群体中罕见的variants

	但是低覆盖度也不能太低，如果低到无法区分出到底是测序错误还是实际变异，也不行。经过摸索，得出结论，每个样本的覆盖度在2-6x之间比较合理

- 如何对低覆盖度的样本得到相对准确的variants discovery？

	首先可以确定的是，如果直接对这个样本进行variants calling，得到的结果不仅不全（因为大部分基因组区域都没有被测序到），而且不准确（即使那些被测序到的区域，depth也很低，往往无法给出准确的变异鉴定）

	目前对于这样的数据，常用的分析策略是将所有的样本混合起来一起鉴定variants——这样做的确能得到更为准确的群体的变异信息，但是无法得到群体中每个样本具体的variants，而准确的个体变异信息是许多群体遗传学分析的基础（比如Hardy–Weinberg equilibrium (HWE) 测试，GWAS等等）

	那么，如何对低覆盖度的样本得到相对准确的variants discovery呢？

	可以利用位点之间的连锁不平衡（LD）来进行基因型填充（genotypes imputation）

	> 例如，对于某一个位点A，某一个样本（称为当前样本）在这个位点上覆盖度太低，无法进行variant calling。但是能在这个群体中找到其他的样本在该位点有足够的覆盖度，且通过分析发现A位点与B位点存在高度的连锁，若当前样本在B位点有高覆盖度，即能给出准确的B位点的genotype，则通过A-B位点之间的连锁关系，可以给出当前样本A位点genotype的相对更加准确的推断

	但是采用imputation策略来辅助variants discovery也不是没有问题

	> - imputation的效率取决于LD的模式（即不同人群之间的LD模式可能不同，要么有连锁，要么无连锁，要么强连锁，要么弱连锁），使用了不太恰当的LD模式作为imputation的reference可能引入潜在的bias；
	>
	> - 目前的imputation的算法效率低速度慢，在大人群的variants discovery中往往成了分析中的瓶颈，若样本更多、采用更为准确的imputation算法，会更慢；

	这使得李恒他们开始思考：对于低覆盖度的样本，imputation真的是目前可选的最优方案吗？

- 相关样本间的 somatic mutation 或 germline mutation 的检测中存在的挑战

	相关样本间的变异率在 $10^{−6} \sim 10^{-7}$ 之间，与测序错误率（不是实际的测序错误率，实际测序错误率在$\sim 10^{-3}$，在高覆盖度的测序数据中通过reads之间的相互校正可以得到更低的错误率： $\sim 10^{-5}$） 在一个相近的数量级，这使得如何区别实际 variants 和测序错误变得十分困难

	对于rare variants的检测本来就是variants discovery中的难题，常规的variants calling工具在面对这样的问题时往往会失效

	从另外一个角度来说，比较相关样本间的rare mutation不是我们的目的，我们的目的是想通过它来分析相关样本间的遗传差异，所以genotype只是遗传差异的一种衡量而已，如果能找到更好的衡量指标，完全可以不进行genotyping

（2）samtools/bcftools进行variants calling的逻辑

符号说明：

| Symbol | Description |
|:---|:---|
| $n$	| Number of samples |
| $m_i$	| Ploidy of the $i$-th sample ($1≤i≤n$) |
| $M$	| Total number of chromosomes in samples: $M=\sum_i m_i$ |
| $d_i$ | Sequencing data (bases and qualities) for the $i$-th sample |
| $g_i$ | Genotype (the number of reference alleles) of the $i$-th sample ( $0≤g_i≤m_i$ ) |
| $\phi_k$ | Probability of observing k reference alleles ( $\sum_{k=o}^M \phi_k=1$ ) |
| $Pr\{A\}$ | Probability of an event A |
| $L_i(\theta)$ | Likelihood function for the $i$-th sample: $L_i(\theta)=Pr\{d_i \mid \theta\}$ |

前提假设：

> - 不同位点间相互独立；
> - 对于同一个位点，不同reads的测序错误或mapping误差相互独立；
> - 只考虑二等位情况；

对于某一个样本的某一个位点，有$k$条 reads 比对上，其中有$l$条 ($0 \le l\le k$) 序列在该位点的碱基组成与 reference 一致，剩余 $k-l$ 条与 reference 不同，其中第 $j$ 条上该碱基的测序错误率为 $\epsilon_j$，则该样本的基因型为 $g \in \{1,2\}$ （$g=1$表示该样本在位点为纯合型，但可能与reference一致，也可能与reference不同；$g=2$表示该位点为杂合型） 的概率为

$$
\begin{aligned}
&\quad L(g) \\
&= Pr(d \mid g) \\
&= \prod_{j=1}^l Pr_j(A)\prod_{j=l+1}^k Pr_j(\overline A) \\
&= \prod_{j=1}^l [Pr_j(B \mid A)+Pr_j(\overline B \mid A)]\prod_{j=l+1}^k [Pr_j(B \mid \overline A)+Pr_j(\overline B \mid \overline A)] \\
&= \prod_{j=1}^l \left[ \frac{g}{m}(1-\epsilon_j) + \frac{m-g}{m}\epsilon_j \right] \prod_{j=l+1}^k \left[  \frac{m-g}{m}(1-\epsilon_j) +  \frac{g}{m}\epsilon_j\right] \\
&= \frac{1}{m^k}\prod_{j=1}^l [g(1-\epsilon_j) + (m-g)\epsilon_j] \prod_{j=l+1}^k [(m-g)(1-\epsilon_j) + g\epsilon_j]
\end{aligned}
$$

其中，$m$ 是该物种的倍性，普通人是二倍体，因此一般 $m=2$

事件$A=\{碱基与ref一致\}$，则$\overline A=\{碱基与ref不一致\}$

事件$B=\{该碱基的测序是正确的\}$，则$\overline B=\{该碱基的测序是错误的\}$



















---

参考资料：

(1) Benjamin Buchfink, Chao Xie & Daniel H. Huson, Fast and Sensitive Protein Alignment using DIAMOND, Nature Methods, 12, 59–60 (2015) doi:10.1038/nmeth.3176.

(2) Chaisson M J, Tesler G. Mapping single molecule sequencing reads using basic local alignment with successive refinement (BLASR): application and theory[J]. BMC bioinformatics, 2012, 13(1): 238.

(3) [生信算法《三代测序序列比对利器-BLASR，更小更快更方便 》](https://mp.weixin.qq.com/s/7xZVvyShcwjPGbgKiLniag)

(4) Alneberg, J. et al. Binning metagenomic contigs by coverage and composition. Nat. Methods 11, 1144–1146 (2014).

(5) Beaulaurier J, Zhu S, Deikus G, et al. Metagenomic binning and association of plasmids with bacterial host genomes using DNA methylation.[J]. Nature Biotechnology, 2017, 36(1).

(6)  Nurk S., Meleshko D., Korobeynikov A., Pevzner P. A. metaSPAdes: a new versatile de novo metagenomics assembler.	Genome Research, 2017 

(7) Nielsen R, Paul JS, Albrechtsen A, Song YS. Genotype and SNP calling from next-generation sequencing data. Nat Rev Genet. 2011;12(6):443–451.

(8) Li H. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 2011;27(21):2987–2993.

(9) [Samtools math notes](http://www.broadinstitute.org/gatk/media/docs/Samtools.pdf)
