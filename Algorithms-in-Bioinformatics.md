<a name="content">目录</a>

[Algorithms in Bioinformatics](#title)
- [Sequence Alignment](#alignment)
	- [DIAMOND](#alignment-diamond)
	- [BLAT](#blat)
	- [BLASR](#blasr)
- [Motif Finding](#motif)
	- [MEME：EM算法](#motif-em)
- [Bining for Metageonome](#bining)
	- [CONCOCT](#bining-concoct)
- [Genome Assembly](#genome-assembly)
	- [构建de Bruijn graph](#assembly-construct-dbg)
	- [简化DBG得到assembly graph：filigree edges](#assembly-simplify-graph)







<h1 name="title">Algorithms in Bioinformatics</h1>

<a name="alignment"><h2>Sequence Alignment [<sup>目录</sup>](#content)</h2></a>

<a name="alignment-diamond"><h3>DIAMOND [<sup>目录</sup>](#content)</h3></a>

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

<a name="blat"><h3>BLAT [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/Algorithms-Bioinf-alignment-BLAT.png width=600/></p>

> Example showing the creation of non-overlapping k-mers from the target database and overlapping k-mers from the query sequence, for k=3. Coordinates of the database sequences are used to clump the matches into larger alignments (full process not shown).

<a name="blasr"><h3>BLASR [<sup>目录</sup>](#content)</h3></a>

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

<a name="motif"><h2>Motif Finding [<sup>目录</sup>](#content)</h2></a>

<a name="motif-em"><h3>MEME：EM算法 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/MeRIP-seq-meme-principle-1.png width=800/></p>

<p align="center"><img src=./picture/MeRIP-seq-meme-principle-2.png width=800/></p>

<p align="center"><img src=./picture/MeRIP-seq-meme-principle-3.png width=800/></p>

<a name="bining"><h2>Bining for Metageonome [<sup>目录</sup>](#content)</h2></a>

<a name="bining-concoct"><h3>CONCOCT [<sup>目录</sup>](#content)</h3></a>

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

<a name="genome-assembly"><h2>Genome Assembly [<sup>目录</sup>](#content)</h2></a>

<a name="assembly-construct-dbg"><h3>构建de Bruijn graph [<sup>目录</sup>](#content)</h3></a>


<a name="assembly-simplify-graph"><h3>简化DBG得到assembly graph：filigree edges [<sup>目录</sup>](#content)</h3></a>

由原始reads推断出的de Bruijn graph 不能就直接用于拼接，需要进行剪枝

传统的拼接工具通过设定一个全局的reads覆盖度阈值来将那些只有较少reads支持的分支当做是测序错误，将这些分支剪去。这种图简化方法在一般的单一物种基因组拼接任务中是行得通的，但是应用在metagenome的拼接中则明显不合适：

metagenome中一种菌往往以strain mixtures形式出现，且不同菌株之间丰度差异很大，当我们想忽视strain mixtures中的菌株间的差异，构造出一个代表该菌的consensus genome时，则我们将图中代表稀有菌株的edge剪去，而保留菌种共有的edge

<p align="center"><img src=./picture/Algorithms-Bioinf-assembly-simplify-graph-1.png height=600 /></p>

采用相邻edge的coverage ratio来剪枝

<p align="center"><img src=./picture/Algorithms-Bioinf-assembly-simplify-graph-2.png ></p>

if ratio * cov (e<sub>i</sub>) < cov (v)，去除该edge



参考资料：

(1) Benjamin Buchfink, Chao Xie & Daniel H. Huson, Fast and Sensitive Protein Alignment using DIAMOND, Nature Methods, 12, 59–60 (2015) doi:10.1038/nmeth.3176.

(2) Alneberg, J. et al. Binning metagenomic contigs by coverage and composition. Nat. Methods 11, 1144–1146 (2014).

(3) Beaulaurier J, Zhu S, Deikus G, et al. Metagenomic binning and association of plasmids with bacterial host genomes using DNA methylation.[J]. Nature Biotechnology, 2017, 36(1).

(4)  Nurk S., Meleshko D., Korobeynikov A., Pevzner P. A. metaSPAdes: a new versatile de novo metagenomics assembler.	Genome Research, 2017 
