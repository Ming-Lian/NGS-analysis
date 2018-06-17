<a name="content">目录</a>

[Algorithms in Bioinformatics](#title)
- [Sequence Alignment](#alignment)
	- [DIAMOND](#diamond)









<h1 name="title">Algorithms in Bioinformatics</h1>

<a name="alignment"><h2>Sequence Alignment [<sup>目录</sup>](#content)</h2></a>

<a name="diamond"><h3>DIAMOND [<sup>目录</sup>](#content)</h3></a>

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








参考资料：

(1) Benjamin Buchfink, Chao Xie & Daniel H. Huson, Fast and Sensitive Protein Alignment using DIAMOND, Nature Methods, 12, 59–60 (2015) doi:10.1038/nmeth.3176.
