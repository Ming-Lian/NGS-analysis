<a name="content">目录</a>

[Algorithms in Bioinformatics](#title)
- [1. Sequence Alignment](#alignment)
	- [1.1. DIAMOND](#alignment-diamond)
	- [1.2. BLAT](#blat)
	- [1.3. BLASR](#blasr)
	- [1.4. BLAST](#blast)
		- [1.4.1. TCR/BCR克隆鉴定](#tcr-bcr-identification)
	- [1.5. 基于后缀树的快速序列比对](#fast-alignment-based-on-suffix-tree)
		- [1.5.1. Trie树](#fast-alignment-based-on-suffix-tree-trie-tree)
        		- [1.5.1.1. Trie树及基本操作](#fast-alignment-based-on-suffix-tree-trie-tree-1)
			- [1.5.1.2. Trie树的应用](#fast-alignment-based-on-suffix-tree-trie-tree-2)
		- [1.5.2. Suffix Trie树](#fast-alignment-based-on-suffix-tree-suffix-trie-tree)
	- [1.6. Subread: seed-and-vote](#subread-alignment-algorithmn)
- [2. Motif Finding](#motif)
	- [2.1. MEME：EM算法](#motif-em)
	- [2.2. position weight matrices (PWMs)](#motif-pwm)
- [3. Bining for Metageonome](#bining)
	- [3.1. CONCOCT](#bining-concoct)
- [4. Genome Assembly](#genome-assembly)
	- [4.1. 构建de Bruijn graph](#assembly-construct-dbg)
	- [4.2. 简化DBG得到assembly graph：filigree edges](#assembly-simplify-graph)
- [5. variants calling](#variants-calling)
	- [5.1. 背景知识](#introduction-to-variants-calling)
		- [5.1.1. base calling过程的错误及校正](#base-calling-error-and-correction)
			- [5.1.1.1. GATK-BQSR](#base-calling-error-gatk-bqsr)
		- [5.1.2. 比对过程的错误及校正](#mapping-error-and-correction)
	- [5.2. snp calling的数学原理](#snp-calling-mathmatic-principle)
		- [5.2.1. 李恒samtools/bcftools](#snp-calling-mathmatic-principle-from-li-heng)
		- [5.2.2. GATK](#snp-calling-mathmatic-principle-from-gatk)
- [6. k-mer based sequence analysis](#k-mer-based-sequence-analysis)
- [补充知识](#supplementary-knowledge)
	- [*1. 哈迪-温伯格平衡(Hardy-Weinberg equilibrium)法则](#hwe)




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

<a name="fast-alignment-based-on-suffix-tree-trie-tree"><h4>1.5.1. Trie树 [<sup>目录</sup>](#content)</h4></a>

<a name="fast-alignment-based-on-suffix-tree-trie-tree-1"><h5>1.5.1.1. Trie树及基本操作 [<sup>目录</sup>](#content)</h5></a>

（1） 定义

Trie 树，也叫“字典树”。顾名思义，它是一个树形结构。它是一种专门处理字符串匹配的数据结构，用来解决**在一组字符串集合中快速查找某个字符串的问题**

假设有 5 个字符串，它们分别是：code，cook，five，file，fat。现在需要在里面多次查找某个字符串是否存在。如果每次查找，都是拿要查找的字符串跟这 5 个字符串依次进行字符串匹配，那效率就比较低，有没有更高效的方法呢？

如果将这 5 个字符串组织成下图的结构，从肉眼上扫描过去感官上是不是比查找起来会更加迅速。

<p align="center"><img src=./picture/Algorithms-Bioinf-alignment-suffix-tree-trie-tree-1.jpg /></p>

通过上图，可以发现 Trie树 的三个特点：

- 根节点不包含字符，除根节点外每一个节点都只包含一个字符
- 从根节点到某一节点，路径上经过的字符连接起来，为该节点对应的字符串
- 每个节点的所有子节点包含的字符都不相同

（2）Trie树构造

在构造过程中的每一步，都相当于往 Trie 树中插入一个字符串。当所有字符串都插入完成之后，Trie 树就构造好了

<p align="center"><img src=./picture/Algorithms-Bioinf-alignment-suffix-tree-trie-tree-2.gif /></p>

（3）Trie树的插入操作

<p align="center"><img src=./picture/Algorithms-Bioinf-alignment-suffix-tree-trie-tree-3.gif /></p>

Trie树的插入操作很简单，其实就是将单词的每个字母逐一插入 Trie树。插入前先看字母对应的节点是否存在，存在则共享该节点，不存在则创建对应的节点

比如要插入新单词cook，就有下面几步：

- 插入第一个字母 c，发现 root 节点下方存在子节点 c，则共享节点 c
- 插入第二个字母 o，发现 c 节点下方存在子节点 o，则共享节点 o
- 插入第三个字母 o，发现 o 节点下方不存在子节点 o，则创建子节点 o
- 插入第三个字母 k，发现 o 节点下方不存在子节点 k，则创建子节点 k
- 至此，单词 cook 中所有字母已被插入 Trie树 中，然后设置节点 k 中的标志位，标记路径 root->c->o->o->k 这条路径上所有节点的字符可以组成一个单词cook

（4）Trie树的查询操作

在 Trie 树中查找一个字符串的时候，比如查找字符串 code，可以将要查找的字符串分割成单个的字符 c，o，d，e，然后从 Trie 树的根节点开始匹配

<p align="center"><img src=./picture/Algorithms-Bioinf-alignment-suffix-tree-trie-tree-4.jpg /></p>

如果要查找的是字符串cod(鳕鱼)呢？还是可以用上面同样的方法，从根节点开始，沿着某条路径来匹配，如图所示，绿色的路径，是字符串cod匹配的路径。但是，路径的最后一个节点「d」并不是橙色的，并不是单词标志位，所以cod字符串不存在。也就是说，cod是某个字符串的前缀子串，但并不能完全匹配任何字符串。

<p align="center"><img src=./picture/Algorithms-Bioinf-alignment-suffix-tree-trie-tree-5.jpg /></p>

（5）删除操作

Trie树的删除操作与二叉树的删除操作有类似的地方，需要考虑删除的节点所处的位置，这里分三种情况进行分析：

> 1. 删除整个单词（比如hi）
> 
> 	<p align="center"><img src=./picture/Algorithms-Bioinf-alignment-suffix-tree-trie-tree-6.gif /></p>
> 
> 	- 从根节点开始查找第一个字符h
> 	- 找到h子节点后，继续查找h的下一个子节点i
> 	- i是单词hi的标志位，将该标志位去掉
> 	- i节点是hi的叶子节点，将其删除
> 	- 删除后发现h节点为叶子节点，并且不是单词标志位，也将其删除
> 
> 	重复以上过程直到遇到内部节点
> 
> 	这样就完成了hi单词的删除操作
> 
> 2. 删除前缀单词（比如cod）
> 
> 	这种方式删除比较简单。
> 
> 	只需要将cod单词整个字符串查找完后，d节点因为不是叶子节点，只需将其单词标志去掉即可。
> 
> 	<p align="center"><img src=./picture/Algorithms-Bioinf-alignment-suffix-tree-trie-tree-7.gif /></p>
> 
> 3. 删除分支单词（比如cook）
> 
> 	<p align="center"><img src=./picture/Algorithms-Bioinf-alignment-suffix-tree-trie-tree-8.gif /></p>
> 
> 	与 **删除整个单词** 情况类似，区别点在于删除到 cook 的第一个 o 时，该节点为非叶子节点，停止删除，这样就完成cook字符串的删除操作

<a name="fast-alignment-based-on-suffix-tree-trie-tree-2"><h5>1.5.1.2. Trie树的应用 [<sup>目录</sup>](#content)</h5></a>

（1）前缀匹配

例如：找出一个字符串集合中所有以 五分钟 开头的字符串。我们只需要用所有字符串构造一个 trie树，然后输出以 五−>分−>钟 开头的路径上的关键字即可。

trie树前缀匹配常用于搜索提示。如当输入一个网址，可以自动搜索出可能的选择。当没有完全匹配的搜索结果，可以返回前缀最相似的可能。

<p align="center"><img src=./picture/Algorithms-Bioinf-alignment-suffix-tree-trie-tree-9.jpg /></p>

（2）字符串检索

给定一组字符串，查找某个字符串是否出现过，思路就是从根节点开始一个一个字符进行比较：

- 如果沿路比较，发现不同的字符，则表示该字符串在集合中不存在。

- 如果所有的字符全部比较完并且全部相同，还需判断最后一个节点的标志位（标记该节点是否代表一个关键字）。

<a name="fast-alignment-based-on-suffix-tree-suffix-trie-tree"><h4>1.5.2. Suffix Trie树 [<sup>目录</sup>](#content)</h4></a>

suffix trie树是一种特殊的trie树：

> - 索引对象：普通的trie树——给定的字符串集合来构建trie树；suffix trie树——给定的某一个字符串的所有后缀子字符串集合
> 
> - 节点字符数：普通的trie树，除了根节点外，任意节点只能包含一个字符；suffix trie树，除了根节点外，任意节点可以包含多个字符组成字符串

以字符串S=banana为例，来说明如何构建后缀树：

（1）首先列出S的所有后缀子字符串：

- banana$
- anana$
- nana$
- na$
- a$
- $

（2）根据上面列出的后缀子字符串，构建suffix trie（后缀字典树）

<p align="center"><img src=./picture/Algorithms-Bioinf-alignment-suffix-tree-1.png /></p>

（3）合并未分支路径，以节省存储空间

<p align="center"><img src=./picture/Algorithms-Bioinf-alignment-suffix-tree-2.png /></p>


<a name="subread-alignment-algorithm"><h3>1.6. Subread: seed-and-vote [<sup>目录</sup>](#content)</h3></a>


Yang Liao, Gordon K. Smyth, Wei Shi. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. Nucleic Acids Res. 2013 May 1;41(10):e108. 















<a name="motif"><h2>2. Motif Finding [<sup>目录</sup>](#content)</h2></a>

Motif Finding问题，实际上包含两个问题：

> （1）给定一系列的结合位点序列，找出一种可以表示这些结合位点序列组成的表示方式，成为一种序列模式（pattern）；
>
> （2）给出一些可能含有结合位点的序列，利用已经总结出来的这种结合位点的序列模式，将这些结合位点在给出的序列中找出来；

转录因子结合位点与限制性内切酶识别位点的差别：

> 限制性内切酶识别位点是相对简单，且确定的，比如像EcoRI的识别位点就可以简单地写成`GAATTC`，或者像HincII，它的识别位点有一些容错，为`GTYRAC`
>
> 限制性内切酶对识别位点要求很严格，只要有一个位点不对，它的识别效率就会显著下降
>
> 相比之下，转录因子的结合位点的可变性就高多了，比如$\lambda$操纵子的结合位点长度为8bp，其中只有两个位点是完全保守，其他位点都有不同程度的多样性

从上面的总结可以看出，限制性内切酶对位点的识别高度严格保守，导致其识别效率要么是0，要么是1

而转录因子则要求比较弱的序列模式即可，允许比较高的序列模式多样性，其识别和结合效率是一个连续变化的范围

这其实和它们的生物学功能是对应的

<a name="motif-em"><h3>2.1. MEME：EM算法 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/MeRIP-seq-meme-principle-1.png width=800/></p>

<p align="center"><img src=./picture/MeRIP-seq-meme-principle-2.png width=800/></p>

<p align="center"><img src=./picture/MeRIP-seq-meme-principle-3.png width=800/></p>

<a name="motif-pwm"><h3>2.2. position weight matrices (PWMs) [<sup>目录</sup>](#content)</h3></a>

$$
\text{TACGAT} \\
\text{TATAAT} \\
\text{TATAAT} \\
\text{GATACT} \\
\text{TATGAT} \\
\text{TATGTT}
$$

<p align="center"><img src=./picture/Algorithms-Bioinf-Motif-Finding-PWM-1.png width=400/></p>

某一个位点i的保守性：

$$2+\sum_{b=A}^T f_{b,i}\log_2 f_{b,i}$$

某一个位点i它的碱基组成为b的分值：

$$f_{b,i}\log_2 f_{b,i}$$

<a name="bining"><h2>3. Bining for Metageonome [<sup>目录</sup>](#content)</h2></a>

<a name="bining-concoct"><h3>3.1. CONCOCT [<sup>目录</sup>](#content)</h3></a>

结合序列组成特征 (sequence composition) 和跨样本覆盖度特征 (coverage across multiple samples) 进行bining

在进行bining之前需要将所有样本的reads进行混拼 (coassembly) 得到contigs

- Sequence composition features

以 k-mer 长度等于5为例

将相互之间成反向互补关系的 5-mers pairs 记做一种，则总共有512种 5-mers，对每一条contig计算其各自 5-mers 的组成频率从而构造出一个长度为v=512的向量 $Z_i$：

$$Z_i=(Z_{i,1},Z_{i,2},...,Z_{i,v})$$


为了保证每个 5-mers 频率的计数非零（为后面的对数转换做准备），进行伪计数处理，然后用该序列的 5-mers的总数进行标准化，得到新的向量 $Z_i'$：

$$Z_{i,j}'=\frac{Z_{i.j}+1}{\sum_{k=1}^v Z_{i,k}}$$

- Contig coverage features

用段序列比对软件，将各个样本（总共有M个样本）的reads比对到contigs上，计算每条contigs在每个样本中的 coverage (Mapped reads * read length / contig length)，得到表示 congtig i 的 coverage 的向量 $Y_i$：

$$Y_i=(Y_{i,1},Y_{i,2},...,Y_{i,M})$$


两轮标准化处理：

> 伪计数处理：额外添加一条比对到该 contig 上的 read；再用该**样本内**所有 contigs（contigs总数为N）的 coverage 进行标准化，得到新的向量 Y<sub>i</sub><sup>'</sup>：
>
> $$Y_{i,j}'=\frac{Y_{i,j}+\text{read|ength}/\text{contigLength}}{\sum_{k=1}^NY_{k,j}}$$
>
> 然后再在**contig内部**进行标准化，得到新的向量 Y<sub>i</sub><sup>''</sup>：
>
> $$Y_{i,j}''=\frac{Y_{i,j}'}{\sum_{k=1}^M Y_{i,k}'}$$

- Combine two features

将表示某一个 contig i 的序列组成特征的向量 Z<sub>i</sub><sup>'</sup> 和 coverage 特征向量 Y<sub>i</sub><sup>''</sup>合并成组成一个新的特征向量 X<sub>i</sub>（向量长度为E=V+M），同时进行对数转换：

$$X_i=\{\log(Z_i'),\log(Y_i'')\}$$

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

- **Phred 质量值与Base-Calling**

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

<a name="base-calling-error-gatk-bqsr"><h5>5.1.1.1. GATK-BQSR [<sup>目录</sup>](#content)</h5></a>

Phred碱基质量值是由测序仪内部自带的base-calling算法评估出来的，而这种base-calling算法由于受专利保护，掌握在测序仪生成商手中，研究人员并不能了解这个算法的细节，它对于人们来说就是一个黑盒子

而测序仪的base-calling算法给出的质量评估并不十分准确，它带有一定程度的系统误差（非随机误差），使得实际测序质量值要么被低估，要么被高估

BQSR试图利用机器学习的方法来对原始的测序质量值进行校正

例如：

> 对于一个给定的Run，我们发现，无论什么时候我在测序一个AA 的子序列时，改子序列后紧接着的一个任意碱基的测序错误率总是要比它的实际错误率高出1%，那么我就可以将这样的碱基找出来，将它的原始测序错误率减去1%来对它进行校正

会影响测序质量评估准确性的因素有很多，主要包括序列组成、碱基在read中的位置、测序反应的cycle等等，它们以类似于叠加的形式协同产生影响，这些可能的影响因素被称作协变量 (covariable)

注意：BQSR只校正碱基质量值而不改变碱基组成，特别是对于那些质量值偏低的碱基，我们只能说它被解析成当前碱基组成的准确性很低，但是我们又无法说明它实际更可能是哪种碱基，所以干脆不改

那么，BQSR的工作原理是怎样的？

BQSR本质上是一种回归模型

前提假设：影响质量评估的因素只有reads group来源，测序的cycle和当前测序碱基的序列组成背景（这里将它上游的若干个连续位点的碱基组成看作它的背景，一般为2~6，BQSR中默认为6）

则基于这个前提假设，我们可以得出以下结论：

> 相同reads group来源，同处于一个cycle，且序列背景相同的碱基，它们具有相同的测序错误率，这样的碱基组成一个bin

则可以建立这样的拟合模型：

$$X_i=(RG_i,Cyc_i,Context_i) \quad \begin{matrix} f \\ \to \end{matrix} \quad y_i$$

其中，i表示当前碱基，$RG_i$表示碱基所属的Reads Group来源，$Cyc_i$表示该碱基所在的测序cycle，$Context_i$表示该碱基的序列组成背景，$y_i$表示该碱基的实际测序质量(emprical quality)

这三个分量可以直接通过输入的BAM文件的记录获得，那如何获得实际的实际测序质量呢？

可以通过BAM文件中的比对结果推出

用给定的大型基因组测序计划得到的人群变异位点作为输入，将样本中潜在变异位点与人群注释位点overlap的部分过滤掉，则剩下的那些位点，我们假设它们都是“假”的变异位点，是测序错误导致的误检

则实际测序质量为：

$$EQ=-10\log \frac{\sharp mismatch + 1}{\sharp bases + 2}$$

注意：emprical quality是以bin为单位计算出来的

这样，有了X和Y，就可以进行拟合模型的训练了，训练好的模型就可以用于碱基质量值的校正

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

<a name="snp-calling-mathmatic-principle"><h3>5.2. snp calling的数学原理 [<sup>目录</sup>](#content)</h3></a>

<a name="snp-calling-mathmatic-principle-from-li-heng"><h4>5.2.1. 李恒samtools/bcftools [<sup>目录</sup>](#content)</h4></a>

基本思路：

> 先要对每个候选的变异位点，计算落在这个位点上的reads的allele的似然，然后基于这多条reads的allele likelihood，推断出出现哪种genotype的概率最大

前提假设：

> - **数据量与数据质量要求**：genotype置信度的推断依据与数据量和数据质量，低覆盖度或低质量的数据都会导致较低的置信度。
>
> 	只采用高测序质量和高比对质量的结果进行后续的分析
>
> - **二倍型 (Ploidy)**：默认提供的物种都是二倍型的，倍型的设定会影响到后续基因分型算法(genotyping algorithmn)的公式推导
>
> - **使用paired end reads**
>
> - **Single-sample vs multi-sample**
>
> 	对于单样本和单样本的情况，使用的genotying model是不一样的

采用李恒在2011年的文章提出的方法计算


<a name="snp-calling-mathmatic-principle-from-gatk"><h4>5.2.2. GATK [<sup>目录</sup>](#content)</h4></a>

GATK进行SNP calling的核心算法为HaplotypeCaller，这个也是GATK中最核心的算法，理解了这个算法基本上就明白了GATK变异检测的原理

HaplotypeCaller它本质上是对贝叶斯原理的应用，只是相于同类算法它有点不同之处

算法思想概述：

> HaplotypeCaller首先是根据所测的数据，先构建这个群体中的单倍体组合（我认为这也是Haplotype这个名字的由来），由于群体中的单倍体是有多个的，所以最好是多个人一起进行HaplotypeCaller这样构建出来的单倍体组合会越接近真实情况
>
> 构建出单倍体的组合之后（每一个单倍体都有一个依据数据得出的后验概率值），再用每个样本的实际数据去反算它们自己属于各个单倍体组合的后验概率，这个组合一旦计算出来了，对应位点上的碱基型（或者说是基因型，genotype）也就跟着计算出来了：
>
> 计算每一个后候选变异位置上的基因型（Genotype）后验概率，最后留下基因型（Genotype）中后验概率最高的哪一个

下面进行详细地说明：

在HaplotypeCaller中变异检测过程被分为以下四个大的步骤

![](./picture/Algorithms-Bioinf-variants-calling-algorithmn-GATK-1.png)

（1）确定候选变异区域（ActiveRegion）

通过read在参考基因组上的比对情况，筛选出潜在的变异区域，这些区域在GATK中被称为ActiveRegion

（2）通过对候选变异区域进行重新组装来确定单倍型

对于每个ActiveRegion，GATK会利用**比对到该区域上的所有read**（如果有多个样本那么是所有这些样本的reads而不是单样本进行）构建一个类似于de Bruijn的图对ActiveRegion进行局部重新组装，构建出该区域中可能的单倍型序列。然后，使用Smith-Waterman算法将每个单倍型序列和参考基因组进行重新比对，重新检测出潜在的变异位点

（3）依据所给定的read比对数据计算各个单倍型的似然值

在步骤2的基础上，我们就得到了在ActiveRegion中所有可能的单倍型序列，接下来需要评估现有数据中对这些单倍型的支持情况

GATK使用PairHMM算法把原本比对于该区域中的每一条read依次和这些单倍型序列进行两两比对，这样我们就可以得出一个read-单倍型序列成对的似然值矩阵，例如以单倍型为列，以read为行，将矩阵记作$(a_{i,j})$

$$
\begin{array}{l|c|c|c|c}
0 & H_1 & H_2 & .. & H_m \newline
\hline
r_1 & a_{11} & a_{12} & .. & a_{1m} \newline
r_2 & a_{21} & a_{22} & .. & a_{2m} \newline
.. & .. & .. & .. & .. \newline
r_n & a_{n1} & a_{n2} & .. & a_{nm} \newline
\hline
\end{array}
$$

则矩阵中的某一个元素$a_{ij}$表示在read i支持单体型为$H_j$的似然

这个似然值矩阵很重要，因为在获得这个矩阵之后，GATK会在每一个潜在的变异位点上把这些似然值相加合并，计算等位基因的边缘概率，这个边缘概率实际上是每一个read在该位点上支持其为变异的似然值

（* 在该步骤中，Pair-HMM这实际上是GATK中最为耗费计算资源的那部分了，GATK的加速也是常常以此为突破口——比如GPU加速或者把Pair-HMM模块烧录到FPGA芯片中，也有人从算法本身出发发表了关于如何更快计算Pair-HMM的文章：`https://journals.sagepub.com/doi/pdf/10.1177/1176934318760543` ）

（4）计算每一个样本在最佳单倍型组合下的基因型（Genotype）

在完成了步骤3之后，我们就知道了**每一条read在每个候选变体位点上支持每一种等位基因（Allele）的概率**了。那么，最后要做的就是通过这些似然值，计算出候选变异位点上最可能的样本基因型，也就是Genotype——这也是发现真正变异的过程。这就需要应用贝叶斯原理来完成这个计算了——GATK这也是到这一步才使用了该原理，通过计算就可以得到每一种Genotype的可能性，最后选择后验概率最高的那一个Genotype作为结果输出至VCF中

后面的分析中，对于每一个变异位点假设只有二等位形式——注意：这和一个ActiveRegion中存在多种单体型不矛盾，若一个ActiveRegion在群体中存在n个变异位点，在只考虑二等位形式的前提下，该区域具有的单体型总共有$2^n$种

下面来推导某个样本中的某一个变异位点最可能的SNP形式

该样本在该位点的genotype为G的后验概率为：

$$P(G \mid D) = \frac{P(G)P(D \mid G)}{\sum_i P(G_i)P(D \mid G_i)} \tag{1}$$

由于分母部分对于任何形式genotype都一样，即它是个定值，所以可以忽略，因此上面的公式可以简化成：

$$P(G \mid D) = P(G)P(D \mid G) \tag{2}$$

其中，$P(G)$为genotype为G的先验概率，理论上为样本来源的群体中allele为G的频率，这个一般需要前期给定，若不给定的话，GATK会默认每种G的频率均等

$P(D \mid G)$表示在已知样本genotype为G的前提下，对样本进行测序得到的测序数据为D（仅考虑该ActiveRegion范围内的）的条件概率，我们假设每条reads之间是相互独立的，所以

$$P(D \mid G)=\prod_j P(D_j \mid G)\tag{3}$$

其中，$D_j$表示该样本测序数据D中的第j条read

由于我们正常人都是二倍体，则对于某一条reads，它既可能来自于同源染色体1，记作$H_1$，也可能开自于同源染色体2，记作$H_2$，所以

$$
\begin{aligned}
&\quad P(D_j \mid G) \newline
&= P(D_j,H_1 \mid G) + P(D_j,H_2 \mid G) \newline
&= P(H_1 \mid G)P(D_j \mid H_1) + P(H_2 \mid G)P(D_j \mid H_2)
\end{aligned} \tag{4}
$$

由于理论上一条read来源于$H_1$还是$H_2$的概率是均等的，都为1/2，即$P(H_1 \mid G)=P(H_2 \mid G)=1/2$，所以

$$P(D_j \mid G)=\frac{P(D_j \mid H_1)}{2} + \frac{P(D_j \mid H_2)}{2} \tag{5}$$

因此(3)可以改写成

$$P(D \mid G)=\prod_j \left( \frac{P(D_j \mid H_1)}{2} + \frac{P(D_j \mid H_2)}{2}\right) \tag{6}$$

现在如果想算出$P(G \mid D)$，就差$P(D_j \mid H_n)$了，那么，如何算$P(D_j \mid H_n)$呢？

上面已经提到，$P(D_j \mid H_n)$表示的是由同源染色体$H_n$产生read $D_j$的条件概率，而每条同源染色体有它各自的单体型，所以这里可以把$H_n$理解为它对应的单体型，则$P(D_j \mid H_n)$可以理解为在特定单体型$H_n$的前提下，产生read $D_j$的条件概率

<a name="snp-calling-using-samtools-bcftools"><h4>5.2.1. samtools/bcftools [<sup>目录</sup>](#content)</h4></a>

（1）李恒对千人基因组计划的思考：

- 为什么选择低覆盖多样本的策略而非高覆盖度少样本的策略？

	如果采用的是高覆盖的的测序策略，那么基本上可以利用reads之间的相互校验来得到更为准确的信息

	但是如果采用的是低覆盖多样本的策略，我们可以减少采样的波动/不稳定性（高覆盖意味着样本量就比较少，采样就很可能不具有代表性，即意味着采样的波动/不稳定性），可以鉴定出存在多个样本中的variants，且能得到许多在群体中罕见的variants

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
| $m_i$	| Ploidy of the $i$-th sample ($1≤i≤n$)，即i样本的倍型 |
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

- 估计某个样本出现特定基因型g的概率：$L(g)$

	对于某一个样本的某一个位点，有$k$条 reads 比对上，其中有$l$条 ($0 \le l\le k$) 序列在该位点的碱基组成与 reference 一致，剩余 $k-l$ 条与 reference 不同，其中第 $j$ 条上该碱基的测序错误率为 $\epsilon_j$，则该样本的基因型与ref一致的有 $g \in [0,m]$ 种的概率为

	$$
	\begin{aligned}
	&\quad L(g) \newline
	&= Pr(d \mid g) \newline
	&= \prod_{i=1}^l Pr_i(A)\prod_{j=l+1}^k Pr_j(\overline A)  & (1)\newline
	&= \prod_{i=1}^l [Pr_i(B , A)+Pr_i(\overline B , A)]\prod_{j=l+1}^k [Pr_j(B , \overline A)+Pr_j(\overline B , \overline A)] & (2)\newline
	&= \prod_{i=1}^l [Pr_i(B,C) + Pr_i(\overline B,\overline C)] \prod_{j=l+1}^k [Pr_j(B,\overline C) + Pr_j(\overline B,C)] & (3)\newline
	\end{aligned}
	$$

	> 其中，$m$ 是该物种的倍性，普通人是二倍体，因此一般 $m=2$
	>
	> 事件$A=\{测序碱基与\text{ref}一致\}$，则$\overline A=\{测序碱基与\text{ref}不一致\}$
	>
	> 事件$B=\{该碱基的测序是正确的\}$，则$\overline B=\{该碱基的测序是错误的\}$
	>
	> 事件$C=\{实际碱基与\text{ref}一致\}$，则$\overline C=\{实际碱基与\text{ref}不一致\}$

	上面公式中，从(2)到(3)的推导涉及到最基本的逻辑常识，这里就不再赘述了

	由于测序错误与基因组的组成无关，即$B \bot C$，因此上面的公式可以向下继续推导：

	$$
	\begin{aligned}
	&=  \prod_{i=1}^l [Pr_i(C)Pr_i(B) + Pr_i(\overline C)Pr_i(\overline B)] \prod_{j=l+1}^k [Pr_j(\overline C)Pr_j(B) + Pr_j(C)Pr_j(\overline B)] & (4)\newline
	&= \prod_{i=1}^l \left[ \frac{g}{m}(1-\epsilon_i) + \frac{m-g}{m}\epsilon_i \right] \prod_{j=l+1}^k \left[  \frac{m-g}{m}(1-\epsilon_j) +  \frac{g}{m}\epsilon_j\right] & (5)\newline
	&= \frac{1}{m^k}\prod_{i=1}^l [g(1-\epsilon_i) + (m-g)\epsilon_i] \prod_{j=l+1}^k [(m-g)(1-\epsilon_j) + g\epsilon_j] & (6)
	\end{aligned}
	$$

	上面公式中，(4)到(5)的推导利用了：

	$$
	\begin{aligned}
	&Pr(B)=1-\epsilon, \quad Pr(\overline B)=\epsilon & (7)\newline
	&Pr(C)=\frac gm , \quad Pr(\overline C)=\frac{m-g}{m} & (8)
	\end{aligned}
	$$

	\(7\)公式很好理解，在这里就不作更多的解释

	对公式(8)，下面作一下简单的解释：

	> 由于上面的前提假设中就已经提到，只考虑双等位情况，ref allele即是双等位中的一种，则对于一个m倍体的个体，它该等位基因座上有m个等位基因，其中与ref allele一致的有g个，则剩下m-g个基因座上的allele与ref allele不一致
	>
	> 则，随机从这m个基因座中抽一个，其基因型与ref一致的概率为$Pr(C)=g/m$，与ref不一致的概率为$Pr(\overline C)=1-g/m=(m-g)/m$

- 估计某一个位点的allele frequency：$\psi$

	假设某一个位点与ref一致的allele的频率为$\psi$，从多个样本的测序数据中估计出这个频率的数值，基本思想为：

	> 推导出$Pr(D\mid \psi)$的表达式
	>
	> 然后基于极大似然估计（如下）来求出$\psi^*$
	>
	> $$\psi^* = \arg \max_{\psi} Pr(D\mid \psi)$$

	对于第i个样本，它的倍型为$m_i$，基因型为$g_i$，测序数据为$d_i$，根据HWE（Hardy-Weinberg，哈迪-温伯格法则，简称哈温平衡）：

	> 哈迪-温伯格(Hardy-Weinberg)法则
	>
	> 核心思想：一个不发生突变、迁移和选择的无限大的相互交配的群体中，基因频率和基因型频率将逐代保持不变
	>
	> 哈迪-温伯格定律可分为3个部分：
	>
	> （1）在一个无穷大的随机交配的群体中，没有进化的压力（突变、迁移和自然选择）；
	> （2）基因频率逐代不变；
	>
	> （3）随机交配一代以后基因型频率将保持平衡：
	>
	>  $p^2$表示AA的基因型的频率，$2pq$表示Aa基因型的频率$q^2$表示aa基因型的频率。其中p是A基因的频率；q是a基因的频率。基因型频率之和应等于1，即$p^2+ 2pq + q^2 = 1$

	在该位点与ref一致的allele的频率为$\psi$情况下，出现我们的测序数据D的概率，即$\psi$的似然为：

	$$
	\begin{aligned}
	&\quad L(\psi)\newline
	&= Pr(D \mid \psi ) \newline
	&=\prod_{i=1}^n Pr_i(d_i \mid \psi) & (1)\newline
	&=\prod_{i=1}^n \sum_{g=0}^{m_i} Pr_i(d_i,g \mid \psi) & (2)\newline
	&= \prod_{i=1}^n \sum_{g=0}^{m_i} Pr_i(d_i \mid g)Pr_i(g\mid \psi) & (3)\newline
	&= \prod_{i=1}^n \sum_{g=0}^{m_i} L_i(g)Binomial(g,m_i,\psi) & (4)\newline
	&= \prod_{i=1}^n \sum_{g=0}^{m_i} L_i(g) \left(\begin{matrix} m_i \newline g\end{matrix}\right) \psi^g(1-\psi)^{m_i-g} & (5)
	\end{aligned}
	$$

	上面公式中第（3）步到第（4）步的推导，利用了$Pr(g \mid \psi)=Binomial(g,m,\psi)$，即对于某一个$m$倍体样本，在群体的ref allele频率为$\psi$的情况下，它与ref allele一致的allele数为g的概率，这相当于在进行m重伯努利实验：

	> 由于只考虑双等位型，则可以把ref allele 记为A，不一致的allel记为a，且$Pr(A)=\psi,Pr(a)=1-\psi$，则对于一个m倍体的个体，其ref allele的数量为g，则它的基因型为$A_1A_2...A_ga_{g+1}a_{g+2}...a_{m}$
	>
	> 在哈温平衡的理想群体中，这样基因型的个体可以看作是由m次伯努实验得到g次结果为1的事件得到的，则出现这样的事件的概率就为
	>
	> $$Binomial(g,m,\psi)=\left(\begin{matrix} m \\ g\end{matrix}\right) \psi^g(1-\psi)^{m-g}$$

	下面可以通过极大似然估计来得到ref allele的频率$\psi$：

	$$\psi^* = \arg \max_{\psi} L(\psi)$$

	由于

	$$L(\psi)=\prod_{i=1}^n \sum_{g=0}^{m_i} L_i(g) \left(\begin{matrix} m_i \\ g\end{matrix}\right) \psi^g(1-\psi)^{m_i-g}$$

	其中含有隐变量g，所以我们要进行的是含有隐变量的极大似然估计

	对于含有隐变量的极大似然估计，首选的最优化方法一般是EM算法，其中第t步到第t+1步的迭代关系式为：

	$$\psi^{(t+1)}=\frac 1M \sum_{i=1}^n \frac{\sum_g gL_i(g)Binomial(g,m_i,\psi^{(t)})}{\sum_g L_i(g)Binomial(g,m_i,\psi^{(t)})}$$

- 估计non-reference alleles的数量

	首先定义一个名词：

	> site reference allele count：在一个位点上与ref allele一致的allele的数量，则对于样本i来说，它的site reference allele count可以记作$G_i$

	定义向量$\vec G=(G_1,G_2,...,G_n)$，表示所有样本在某一个位点的site reference allele count

	$X=\sum_i G_i$表示所有样本的site reference allele count的总和

	$$Pr(\vec G = \vec g \mid X=k)= \mathbb{I}(k=\sum_{i=1}^n g_i) \frac{\prod_{i=1}^n \left( \begin{matrix} m_i \\ g_i\end{matrix} \right)}{\left( \begin{matrix} M \\ k\end{matrix} \right)}$$

	其中，$\mathbb{I}(k=\sum_{i=1}^n g_i)$中的$\mathbb{I}(·)$表示布尔函数，即$\mathbb{I}(True)=1 \, \text{or} \, \mathbb{I}(False)=0$，也可以用克罗内克函数(Kronecker delta function) $\delta_{··}$ 表示，不过为了避免引入不必要的复杂概念对后续推导过程的干扰，这里还是使用常规的布尔函数表示

	则群体中site reference allele count 为k似然为：

	$$
	\begin{aligned}
	&\quad L(k) \newline
	&=Pr(D \mid X=k) & (1)\newline
	&= \sum_{\vec g} Pr(D,\vec G = \vec g \mid X=k)  & (2)\newline
	&= \sum_{\vec g} Pr(G = \vec g \mid X=k) Pr(D \mid \vec G = \vec g) & (3)\newline
	&= \sum_{\vec g} \left( \mathbb{I}(k=\sum_{i=1}^n g_i) \frac{\prod_{i=1}^n \left( \begin{matrix} m_i \newline g_i\end{matrix} \right)}{\left( \begin{matrix} M \\ k\end{matrix} \right)} \right)\left(\prod_{i=1}^n L_i(g_i)\right) & (4)\newline
	&= \frac{1}{\left( \begin{matrix} M \\ k\end{matrix} \right)} \sum_{\vec g} \left(\mathbb{I}(k=\sum_{i=1}^n g_i)\prod_{i=1}^n \left( \begin{matrix} m_i \newline g_i\end{matrix} \right) L_i(g_i)\right) & (5)
	\end{aligned}
	$$

	其中，$k$要满足$0 \le k \le M=\sum_{i=1}^n m_i$，即群体的site reference allele count来自于群体的总染色体数M，不能超过M，否则这就是一个不可能事件，则$L(k)=0$

	为了方便后续的推导，定义一个变量：

	$$z_{jl}=\sum_{\vec g = (g_1,...,g_j)} \left(\mathbb{I}(l=\sum_{i=1}^j g_i)\prod_{i=1}^j \left( \begin{matrix} m_i \\ g_i\end{matrix} \right) L_i(g_i)\right)$$

	其中$j\le n$，即上式只考虑前$j$个样本的情况，相当于上面公式(5)的右半部分

	假设$j$个样本的最后一个，即第$j$个样本的为$m_j$倍体，rel allele数为$g_j$，则扣除最后一个样本的值为$Z_{j-1,l-g_j}$，那么$Z_{j,l}$与$Z_{j-1,l-g_j}$之间的有什么关系呢？

	先写出$Z_{j-1,l-g_j}$的公式：

	$$Z_{j-1,l-g_j}=\sum_{\vec g = (g_1,...,g_{j-1})} \left(\mathbb{I}(l-g_j=\sum_{i=1}^{j-1} g_i)\prod_{i=1}^{j-1} \left( \begin{matrix} m_i \\ g_i\end{matrix} \right) L_i(g_i)\right)$$

	可以看出：

	$$Z_{jl}=\sum_{g_j=0}^{m_j} Z_{j-1,l-g_j} \left( \begin{matrix} m_j \\ g_j\end{matrix} \right) L_j(g_j)$$

	根据上面的公式的特点，可以采用迭代的方法最终算出$Z_{nk}$，则可以进一步算出$L(k)$：

	$$L(k)=\frac{Z_{nk}}{\left( \begin{matrix} M \\ k\end{matrix} \right)}$$

<a name="k-mer-based-sequence-analysis"><h2>6. k-mer based sequence analysis [<sup>目录</sup>](#content)</h2></a>

<a name="k-mer-based-sequence-analysis"><h3>6.1. SEEKR: noncoding RNA similarity [<sup>目录</sup>](#content)</h3></a>

方法思想：

<p align='center'><img src=./picture/Algorithms-Bioinf-k-mer-SEEKR-1.png width=600/></p>

方法流程：

<p align='center'><img src=./picture/Algorithms-Bioinf-k-mer-SEEKR-2.png width=600/></p>

性能表现：在较低信噪比下，仍能给出相对准确的预测

<p align='center'><img src=./picture/Algorithms-Bioinf-k-mer-SEEKR-3.png width=400/></p>

<a name="supplementary-knowledge"><h2>补充知识 [<sup>目录</sup>](#content)</h2></a>

<a name="hwe"><h3>*1. 哈迪-温伯格平衡(Hardy-Weinberg equilibrium)法则 [<sup>目录</sup>](#content)</h3></a>

哈温平衡针对的是理想群体，即的前提假设：

> - 群体无穷大
> - 随机交配
> - 没有进化的力量

分析讨论哈温平衡的三个假设：

**（1）群体无穷大**

若一个群体的大小有限，可能导致基因频率和预期的比例随机发生偏差。这种基因频率的改变就称遗传漂变（genetic draft）

所谓的无穷大完全是设想的。没有任何群体具有无穷的个体。然而样本的误差仅对一个相当小群体的基因频率有明显的影响。实际应用时群体不需无穷大，只要不至太小即可

**（2）随机交配**

随机交配（random mating）是指各基因型之间的交配和群体中这些基因型的频率成正比。更为特别的是两个基因型之间交配的概率等于两个基因型频率的乘积

为了说明随机交配，现以人类的M-N血型为例来解释：

> M-N血型是由于在细胞的表面上存在一种抗原，与ABO系统的抗原相似。但在M-N系统中除了产生不相容以外，在输血时并不会产生凝血。M-N血型是由带有两共显性等位基因$L^M$和$L^N$的座位决定的。在爱斯基摩人的群体中，3种M-N基因型频率分别为$L^M/L^M$= 0.835, $L^M/L^N$ =0.156，$L^N/L^N$= 0.009。若爱期基摩人的婚配是随机的，那么$L^M/L^M$男人和$L^M/L^M$女人婚配的率就应等于$L^M/L^M$基因型频率乘以$L^M/L^M$的频率，即0.835×0.835= 0.679。其它基因型之间的婚配率的计算也与此相似

注意：**随机交配是针对所有性状**

> 若随机交配是针对所有性状，那么人类群体就不能符合哈迪-温伯格定律的要求了：人类择偶并不是随机的，而是对智商、外貌、性格、身高、肤色、学历以及社会地位等都有一定的要求
>
> 虽然对某些性状的要求不是随机的，但大部分人对血型等并无要求，甚至有的人并不知道自己的M-N系统的具体血型。因此哈迪-温伯格定律要求的随机性是指诸如像血型这样一些性状，而不是那样非随机性状的座位

**（3）没有进化的力量**

在哈迪-温伯格定律中人们只关注遗传是否会改变基因频率以及繁殖怎样会影响到基因型频率。因此其它的进化力量可被排除

在没有进化力量作用在群体上时，只适用于某些座位，其它的座位可能照样受到进化力量的影响


---

参考资料：

(1) Benjamin Buchfink, Chao Xie & Daniel H. Huson, Fast and Sensitive Protein Alignment using DIAMOND, Nature Methods, 12, 59–60 (2015) doi:10.1038/nmeth.3176.

(2) Chaisson M J, Tesler G. Mapping single molecule sequencing reads using basic local alignment with successive refinement (BLASR): application and theory[J]. BMC bioinformatics, 2012, 13(1): 238.

(3) [生信算法《三代测序序列比对利器-BLASR，更小更快更方便 》](https://mp.weixin.qq.com/s/7xZVvyShcwjPGbgKiLniag)

(4) Stormo, G. DNA binding sites: representation and discovery. Bioinformatics 16, 16–23 (2000).

(5) Alneberg, J. et al. Binning metagenomic contigs by coverage and composition. Nat. Methods 11, 1144–1146 (2014).

(6) Beaulaurier J, Zhu S, Deikus G, et al. Metagenomic binning and association of plasmids with bacterial host genomes using DNA methylation.[J]. Nature Biotechnology, 2017, 36(1).

(7)  Nurk S., Meleshko D., Korobeynikov A., Pevzner P. A. metaSPAdes: a new versatile de novo metagenomics assembler.	Genome Research, 2017 

(8) Nielsen R, Paul JS, Albrechtsen A, Song YS. Genotype and SNP calling from next-generation sequencing data. Nat Rev Genet. 2011;12(6):443–451.

(9) Li H. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 2011;27(21):2987–2993.

(10) [Samtools math notes](http://www.broadinstitute.org/gatk/media/docs/Samtools.pdf)

(11) 黄树嘉·知识星球《达尔文生信星球》

(12) [GATK官方文档《Methods and Algorithms: HaplotypeCaller in a nutshell》](https://software.broadinstitute.org/gatk/documentation/article?id=11068)

(13) [GATK官方文档《Methods and Algorithms: Assigning per-sample genotypes (HaplotypeCaller)》](https://software.broadinstitute.org/gatk/documentation/article?id=11079)

(14) Kirk JM, Kim SO, Inoue K, et al. Functional classification of long non-coding RNAs by k-mer content. Nat Genet. 2018 Oct;50(10):1474-1482.


(15) [CSDN · Marphy11《哈迪-温伯格平衡(Hardy-Weinberg equilibrium)法则》](https://blog.csdn.net/lj695242104/article/details/41014339)

(16) [【 五分钟学算法】看动画轻松理解「Trie树」](https://mp.weixin.qq.com/s/Y5_r4C5a9gU0FDtqXD9bkQ)
