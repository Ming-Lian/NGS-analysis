<a name="content">目录</a>

[连锁不平衡：linkage disequilibrium](#title)
- [什么是连锁不平衡](#introduction)
- [定量连锁不平衡](#quantify-ld)
	- [r2 和D'的比较](#compare-r2-with-d-plus)
- [可视化连锁不平衡](#visualize-ld)
	- [热图可视化haplotype block](#heatmap-for-haplotype-block)
	- [LD衰减曲线](#ld-decay-curve)





<h1 name="title">连锁不平衡：linkage disequilibrium</h1>

<a name="introduction"><h2>什么是连锁不平衡 [<sup>目录</sup>](#content)</h2></a>

连锁不平衡：指的是在某一群体中，两个基因同时遗传的频率大于随机组合的频率

举个例子：

> 基因A的两个allel 分别用A和a 表示，基因B的allel分别用B和b表示，如果这这两个基因完全独立遗传，也就是说其allel 完全随机组合，那么后代中会出现4种单倍型，AB， Ab, aB， ab, 而且出现的概率都是相同的，都是0.25；如果这两个基因在遗传时不是独立的，意味着后代中单倍型出现的概率不是完全随机的了，我们就可以说两个基因是存在连锁关系的，基因在遗传时出现连锁的现象就叫做连锁不平衡。

那么，如何定量描述连锁不平衡？

<a name="quantify-ld"><h2>定量连锁不平衡 [<sup>目录</sup>](#content)</h2></a>

1. 用D值衡量：

	<p align="center">D = P(AB) -  P(A) X P(B)</p>

	独立遗传时，单倍型AB出现的概率为 **P(A) X P(B)**, 这个概率我们暂且称之为**理论概率**；当出现了连锁不平衡时，单倍型AB出现的概率用**P(AB)**表示，我们暂且称之为**实际概率**；这两个概率之间的差，就反应了连锁不平衡的程度
	
	D值不等于0，就可以说两个基因之间是连锁不平衡的，D绝对值大小直接反应了两个基因之间的连锁程度的大小，绝对值越大，连锁程度越大

	不足：其严格依赖于等位基因频率（allele frequency），故不适合应用于表述实际的LD强度，尤其是进行不同研究的LD值的相互比较

2. 用`D'`衡量：

	<p align="center">D’ = D / Dmax</p>

	Dmax 的计算方式如下：

	<p align="center"><img src=./picture/LD-Dmax.png width=600 /></p>

	D'值可以看做是归一化之后的D值，归一化之的值可以用于比较不同基因连锁程度的大小。D’的取值范围为0到1，`D’ = 0` 表示完全连锁平衡，独立遗传；`D’ = 1` 表示完全连锁不平衡

3. 用r衡量：

	<p align="center"><img src=./picture/LD-r-formula.png width=400 /></p>

	通常情况下，会通过r值的平方来表征连锁不平衡程度，r平方等于0时，表示完全连锁平衡，独立遗传；r平方等于1时, 表示完全连锁不平衡

在实际分析中，我们通常会拿到样本的基因分型文件，通过这个文件我们可以非常容易的计算出allel的频率，但是对于单倍型的频率是不能直接计算得到的，都是借助算法估算出单倍型的概率，然后进行计算

即，我们一般能拿到的是VCF格式的基因分型文件，得到的只能是一个个离散的位点，而无法得到单体型结果，此时只能进行allele频率的计算，也就无法计算单体型的频率，但是可以通过单体型分型(phasing)获得推断出来的单体型从而来估计单体型频率

<a name="compare-r2-with-d-plus"><h3>r2 和D'的比较 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/LD-relation-Dplus-r2.jpg width=800 /></p>

r2和D'反映了LD的不同方面。r2包括了重组和突变，而D'只包括重组史

D'能更准确地估测重组差异，但样本较小时，低频率等位基因组合可能无法观测到，导致LD强度被高估，所以D'不适合小样本群体研究；

LD衰减作图中通常采用r2来表示群体的LD水平；Haplotype Block中通常采用D'来定义Block；

<a name="visualize-ld"><h2>可视化连锁不平衡 [<sup>目录</sup>](#content)</h2></a>

理论上来说任意两个基因之间都可能存在连锁不平衡，但是实际操作中，认为只有一定区间范围内的基因会存在连锁不平衡，距离大于区间的基因，两者出现连锁不平衡的概率非常小，所以就不去计算。

<a name="heatmap-for-haplotype-block"><h3>热图可视化haplotype block [<sup>目录</sup>](#content)</h3></a>

对于连锁不平衡的结果，通常采用heatmap热图的形式进行展示，`haploview` 给出的LD heatmap 示例如下：
	
<p align="center"><img src=./picture/LD-visulize-LD-heatmap.jpg width=800 /></p>
	
颜色从白色到红色，代表连锁程度从低到高，方框中的数值为r2，为了美观，这里将 r2 乘以了100

**haplotype block**，即单体型块，即连锁不平衡区域，是指同一条染色体上处于连锁不平衡状态的一段连续的区域。单体型块分析可以用于筛选tag SNP、确定候选基因的范围等

<p align="center"><img src=./picture/LD-haplotype-block-1.jpg width=600 /></p>

<p align="center">Nature Genetics 48, 927–934 (2016) doi:10.1038/ng.3596</p>

<p align="center"><img src=./picture/LD-haplotype-block-2.jpg width=600 /></p>

<p align="center">Plant Biotechnol J. 2017 Mar 29. doi: 10.1111/pbi.12734</p>

<a name="ld-decay-curve"><h3>LD衰减曲线 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/LD-visulize-LD-curve.jpg width=600 /></p>

横坐标为基因之间的距离，纵坐标为衡量连锁不平衡的R2值。从图中可以看出连锁不平衡的规律，在一定距离内存在连锁不平衡程度较高，大于一定距离后，出现连锁不平衡的概率就大大降低了。这就是为什么在实际操作中只计算一定范围内的连锁不平衡的原因。

通常会使用1个标准——**“LD衰减距离”**来描述LD衰减速度的快慢。

<p align="center"><img src=./picture/LD-LD-Decay-distance.jpg width=800 /></p>

<p align="center">Nature Biotechnology 30, 105–111 (2012)  doi:10.1038/nbt.2050</p>

**LD衰减距离**：当平均LD系数r2 衰减到一定大小的时候，对应的物理距离。“一定大小”是这个定义的关键点，但没有特别统一的标准，在不同文章中标准不同。常见的标准包括：

> a）LD系数降低到最大值的一半；
> 
> b）LD系数降低到0.5以下；
> 
> c）LD系数降低到0.1以下；
> 
> d）LD系数降低到基线水平（注意，不同物种的基线值是不同的）;

下次你在文章中看到“LDdecay distance is XXkb”的时候，不要忘了看看文章使用的标准是什么




---

参考资料：

(1) [【简书】连锁不平衡：linkage disequilibrium](https://www.jianshu.com/p/477bdb7f57fa)

(2) [【百迈客基因】群体遗传分析—LD连锁不平衡](https://mp.weixin.qq.com/s/1y27LMskwKXx_PfPSNr1hw)
