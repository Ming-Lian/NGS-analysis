<a name="content">目录</a>

[群体遗传学知识点](#title)
- [1. Fst及计算方法](#fst)
- [2. 连锁不平衡：Linkage Disequilibrium](#linkage-disequilibrium)
	- [2.1. 什么是连锁不平衡](#introduction)
	- [2.2. 定量连锁不平衡](#quantify-ld)
		- [2.2.1. r2 和D'的比较](#compare-r2-with-d-plus)
	- [2.3. 可视化连锁不平衡](#visualize-ld)
		- [2.3.1. 热图可视化haplotype block](#heatmap-for-haplotype-block)
		- [2.3.2. LD衰减曲线](#ld-decay-curve)






<h1 name="title">群体遗传学知识点</h1>

<a name="fst"><h2>1. Fst及计算方法 [<sup>目录</sup>](#content)</h2></a>

Fst是什么？

> Fst：群体间遗传分化指数，是种群分化和遗传距离的一种衡量方法，分化指数越大，差异越大。适用于亚群体间多样性的比较。
> 
> 用于衡量种群分化程度，取值从0到1，为0则认为两个种群间是随机交配的，基因型完全相似；为1则表示是完全隔离的，完全不相似。它往往从基因的多样性来估计，比如SNP或者microsatellites(串联重复序列一种，长度小于等于10bp)。是一种以哈温平衡为前提的种群遗传学统计方法。

Fst的计算公式如下：

<p align="center"><img src=./picture/PopulationGenetics-Fst-formula.png width=200 /></p>

Hs：亚群体中的平均杂合度
Ht：复合群体中的平均杂合度

在遗传学中，F一词通常代表“近亲繁殖”，它倾向于减少群体中的遗传变异。遗传变异可以用杂合度来衡量，所以F一般表示群体中杂合性的减少。 FST是与它们所属的总群体相比，亚群体中杂合性的减少量。

**如何计算Fst？**

以一个例子进行说明：

> 基因SLC24A5是黑色素表达途径的关键部分，其导致皮肤和毛发色素沉着。与欧洲较轻的皮肤色素密切相关的SNP是rs1426654。 SNP有两个等位基因A和G，其中G与轻度皮肤相关，在犹他州的欧裔美国人中，频率为100％。
> 
> 美洲印第安人与美国印第安人混血儿的SNP在频率上有所不同。
> - 墨西哥的样本有38％A和62％G;
> - 在波多黎各，频率分别为59％A和41％G；
> - 查尔斯顿的非裔美国人样本中有19％A和81％G；
> 
> 这个例子中的FST是什么？

<p align="center"><img src=./picture/PopulationGenetics-Fst-example.png width=800 /></p>

Fst值的计算可以使用VCFtools实现

```
# 对每一个SNP变异位点进行计算

vcftools --vcf test.vcf --weir-fst-pop 1_population.txt --weir-fst-pop 2_population.txt  --out p_1_2—single

# 按照区域来计算

vcftools --vcf test.vcf --weir-fst-pop 1_population.txt --weir-fst-pop 2_population.txt  --out p_1_2_bin --fst-window-size 500000 --fst-window-step 50000
```
参数说明：

> - `--vcf`：SNP calling 过滤后生成的vcf 文件；
> - `--weir-fst-pop`：群体组成文件，包含同一个群体中所有个体，一般每行一个个体。个体名字要和vcf的名字对应；
> - `--out`：生成结果的prefix；
> - `--fst-window-size`与`--fst-window-step`：在按照区域计算Fst时，指定计算窗口和步长

<p align="center"><img src=./picture/PopulationGenetics-Fst-plot-1.png width=800 /></p>

<p align="center">对每个SNP位点计算Fst的散点分布图</p>

<p align="center"><img src=./picture/PopulationGenetics-Fst-plot-2.png width=800 /></p>

<p align="center">对每个区块计算Fst的散点分布图</p>

<a name="linkage-disequilibrium"><h2>2. 连锁不平衡：Linkage Disequilibrium [<sup>目录</sup>](#content)</h2></a>

<a name="introduction"><h3>2.1. 什么是连锁不平衡 [<sup>目录</sup>](#content)</h3></a>

连锁不平衡：指的是在某一群体中，两个基因同时遗传的频率大于随机组合的频率

举个例子：

> 基因A的两个allel 分别用A和a 表示，基因B的allel分别用B和b表示，如果这这两个基因完全独立遗传，也就是说其allel 完全随机组合，那么后代中会出现4种单倍型，AB， Ab, aB， ab, 而且出现的概率都是相同的，都是0.25；如果这两个基因在遗传时不是独立的，意味着后代中单倍型出现的概率不是完全随机的了，我们就可以说两个基因是存在连锁关系的，基因在遗传时出现连锁的现象就叫做连锁不平衡。

那么，如何定量描述连锁不平衡？

<a name="quantify-ld"><h3>2.2. 定量连锁不平衡 [<sup>目录</sup>](#content)</h3></a>

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

<a name="compare-r2-with-d-plus"><h4>2.2.1. r2 和D'的比较 [<sup>目录</sup>](#content)</h4></a>

<p align="center"><img src=./picture/LD-relation-Dplus-r2.jpg width=800 /></p>

r2和D'反映了LD的不同方面。r2包括了重组和突变，而D'只包括重组史

D'能更准确地估测重组差异，但样本较小时，低频率等位基因组合可能无法观测到，导致LD强度被高估，所以D'不适合小样本群体研究；

LD衰减作图中通常采用r2来表示群体的LD水平；Haplotype Block中通常采用D'来定义Block；

<a name="visualize-ld"><h3>2.3. 可视化连锁不平衡 [<sup>目录</sup>](#content)</h3></a>

理论上来说任意两个基因之间都可能存在连锁不平衡，但是实际操作中，认为只有一定区间范围内的基因会存在连锁不平衡，距离大于区间的基因，两者出现连锁不平衡的概率非常小，所以就不去计算。

<a name="heatmap-for-haplotype-block"><h4>2.3.1. 热图可视化haplotype block [<sup>目录</sup>](#content)</h4></a>

对于连锁不平衡的结果，通常采用heatmap热图的形式进行展示，`haploview` 给出的LD heatmap 示例如下：
	
<p align="center"><img src=./picture/LD-visulize-LD-heatmap.jpg width=800 /></p>
	
颜色从白色到红色，代表连锁程度从低到高，方框中的数值为r2，为了美观，这里将 r2 乘以了100

**haplotype block**，即单体型块，即连锁不平衡区域，是指同一条染色体上处于连锁不平衡状态的一段连续的区域。单体型块分析可以用于筛选tag SNP、确定候选基因的范围等

<p align="center"><img src=./picture/LD-haplotype-block-1.jpg width=600 /></p>

<p align="center">Nature Genetics 48, 927–934 (2016) doi:10.1038/ng.3596</p>

<p align="center"><img src=./picture/LD-haplotype-block-2.jpg width=600 /></p>

<p align="center">Plant Biotechnol J. 2017 Mar 29. doi: 10.1111/pbi.12734</p>

<a name="ld-decay-curve"><h4>2.3.2. LD衰减曲线 [<sup>目录</sup>](#content)</h4></a>

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

(1) [【简书】Fst的计算原理与实战](https://www.jianshu.com/p/b73a8d6233be)

(2) [【简书】连锁不平衡：linkage disequilibrium](https://www.jianshu.com/p/477bdb7f57fa)

(3) [【百迈客基因】群体遗传分析—LD连锁不平衡](https://mp.weixin.qq.com/s/1y27LMskwKXx_PfPSNr1hw)
