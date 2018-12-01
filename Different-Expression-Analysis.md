<a name="content">目录</a>

[差异表达分析](#title)
- [统计学原理](#statistic-principle)
	- [泊松分布 or 负二项分布？](#poisson-or-negative-binomial-distribution)
- [数据预处理](#data-preprocess)
	- [1. 过滤低表达的基因](#filt-low-exp-genes)
	- [2. 标准化](#normalization)
		- [2.1. CPM](#normalization-cpm)
		- [2.2. TCS](#normalization-tcs)
		- [2.3. Quantile](#normalization-quantile)
		- [2.4. Median of Ration](#normalization-deseq2)
		- [2.5. TMM](#normalization-tmm)





<h1 name="title">差异表达分析</h1>

<a name="statistic-principle"><h2>统计学原理 [<sup>目录</sup>](#content)</h2></a>

<a name="poisson-or-negative-binomial-distribution"><h3>泊松分布 or 负二项分布？ [<sup>目录</sup>](#content)</h3></a>

从统计学的角度出发，进行差异分析肯定会需要假设检验，通常对于分布已知的数据，运用参数检验结果的假阳性率会更低。转录组数据中，raw count值符合什么样的分布呢？

count值本质是reads的数目，是一个非零整数，而且是离散的，其分布肯定也是离散型分布。对于转录组数据，学术界常用的分布包括**泊松分布 (poisson)**和**负二项分布 (negative binomial)**两种。

在数据分析的早期，确实有学者采用泊松分布进行差异分析，但是发展到现在，几乎全部都是基于负二项分布了，究竟是什么因素导致了这种现象呢？为了解释这个问题，我们必须提到一个概念 **overdispersion**。

dispersion指的是离散程度，研究一个数据分布的离散程度，我们常用方差这个指标。**对于泊松分布而言，其均值和方差是相等的，但是我们的数据确不符合这样的规律**。通过计算所有基因的均值和方差，可以绘制如下的图片：

<p align="center"><img src=./picture/DiffExpAna-statistic-principle-1.png width=500 /></p>

横坐标为基因在所有样本中的均值，纵坐标为基因在所有样本中的方差，直线的斜率为1，代表泊松分布的均值和方差的分布。可以看到，真实数据的分布是偏离了泊松分布的，方差明显比均值要大。

如果假定总体分布为泊松分布， 根据我们的定量数据是无法估计出一个合理的参数，能够符合上图中所示分布的，这样的现象就称之为overdispersion。

正是由于真实数据与泊松分布之间的overdispersion， 才会选择负二项分布作为总体的分布。

<a name="data-preprocess"><h2>数据预处理 [<sup>目录</sup>](#content)</h2></a>

<a name="filt-low-exp-genes"><h3>1. 过滤低表达的基因 [<sup>目录</sup>](#content)</h3></a>

> 所有数据集将包括表达的基因和不表达的基因的组合。 虽然检查在一种条件下表达但不在另一种条件下表达的基因是有意义的，但是一些基因在所有样品中都未表达

<p align="center"><img src=./picture/DiffExpAna-filtData.png width=800 /></p>

<a name="normalization"><h3>2. 标准化 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/DiffExpAna-normalization.png width=800 /></p>

<a name="normalization-cpm"><h4>2.1. CPM [<sup>目录</sup>](#content)</h4></a>

CPM(count-per-million)

<p align="center"><img src=./picture/DiffExpAna-normalization-CPM.png height=100 /></p>

<a name="normalization-tcs"><h4>2.2. TCS (Total Count Scaling) [<sup>目录</sup>](#content)</h4></a>

简单来说，就是找出多个样本中library size为中位数的样本，作为参考样本，将所有的样本的library size按比例缩放到参考样本的水平

选择一个library size为中位数的sample，以它为baseline，计算出其它sample对于baseline的normalization factor，即一个缩放因子：

<p align="center"><img src=./picture/DiffExpAna-normalization-TCS-1.png height=100 /></p>

然后基于该缩放因子对特定的sample中的每个基因的read count进行标准化（缩放）：

<p align="center"><img src=./picture/DiffExpAna-normalization-TCS-2.png height=50 /></p>

<a name="normalization-quantile"><h4>2.3. Quantile [<sup>目录</sup>](#content)</h4></a>

简单来说，就是排序后求平均，然后再回序

<p align="center"><img src=./picture/DiffExpAna-normalization-quantile.png width=600 /></p>

在R里面，推荐用preprocessCore 包来做quantile normalization，不需要自己造轮子啦！
但是需要明白什么时候该用quantile normalization，什么时候不应该用，就复杂很多了

<a name="normalization-deseq2"><h4>2.4. Median of Ratio (DESeq2) [<sup>目录</sup>](#content)</h4></a> 

该方法基于的假设是，即使处在不同条件下的不同个样本，大多数基因的表达是不存在差异的，实际存在差异的基因只占很小的部分那么我们只需要将这些稳定的部分找出来，作为标准化的内参，依据内参算出各个样本的标准化因子

（1）对每个基因计算几何平均数，得到一个假设的参考样本(pseudo-reference sample)

| gene | sampleA | sampleB | pseudo-reference sample |
|:---|:---|:---|:---|
| EF2A | 1489 | 906 | sqrt(1489 * 906) = 1161.5 |
| ABCD1 | 22 | 13 | sqrt(22 * 13) = 17.7 |
| … | … | … | … |

（2）对每个样本的每个基因对于参考样本计算Fold Change

|	gene	|	sampleA	|	sampleB	|	pseudo-reference sample	|	ratio of sampleA/ref	|	ratio of sampleB/ref	|
|:---|:---|:---|:---|:---|:---|
|	EF2A	|	1489	|	906	|	1161.5	|	1489/1161.5 = 1.28	|	906/1161.5 = 0.78	|
|	ABCD1	|	22	|	13	|	16.9	|	22/16.9 = 1.30	|	13/16.9 = 0.77	|
|	MEFV	|	793	|	410	|	570.2	|	793/570.2 = 1.39	|	410/570.2 = 0.72	|
|	BAG1	|	76	|	42	|	56.5	|	76/56.5 = 1.35	|	42/56.5 = 0.74	|
|	MOV10	|	521	|	1196	|	883.7	|	521/883.7 = 0.590	|	1196/883.7 = 1.35	|
|	…	|	…	|	…	|	…	|	…	|	…	|


<p align="center"><img src=./picture/DiffExpAna-normalization-DESeq2.png width=600 /></p>

（3）获取每个样本中Fold Change的中位数，我们就得到了非DE基因代表的Fold Change，该基因就是我们选择的该样本的内参基因，它的Fold Change就是该样本的标准化因子

```
normalization_factor_sampleA <- median(c(1.28, 1.3, 1.39, 1.35, 0.59))

normalization_factor_sampleB <- median(c(0.78, 0.77, 0.72, 0.74, 1.35))
```

<a name="normalization-tmm"><h4>2.5. TMM (Trimmed Mean of M value, edgeR) [<sup>目录</sup>](#content)</h4></a> 

该方法的思想与DESeq2的Median of Ratio相同，假设前提都是：大多数基因的表达是不存在差异的

它与DESeq2的不同之处在于对内参的选择上：

> - DESeq2选择一个内参基因，它的Ratio/Fold-Change就是标准化因子
> 
> - edgeR选择一组内参基因集合，它们对标准化因子的计算均有贡献：加权平均

（1）移除所有未表达基因

（2）从众多样本中找出一个数据趋势较为平均的样本作为参考样本

> - 对所有样本求总Read数；
> 
> - 各样本除以各自的总Read数，得到修正Read数；
> 
> - 求出各自样本修正Read数的Q3值（第3个四分位数）；
> 
> - 所有的Q3值求平均，与平均Q3相差最小的样本即是参考样本。

<p align="center"><img src=./picture/DiffExpAna-normalization-TMM-1.jpg width=600 /></p>

（3）找出每个样本中的代表基因集，参考这些代表基因集的fold change，计算出该样本的标准化因子

寻找样本的代表基因集：依据基因的**偏倚程度**和**Reads数**大小选出——偏倚程度小、reads数居中的基因

> - **衡量偏倚度的量：LFC (log fold change)**
> 
> <p align="center"><img src=./picture/DiffExpAna-normalization-TMM-2.png height=80 /></p>
> 
> LFC过大或过小都表示具有偏倚性，LFC越大表示reads数在sample<sub>i</sub>中越高，即偏向sample<sub>i</sub>；LFC越小表示reads数在ref中越高，即偏向ref
> 
> - **衡量reads数的量：read的几何平均数 (read geometric mean, RGM)**
> 
> RGM越大表示基因reads越多，RGM越小表示基因reads越少
> 
> 结合偏倚度、read数画出散点图：
> 
> <p align="center"><img src=./picture/DiffExpAna-normalization-TMM-3.jpg width=600 /></p>
> 
> 偏倚度小、表达量居中的那些基因落在图中的红线附近

由参考代表基因集计算样本的标准化因子：

> 对这些代表基因集计算加权平均数：
> 
> <p align="center"><img src=./picture/DiffExpAna-normalization-TMM-4.png height=100 /></p>
> 
> 该加权平均数就能代表该样本的标准化因子，只是经过了log变换，所以需要恢复为它的正值：
> 
> <p align="center">Scaling-Factor = 2 <sup>weight-average</sup></p>


---

参考资料：

(1) [【生信修炼手册】负二项分布在差异分析中的应用](https://mp.weixin.qq.com/s/m2ydqpKofYo2bK61A9hZWw)

(2) [【生信菜鸟团】quantile normalization到底对数据做了什么？](http://www.bio-info-trainee.com/2043.html)

(3) [Introduction to DGE](https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html)

(4) [生信菜鸟团：StatQuest生物统计学专题 - library normalization进阶之edgeR的标准化方法 ](https://mp.weixin.qq.com/s?__biz=MzUzMTEwODk0Ng==&mid=2247485369&idx=1&sn=791cb8c26b19a1181ceccf586787f078&scene=21#wechat_redirect)


