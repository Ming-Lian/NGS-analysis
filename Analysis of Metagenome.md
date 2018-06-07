<a name="content">目录</a>

[Analysis of Metagenome](#title)
- [宏基因组binning](#binning)
	- [binning原理](#principle-of-binning)
	- [binning具体操作](#how-to-binning)
	- [目前binning工具存在的问题](#problems-in-binning)



<h1 name="title">Analysis of Metagenome</h1>

<a name="binning"><h2>宏基因组binning [<sup>目录</sup>](#content)</h2></a>

Binning的含义是分箱、聚类，指从微生物群体序列中将不同个体的序列（reads或contigs等）分离开来的过程。简单来说就是把宏基因组数据中来自同一菌株的序列聚到一起，得到一个菌株的基因组。是的，可以达到菌株水平。

宏基因组binning的两方面的重要应用：
> - 关联分析
>
>     通过binning得到的bins（暂且简称为bins，更确切的说是strain-level clusters 或strain-level taxonomic units）可以进行宏基因组关联分析以及多组学联合分析，将特定功能代谢产物与特定物种、特定基因进行关联研究，推动其因果机制的探究，为疾病监控、环境监测提供了菌株水平的生物靶标。
> - 单菌组装
>
>     对binning得到的bins进行后续组装，可以得到很多不能在实验室里培养的细菌、古菌、病毒的基因组草图，然后根据单菌组装结果进行菌株水平的基因和功能注释、比较基因组分析、进化分析

<a name="principle-of-binning"><h3>binning原理 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/Metagenome-gene-binning-1.gif width=800 /></p>

- **根据核酸组成信息来进行binning**

依据：来自同一菌株的序列，其核酸组成是相似的

例如根**据核酸使用频率**（oligonucleotide frequency variations），通常是四核苷酸频率（tetranucleotide frequency），**GC含量**和**必需的单拷贝基因**等

优势：即便只有一个样品的宏基因组数据也可以进行binning，这在原理上是可操作的

不足：由于很多微生物种内各基因型之间的基因组相似性很高，想利用1个样品的宏基因组数据通过核酸组成信息进行binning，效果往往并不理想或难度很大。利用核酸组成信息进行binning，基本上只适合那些群落中物种基因型有明显核酸组成差异的，例如低GC含量和一致的寡核苷酸使用频率

- **根据丰度信息来进行binning**

依据：来自同一个菌株的基因在不同的样品中 ( 不同时间或不同病理程度 ) 的丰度分布模式是相似的【PMID: 24997787】

优势：这种方法更有普适性，一般效果也比较好，能达到菌株的水平

不足：必须要大样本量，一般至少要50个样本以上，至少要有2个组能呈现丰度变化 ( 即不同的处理，不同的时间，疾病和健康，或者不同的采样地点等 ) ，每个组内的生物学重复也要尽量的多

- **同时依据核酸组成和丰度变化信息**

将核酸组成信息和丰度差异信息创建一个综合的距离矩阵，既能保证binning效果，也能相对节约计算资源，现在比较主流的binning软件多是同时依据核酸组成和丰度变化信息

<a name="how-to-binning"><h3>binning具体操作 [<sup>目录</sup>](#content)</h3></a>

**Q1：从哪些序列下手进行binning呢？**

从原始的clean reads，还是从组装成的contig，还是从预测到的gene，都可以。根据基于聚类的序列类型的不同，暂且分为reads binning， contig binning和 genes binning

比较这三种binning的优劣：

> - **contig binning**
> 
>    由于核酸组成和物种丰度变化模式在越长的序列中越显著和稳定，基于contig binning效果可能更好
> - **reads binning**
> 
>    基于reads binning的优势是可以聚类出宏基因组中丰度非常低的物种
>    
>    考虑到在宏基因组组装中reads利用率很低，单样品5Gb测序量情况下，环境样品组装reads利用率一般只有10%左右，肠道样品或极端环境样品组装reads利用率一般能达到30%，这样很多物种，尤其是低丰度的物种可能没有被组装出来，没有体现在gene 或者contig 中，因此基于reads binning 才有可能得到低丰度的物种
>
> - **genes binning**
> 
>    应用非常广泛
>    
>    原因可能是（1）基于genes丰度变化模式进行binning可操作性比较强，宏基因组分析中肯定都会计算gene丰度，一般不会计算contig丰度，gene丰度数据可以信手拈来；（2）基于genes binning有很多可参考的文献，过程也并不复杂，可复制性强；（3）对计算机资源消耗比较低

总体来说应用最广泛的就是基于genes binning 和 contig binning

**Genes binning的一般流程**

在宏基因组做完组装和基因预测之后，把所有样品中预测到的基因混合在一起，去冗余得到unique genes集合，对这个unique genes集合进行binning，主要是根据gene在各个样品中的丰度变化模式，计算gene之间的相关性，利用这种相关性进行聚类

<p align="center"><img src=./picture/Metagenome-gene-binning-2.png width=800 /></p>

该图中的聚类过程类似于**K-means聚类**：随机选择几个seed genes作为诱饵，计算其他基因丰度分布模式与seed genes的相关性，按照固定的相关性值PCC>0.9，将它们归属于不同seed genes所代表的类，然后在聚好的类内重新选择seed genes，进行迭代

<a name="problems-in-binning"><h3>目前binning工具存在的问题 [<sup>目录</sup>](#content)</h3></a>

- 还有很多可提升的空间

比如对核酸组成信息的利用，开发得就不够充分，四碱基使用频率因简单而被广泛使用和接受，但现在已有研究表明k-mer丰度信息也是很好的种系特征，同时越长的k-mer含有越多的信息，还有基因和参考基因组间的同源关系也是有价值的种系信号，但这些都还没有被自动化的binning软件整合

- 对于参数设置是很敏感的，且只有有限的可调整的参数

想要获得高质量的bins经常需要手动调整




参考资料：

(1) [句句干货！一文读懂宏基因组binning](http://baijiahao.baidu.com/s?id=1577425474036936057&wfr=spider&for=pc)

(2) Nielsen H B, Almeida M, Juncker A S, et al. Identification and assembly of genomes and genetic elements in complex metagenomic samples without using reference genomes[J]. Nature Biotechnology, 2014, 32(8):822-828.

(3) Sangwan N, Xia F, Gilbert J A. Recovering complete and draft population genomes from metagenome datasets[J]. Microbiome, 2016, 4(1):8.

