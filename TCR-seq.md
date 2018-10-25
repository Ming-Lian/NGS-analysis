<a name="content">目录</a>

[TCR-seq文献调研](#title)
- [背景介绍](#introduction)
- [健康个体的免疫组库](#properties-of-a-healthy-repertoire)
- [对低丰度的T细胞克隆具有极高的灵敏度](#ultra-sensitive-detection-of-rare-T-cell-lones)



<h1 name="title">TCR-seq文献调研</h1>

<a name="introduction"><h2>背景介绍 [<sup>目录</sup>](#content)</h2></a>

TCR与BCR的结构：

<table>
<tr>
	<td><p align="center"><img src=./picture/immuSeq-paper-survey-introduction-1.png hegiht=600 /></p></td>
	<td><p align="center"><img src=./picture/immuSeq-paper-survey-introduction-2.png hegiht=600 /></p></td>
</tr>
</table>

TCR与APC的互作：

<p align="center"><img src=./picture/immuSeq-paper-survey-introduction-3.jpg width=600 /></p>

<p align="center">(a)(b) TCR与APC的互作; (c)T细胞中的VDJ基因重排</p>

CDR3区域：

>  互补决定区域（complementarity determining region 3 (CDR3) domain），一般长度为45nt，有VJ junction形式（TCR-α）和VDJ junction形式（TCR-β）

材料选择：用αβ还是γδ？

> TCR有 αβ（外周血中90%）和 γδ（外周血中10%）二聚体形式，一般对TCR的研究都是针对占主体的αβ二聚体形式，且α的CDR3区域是VJ junction，β的CDR3区域是VDJ junction，所以β链与α链相比具有更多的组合形式和连接的多样性，因此TCR-seq一般都是针对β链的CDR3区域

材料选择：用DNA还是RNA？

> DNA
> - 优点：丰富，容易提取且能保持长时间的稳定，且对于每一个TCR subunit，一个细胞只有两个位置有，或者说只有固定的两份拷贝，因此DNA模板分子的数量能反映T细胞的数量
> - 缺点：必须进行PCR扩增来达到足够的测序量，而为了得到进尽可能全的TCR库的组成，使用了多套PCR引物多重PCR方法，使得很容易在PCR过程中引入PCR bias
> 
> RNA
> 
> 使用5' RACE方法进行cDNA的扩增，因此只需要使用一套PCR引物即可
> - 优点：只使用一套PCR引物，极大地降低了PCR bias
> - 缺点：TCR表达水平的变异很大，不能准确地反映T细胞的数量

测序错误的影响及处理方法：

> TCR-seq对测序错误十分敏感，因为只要有一个碱基不同，一条TCR β链就能区别于其他的克隆，一个碱基的测序错误可能在后续的分析中会被错误地鉴定出一个低丰度的新克隆，因此
> 
> **（1）**在进入后续分析之前需要执行严格的质控，但是**（2）**对于深度的TCR-seq则没有这个必要，因为错误的TCR序列总是表现出低丰度的特征，因此通过一个丰度的阈值筛选就可以比较轻松且准确地将这些错误的TCR克隆过滤掉；还有另外一种解决方法**（3）**假设每一种低丰度的克隆都是由测序错误产生的，将它们分别与高丰度的克隆依据序列相似性进行聚类，将高丰度的克隆的序列作为它的正确的序列

详尽彻底的测序不是T细胞库分析的目标，对于旨在阐明样本间差异的比较研究，极端深的TCR-seq也是没有必要的。对应TCR-seq来说，人们主要关心的是：**当取样不完整，测序深度差异比较大时，怎么鉴定一个给定的样本它的TCR组成与其他的样本不同？**

可以计算一个**置信度（confidence)**，来表示一个克隆在给定样本中差异于另一个样本是偶然产生的概率，概率越低则说明越不可能是由于偶然因素导致的，也就是说是真实的差异的可能性比较大

为了能够将两个或多个样本进行比较，需要先将它们的input data进行标准化处理：

> - 以多个样本中的数据量**最少的那个样本为基准**，对其他样本的reads进行无偏好的随机抽样，将它们的input data砍到同一水平 —— 这是在**生态学**研究中常用的方法
> 
> 这是目前免疫组库测序领域常用的标准化方法，但是该标准化方法是否合理？是否还有其他可选的方法？



比较样本间差异或相似度的几个指标：

- **Simpson diversity index**：样本间的多样性的比较
- **Morisita-Horn similarity index**：样本间相似度的比较

<a name="properties-of-a-healthy-repertoire"><h2>健康个体的免疫组库 [<sup>目录</sup>](#content)</h2></a>

- **TCR的多样性/克隆种类**

	TCR β的CDR3序列，长度为45bp，其最大可能的容量为4<sup>45</sup>，考虑一些已知的限制因素，它理论上的多样性也能达到10<sup>11</sup>，然而实际上T细胞在胸腺成熟的过程中要经历两个选择过程：
	
	> - 阳性选择：留下那些能与自身MHC结合的T细胞克隆
	> - 阴性选择：消除那些与自身MHC结合能力过强的T细胞克隆
	
	只有经过这两步筛选，才能产生成熟的且具有功能的T细胞，并进入外周血然后分散到各个组织器官中，因此实际上产生的成熟的T细胞克隆的种类要远远少于理论值
	
	在1999年的Science文章中，有人基于Vβ18 和 Jβ1.4 subset抽样推断，认为TCR β的克隆种类大概为10<sup>6</sup>
	
	2009年，基于深度的TCR-seq和unseen species model，推断TCR β的克隆种类大概为3-4 million，目前对个体进行全面深度的TCR-seq的研究得出结论，一个健康个体大概有**1.3 million**的distinct TCR β chain sequences

- **个体内不同的TCR克隆，其丰度有数量级上的差异**

- **样本间共享的TCR克隆**

	**Public T cells**：个体间共有的相同的T细胞克隆型，由于不同个体偶然产生相同的TCR的可能性极低，一段时间以来，它们一直是一种稀奇的事物。
	
	TCR-seq研究表明，公共T细胞实际上是常见的，这是由于这些跨个体共享的TCR特异性的产生概率增加，以及由于遗传密码的简并，不同的TCR核苷酸序列可以编码相同的TCR氨基酸序列。一个人共有的TCR库所占比例已被证明高达14%，而共有TCR库的真实程度可能还要高得多

	简单来说，就是任意两个健康个体，它们之间共享的TCR克隆种类大约占到各自总克隆种类的14%，两个个体间的共享克隆常见但不多，但是两个以上个体间的共享克隆则少之又少，甚至根本没有

	因此想要通过**比较两组样本间某个克隆的丰度是否存在显著差异基本是不可能的**：若control组和case组各有3个样本，要比较的克隆只在其中的一个或两个样本中有检测到（绝大多数克隆都是这种情况），此时根本无法进行比较！

<a name="ultra-sensitive-detection-of-rare-T-cell-clones"><h2>对低丰度的T细胞克隆具有极高的灵敏度 [<sup>目录</sup>](#content)</h2></a>


在相同的T细胞克隆的混合背景（1 million) 中添加不同量的已知的T细胞克隆作为spike-in

<p align="center"><img src=./picture/immuSeq-paper-survey-01.png width=600 /></p>

其中D克隆在Mix3和G克隆在Mix1中的量最少，都只有10个，但都在后续的分析中成功检测到，说明：**免疫组库测序对低丰度的T细胞克隆具有极高的灵敏**

实际检测到频率与期望的频率基本都十分相近

<p align="center"><img src=./picture/immuSeq-paper-survey-02.jpg width=600 /></p>



---

参考资料：

(1) Harlan Robins, Cindy Desmarais, Jessica Matthis, et al. Ultra-sensitive detection of rare T cell clones[J]. Journal of Immunological Methods, 2012, 375(1-2):14-19.

(2) Woodsworth DJ, Castellarin M, Holt RA. Sequence analysis of T-cell repertoires in health and disease. Genome Med. 2013;5(10):98. Published 2013 Oct 30. doi:10.1186/gm502
