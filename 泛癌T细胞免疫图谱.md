<a name="content">目录</a>

[泛癌T细胞免疫图谱](#title)
- [肺癌](#lung-cancer)
	- [1. 简介](#lung-cancer-introduction)
	- [2. 结果](#lung-cancer-result)
		- [2.1. 亚群分类](#identify-t-cell-clusters)
		- [2.2. 不同的细胞亚群存在明显的组织来源特异性](#distinct-tissue-distributions)
		- [2.3. 克隆扩张](#clonal-expansion)







<h1 name="title">泛癌T细胞免疫图谱</h1>

<a name="lung-cancer"><h2>肺癌 [<sup>目录</sup>](#content)</h2></a>

<p align="center"><img src=/picture/T-cell-sequencing-in-cancers-lung-cancer.png></p>

<a name="lung-cancer-introduction"><h3>1. 简介 [<sup>目录</sup>](#content)</h3></a>

肺癌是一种非常常见的恶性肿瘤，其发病率和致死率在男性中一直位居恶性肿瘤之首。一般来讲，肺癌主要分为小细胞肺癌（约占15%）和非小细胞肺癌（Non-small-cell lung cancer，NSCLC）两大类，而其中非小细胞肺癌又可以分为腺癌（约占40%）、鳞状细胞癌（约占30%）和大细胞癌（约占15%）

Material：

> 研究人员对来自14个药物治疗前非小细胞肺癌患者的**外周血**、**癌旁组织**和**癌组织**的12,346个T细胞进行了**单细胞转录组测序**，全面描绘和解析了肺癌T细胞群体的组成、谱系以及功能状态图谱
> 
> <p align="center"><img src=/picture/T-cell-sequencing-in-cancers-lung-cancer-2.png width=800 /></p>

Result：

> - 跨组织分布的T细胞类群
> 
> - 肿瘤浸润T细胞的组成
> 
> - 肿瘤浸润T细胞亚群间潜在的状态转换关系
> 
> - 提出了新的肺腺癌临床标志物

<a name="lung-cancer-result"><h3>2. 结果 [<sup>目录</sup>](#content)</h3></a>

<a name="identify-t-cell-clusters"><h4>2.1. 亚群分类 [<sup>目录</sup>](#content)</h4></a>

对所有样本的原始的表达谱使用t-SNE方法进行降维
					
<p align="center"><img src=/picture/T-cell-sequencing-in-cancers-lung-cancer-4.png width=600 /></p>
					
得到的聚类结果可以看到与样本的组织来源和亚型相关
					
为了得到内部的亚群，对上一步得到的降维结果进行无监督聚类（densityClust，一种基于密度的聚类算法，类似于Mean-shift）

> **Mean-shift算法**：
>
> Mean-shift 聚类是一个基于滑窗的算法，尝试找到数据点密集的区域。它是一个基于质心的算法，也就是说他的目标是通过更新中心点候选者定位每个组或类的中心点，将中心点候选者更新为滑窗内点的均值。这些候选滑窗之后会在后处理阶段被过滤，来减少临近的重复点，最后形成了中心点的集合和他们对应的组。查看下面的说明图：
>
> <p align="center"><img src=/picture/T-cell-sequencing-in-cancers-lung-cancer-mean-sift-algorithmn.gif width=300 /></p>
> 
> <p align="center">单滑窗的 Mean-Shift 聚类</p>
	
<p align="center"><img src=/picture/T-cell-sequencing-in-cancers-lung-cancer-5.png width=600 /></p>
	
得到16个主要的亚群（亚群名称的最后部分为与该类聚类中心重合的基因的名称），其中：

> - 7个位 CD8+ T 细胞亚群；
> - 7个为CD4+ T 细胞亚群（Tconvs; C1, C2, C3, C4, C5, C6 and C7 of CD4 clusters）；
> - 2个调节性T细胞亚群（Tregs; C8 and C9 of CD4 clusters）；

<a name="distinct-tissue-distributions"><h4>2.2. 不同的细胞亚群存在明显的组织来源特异性 [<sup>目录</sup>](#content)</h4></a>

使用卡方检验来定量描述这种组织来源偏好性

> **卡方检验**（chi-square test）：
> 
> 描述是统计样本的实际观测值与理论推断值之间的偏离程度，实际观测值与理论推断值之间的偏离程度就决定卡方值的大小，如果卡方值越大，二者偏差程度越大；反之，二者偏差越小；若两个值完全相等时，卡方值就为0，表明理论值完全符合。
> 
> 在本研究中采用样本中某类亚群细胞的实际检测数与其随机检测到的期望个数的比值的卡方值R<sub>O/E</sub>来描述这种偏差

<p align="center"><img src=/picture/T-cell-sequencing-in-cancers-lung-cancer-6.png width=800 /></p>

例如，CD8-C1-LEF1 与 CD4-C1-CCR7 都是 naïve T 细胞，它们都富集在**外周血**中；而CD8-C6-LAYN (CD8+ exhausted T cells) 和 CD4-C9-CTLA4 (CD4+Tregs)只富集在**肿瘤组织**中； CD8-C3-CX3CR1 和 CD4-C3-GNLY 是效应细胞，富集于外周血和癌旁组织中，其高表达趋化因子受体和细胞毒性因子，而低表达T细胞“耗竭”（exhaustion）标记基因

> 肿瘤微环境中的杀伤性CD8 T细胞由于长期接受抗原刺激，会出现被称为“耗竭”（exhaustion）的失能状态

<a name="clonal-expansion"><h4>2.3. 克隆扩张 [<sup>目录</sup>](#content)</h4></a>

对16个亚群的8,038个T细胞进行TCR全长测序，得到5,015条 unique TCRs 和 3,023 条重复的TCRs，这意味着存在克隆扩张的现象
	
<p align="center"><img src=/picture/T-cell-sequencing-in-cancers-lung-cancer-7.png width=500 /></p>

扩张范围在2~75之间（即2<sup>1</sup>~2<sup>6</sup>），根据克隆的分布位置可以分为以下两种克隆扩张类型：

> - 组织内部（ intra-tissue）的克隆扩张，某种TCR克隆只在该克隆所在的组织中找到其他相同克隆；
> - 组织间（inter-tissue）的克隆扩张，某种TCR克隆在该克隆所在的组织之外找到其他相同克隆；

大多为组织内部（intra-tissue）的克隆扩张，少数为组织间（inter-tissue）的扩张

<p align="center"><img src=/picture/T-cell-sequencing-in-cancers-lung-cancer-8.png width=500 /></p>

不同细胞亚群的克隆扩张的比例

<p align="center"><img src=/picture/T-cell-sequencing-in-cancers-lung-cancer-9.png width=500 /></p>

可以看到属于naïve T 细胞的两个亚群 CD8-C1-LEF1 和 CD4-C1-CCR7，它们的克隆扩张比例很低。而属于效应T细胞的两个亚群 CD8-C3-CX3CR1 和 CD4-C3-GNLY，它们不仅克隆扩张比例高，而且有很大比例在外周血、癌旁组织和癌组织中都有分布，如下图：

<p align="center"><img src=/picture/T-cell-sequencing-in-cancers-lung-cancer-10.png width=500 /></p>

这两个细胞亚群这么显著的组织间克隆共享，意味着它们很可能来自共同的祖先细胞，并且在外周血与实质组织间迁移，进一步研究发现一些涉及细胞粘附与迁移的基因的相对高表达，也从侧面印证了它们存在的迁移特性

<p align="center"><img src=/picture/T-cell-sequencing-in-cancers-lung-cancer-11.png width=500 /></p>




除了耗竭细胞，肺癌的浸润CD8 T细胞群体还包含两群与耗竭细胞可能存在状态转换关系的“耗竭前”细胞

“耗竭前”细胞相对于耗竭细胞的比例与肺腺癌病人的预后相关，这就为肺腺癌提供了新的临床标志物

<p align="center"><img src=/picture/T-cell-sequencing-in-cancers-lung-cancer-3.png width=600 /></p>




---

参考资料：

(1)  Guo X , Zhang Y , Zheng L , et al. Global characterization of T cells in non-small-cell lung cancer by single-cell sequencing[J]. Nature Medicine, 2018, 24(7).

(2) [Nature Medicine| 张泽民组在单细胞水平绘制肺癌T细胞免疫图谱](http://www.sohu.com/a/237814402_650136)

(3) [AI研习社《数据科学中必须熟知的5种聚类算法》](https://mp.weixin.qq.com/s/VzZ9uoVDbgllFxVuMszD_g)
