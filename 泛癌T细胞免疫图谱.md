<a name="content">目录</a>

[泛癌T细胞免疫图谱](#title)
- [肺癌](#lung-cancer)
	- [1. 简介](#lung-cancer-introduction)
	- [2. 结果](#lung-cancer-result)
		- [2.1. 亚群分类](#lung-cancer-result-identify-t-cell-clusters)
		- [2.2. 不同的细胞亚群存在明显的组织来源特异性](#lung-cancer-result-distinct-tissue-distributions)
		- [2.3. 克隆扩张](#lung-cancer-result-clonal-expansion)
	- [3. 可借鉴的分析方法](#lung-cancer-referable-methods)
		- [3.1. 单细胞测序数据处理](#lung-cancer-single-cell-sequencing-data-processing)
		- [3.2. TCR分析](#lung-cancer-tcr-analysis)
		- [3.3. 组织分布偏好性或其他bias类型的分析](#bias-analysis)
- [结直肠癌](#colorectal-carcinoma)
	- [3. 可借鉴的分析方法](#colorectal-carcinoma-referable-methods)
		- [3.1. 单细胞测序数据处理（R包）](#colorectal-carcinoma-single-cell-sequencing-data-processing)
		- [3.2. 单细胞表达谱的无监督聚类](#unsupervised-clustering-analysis)




<h1 name="title">泛癌T细胞免疫图谱</h1>

<a name="lung-cancer"><h2>肺癌 [<sup>目录</sup>](#content)</h2></a>

<p align="center"><img src=/picture/T-cell-sequencing-in-cancers-lung-cancer.png /></p>

<a name="lung-cancer-introduction"><h3>1. 简介 [<sup>目录</sup>](#content)</h3></a>

肺癌是一种非常常见的恶性肿瘤，其发病率和致死率在男性中一直位居恶性肿瘤之首。一般来讲，肺癌主要分为小细胞肺癌（约占15%）和非小细胞肺癌（Non-small-cell lung cancer，NSCLC）两大类，而其中非小细胞肺癌又可以分为腺癌（约占40%）、鳞状细胞癌（约占30%）和大细胞癌（约占15%）

Material：

> 研究人员对来自14个药物治疗前非小细胞肺癌患者的**外周血**、**癌旁组织**和**癌组织**的12,346个T细胞进行了**单细胞转录组测序**，全面描绘和解析了肺癌T细胞群体的组成、谱系以及功能状态图谱
> 
> <p align="center"><img src=/picture/T-cell-sequencing-in-cancers-lung-cancer-2.png width=800 /></p>
> 
> 使用 CD3, CD4, CD8 和 CD25 的抗体分离富集对应的T细胞：
> - 细胞毒性T细胞：CD3<sup>+</sup>CD8<sup>+</sup>
> - 辅助T细胞：CD3<sup>+</sup>CD4<sup>+</sup>CD25<sup>low/int</sup>
> - 调节性T细胞：CD3<sup>+</sup>CD4<sup>+</sup>CD25<sup>high</sup>

Result：

> - 跨组织分布的T细胞类群
> 
> - 肿瘤浸润T细胞的组成
> 
> - 肿瘤浸润T细胞亚群间潜在的状态转换关系
> 
> - 提出了新的肺腺癌临床标志物

<a name="lung-cancer-result"><h3>2. 结果 [<sup>目录</sup>](#content)</h3></a>

<a name="lung-cancer-result-identify-t-cell-clusters"><h4>2.1. 亚群分类 [<sup>目录</sup>](#content)</h4></a>

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

<a name="lung-cancer-result-distinct-tissue-distributions"><h4>2.2. 不同的细胞亚群存在明显的组织来源特异性 [<sup>目录</sup>](#content)</h4></a>

使用卡方检验来定量描述这种组织来源偏好性

> **卡方检验**（chi-square test）：
> 
> 描述是统计样本的实际观测值与理论推断值之间的偏离程度，实际观测值与理论推断值之间的偏离程度就决定卡方值的大小，如果卡方值越大，二者偏差程度越大；反之，二者偏差越小；若两个值完全相等时，卡方值就为0，表明理论值完全符合。
> 
> 在本研究中采用样本中某类亚群细胞的实际检测数与其随机检测到的期望个数的比值的卡方值R<sub>O/E</sub>来描述这种偏差

<p align="center"><img src=/picture/T-cell-sequencing-in-cancers-lung-cancer-6.png width=800 /></p>

例如，CD8-C1-LEF1 与 CD4-C1-CCR7 都是 naïve T 细胞，它们都富集在**外周血**中；而CD8-C6-LAYN (CD8+ exhausted T cells) 和 CD4-C9-CTLA4 (CD4+Tregs)只富集在**肿瘤组织**中； CD8-C3-CX3CR1 和 CD4-C3-GNLY 是效应细胞，富集于外周血和癌旁组织中，其高表达趋化因子受体和细胞毒性因子，而低表达T细胞“耗竭”（exhaustion）标记基因

> 肿瘤微环境中的杀伤性CD8 T细胞由于长期接受抗原刺激，会出现被称为“耗竭”（exhaustion）的失能状态

<a name="lung-cancer-result-clonal-expansion"><h4>2.3. 克隆扩张 [<sup>目录</sup>](#content)</h4></a>

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

克隆发生扩张，且在三个组织中都有分布——侧面反映了T细胞克隆在从发生位置到病灶位置转移过程中，**边迁移边扩增**

进一步研究发现一些涉及细胞粘附与迁移的基因的相对高表达，也从侧面印证了它们存在的迁移特性

<p align="center"><img src=/picture/T-cell-sequencing-in-cancers-lung-cancer-11.png width=800 /></p>


除了耗竭细胞，肺癌的浸润CD8 T细胞群体还包含两群与耗竭细胞可能存在状态转换关系的“耗竭前”细胞

“耗竭前”细胞相对于耗竭细胞的比例与肺腺癌病人的预后相关，这就为肺腺癌提供了新的临床标志物

<p align="center"><img src=/picture/T-cell-sequencing-in-cancers-lung-cancer-3.png width=600 /></p>

<a name="lung-cancer-referable-methods"><h3>3. 可借鉴的methods [<sup>目录</sup>](#content)</h3></a>

<a name="lung-cancer-single-cell-sequencing-data-processing"><h4>3.1. 单细胞测序数据处理 [<sup>目录</sup>](#content)</h4></a>

1. 去除核糖体RNA序列

	从RFam数据库中下载核糖体RNA序列，将原始测序数据与其比对，过滤出未比对上的reads

2. 舍弃低质量的样本数据

	若单细胞的文库太小或检测到的表达的基因的数量太少（阈值设为所有细胞表达基因数的中位数减去3倍的方差），则将整个样本的数据舍弃

	若检测到的CD3的TPM（CD3D, CD3E 和 CD3G 的均值）太低，该样本也舍弃

3. 鉴定T细胞的克隆型

	TPM<sub>CD8</sub> > 30 ：CD8<sup>+</sup>

	TPM<sub>CD8</sub> < 3 ：CD8<sup>-</sup>

	TPM<sub>CD4</sub> > 30 ：CD4<sup>+</sup>

	TPM<sub>CD4</sub> > 30 ：CD4<sup>-</sup>

<a name="lung-cancer-tcr-analysis"><h4>3.2. TCR分析 [<sup>目录</sup>](#content)</h4></a>

使用TraCeR工具，从T细胞转录组数据中组装出完整的TCR序列

<p align="center"><img src=/picture/T-cell-sequencing-in-cancers-lung-cancer-12.png width=500 /></p>

<a name="bias-analysis"><h4>3.3. 组织分布偏好性或其他bias类型的分析 [<sup>目录</sup>](#content)</h4></a>

可以利用独立性检验来考察两个变量是否有关系，若两个变量之间存在关联性，它们就会以比较大的概率一起出现，从而表现出偏好性

适用于独立性检验（又称关联分析）的统计学方法有

> - 卡方检验（Chi-Square test）
> 
> - Fisher精确检验（Fisher's exact test）

- **卡方检验**

	<p align="center"><img src=/picture/T-cell-sequencing-in-cancers-lung-method-Chi-Square-test-1.png height=70 /></p>
	
	<p align="center"><img src=/picture/T-cell-sequencing-in-cancers-lung-method-Chi-Square-test-2.png width=500 /></p>

- **Fisher精确检验**

	分析男人女人节食是否有显著区别：
	
	<p align="center"><img src=/picture/T-cell-sequencing-in-cancers-lung-method-Fisher-test-1.png width=500 /></p>

	出现上述情况的概率是：

	<p align="center"><img src=/picture/T-cell-sequencing-in-cancers-lung-method-Fisher-test-2.png height=100 /></p>







<a name="colorectal-carcinoma"><h2>结直肠癌 [<sup>目录</sup>](#content)</h2></a>

<p align="center"><img src=/picture/T-cell-sequencing-in-cancers-colorectal-cancer-1.png width=800 /></p>

<a name="colorectal-carcinoma-referable-methods"><h3>3. 可借鉴的分析方法 [<sup>目录</sup>](#content)</h3></a>

<a name="colorectal-carcinoma-single-cell-sequencing-data-processing"><h4>3.1. 单细胞测序数据处理（R包） [<sup>目录</sup>](#content)</h4></a>

使用R包`HTSeqGenie`完成NGS数据的分析，其本质是整合了多个其他R包完成各个分析步骤：

> - `ShortRead`：quality control reporting
> 
> - `gmapR`
>	- detection of adapter contamination
>	- read alignment versus a reference genome
> - `GenomicRanges`
> 	- counting reads in genomic regions
> 	- read-depth coverage computation

<a name="unsupervised-clustering-analysis"><h4>3.2. 单细胞表达谱的无监督聚类 [<sup>目录</sup>](#content)</h4></a>

该聚类方法名为SC3（Single-Cell Consensus Clustering）

该方法本质上就是K-means聚类，不过在执行K-means聚类的前后进行了一些特殊的操作：

> - **k-means聚类前**：进行了数据预处理，即特征的构造，称为特征工程，该方法中是对输入的原始特征空间进行PCA变换或拉普拉斯矩阵变换，对变换后的新特征矩阵逐渐增加提取的主成分数，来构造一系列新特征；
> - **k-means聚类后**：特征工程构造出来的一系列新特征集合，基于这些新特征集合通过k-means聚类能得到一系列不同的聚类结果，尝试对这些聚类结果总结出consensus clustering

<p align="center"><img src=/picture/T-cell-sequencing-in-cancers-colorectal-cancer-2.png width=800 /></p>

本人比较好奇的地方是：**怎么从一系列不同的聚类结果中总结出consensus clustering？**

> 使用CSPA算法（cluster-based similarity partitioning algorithm）
> 
> （1）对每一个聚类结果按照以下方法构造二值相似度矩阵S：如果两个样本i和j在该聚类结果中被聚到同一个集合中，则它们之间的相似度为1，在二值相似度矩阵中对应的值 S<sub>i,j</sub> = 1，否则S<sub>i,j</sub> = 0；
> 
> （2）对所有的聚类结果的二值相似度矩阵S取平均，得到consensus matrix；
> 
> （3）基于consensus matrix进行层次聚类，得到最终的consensus clustering；







---

参考资料：

(1)  Guo X , Zhang Y , Zheng L , et al. Global characterization of T cells in non-small-cell lung cancer by single-cell sequencing[J]. Nature Medicine, 2018, 24(7).

(2) [Nature Medicine| 张泽民组在单细胞水平绘制肺癌T细胞免疫图谱](http://www.sohu.com/a/237814402_650136)

(3) [AI研习社《数据科学中必须熟知的5种聚类算法》](https://mp.weixin.qq.com/s/VzZ9uoVDbgllFxVuMszD_g)

(4) Stubbington, M. J. T. et al. T cell fate and clonality inference from single-cell transcriptomes[J]. Nat. Methods 13, 329–332 (2016).

(5) [简书·Yan文怡《结合日常生活的例子，了解什么是卡方检验》](https://www.jianshu.com/p/807b2c2bfd9b)

(6) [CSDN·joey周琦《Fisher's exact test( 费希尔精确检验)》](https://blog.csdn.net/u011467621/article/details/47971909)

(7) Zhang L, Yu X, Zheng L, et al.Lineage tracking reveals dynamic relationships of T cells in colorectal cancer[J]. Nature, 2018 Dec;564(7735):268-272.

(8) [HTSeqGenie's Documentation](http://www.bioconductor.org/packages/release/bioc/html/HTSeqGenie.html)

(9) Kiselev, V. Y. et al. SC3: consensus clustering of single-cell RNA-seq data[J]. Nat. Methods 14, 483–486 (2017).
