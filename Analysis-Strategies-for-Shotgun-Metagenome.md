<a name="content">目录</a>

[宏基因组shotgun研究套路](#title)
- [1. MGWAS](#mgwas)
	- [1.1. 构建基因集 (gene catalogue)](#gene-catalogue-construction)
	- [1.2. 宏基因组组成定量](#quantification-of-metagenome-content)





<h1 name="title">宏基因组shotgun研究套路</h1>

<a name="mgwas"><h2>1. MGWAS [<sup>目录</sup>](#content)</h2></a>

宏基因组关联分析（Metagenome-Wide Association Study，MWAS）

> 类似于全基因组关联分析（Genome-Wide Association Analysis, GWAS），宏基因组关联分析(MWAS)是在宏基因组范围内检测复杂疾病相关的微生物变化的一种群体关联分析方法。MWAS研究不仅能鉴定疾病相关的微生物种类差异、丰度差异，还能够发现相关微生物功能的增强或减弱，可以高深度、多角度研究微生物组与复杂疾病的关联，具有重要的科研和临床意义。

<a name="gene-catalogue-construction"><h3>1.1. 构建基因集 (gene catalogue) [<sup>目录</sup>](#content)</h3></a>

（1）采用MetaHIT项目中相同的策略来构建肠道微生物基因集

> - **序列拼接**：用SOAPdenovo进行序列拼接
> 
> - **基因预测与去冗余**：用GeneMark进行基因预测，然后对所有预测出来的基因用BLAST进行双序列比对，将那些90%序列长度能比对上其他序列的基因丢弃来去除冗余
> 
> - **合并MetaHIT中的基因集并去冗余**：

（2）确定基因的物种归属 (Taxonomic assignment of genes)

> 将上一步得到的基因序列与 IMG 数据库中的微生物参考基因组进行比对，来确定基因的物种归属
> 
> 在MetaHIT发表的关于人的微生物肠型 (enterotype) 的文章中，对从系统发生角度进行物种鉴定的相似性参数进行了比较全面的探索，采样该文章推荐的参数：对于**种 (genus) 水平的鉴定**采用85%的 identity 阈值，以及80%的 alignment coverage；对于**门 (phylum) 水平的鉴定**采用65%的identity阈值

（3）基因的功能注释

> 将之前得到的基因序列进行翻译，得到其假定的氨基酸序列，然后用 BLASTP 与 eggNOG 和 KEGG 数据库的 proteins/domains 进行比对 (e-value ≤ 1e-5)
> 
> 依据 BLASTP 比对的 score 值最高的 hit(s)，将每个蛋白质归属到对应的 KEGG orthologue group (KO) 或 eggNOG orthologue group (OG) 中
> 
> 对于那些在eggNOG数据库中找不到注释信息的基因（孤儿基因），使用MCL进行新基因家族的鉴定 (inflation factor of 1.1, bit-score cutoff of 60)

<a name="quantification-of-metagenome-content"><h3>1.2. 宏基因组组成定量 [<sup>目录</sup>](#content)</h3></a>








---

参考资料：

(1) [【华大科技BGITech】Nature子刊：这些疾病都与肠道微生物相关（上）](https://mp.weixin.qq.com/s?__biz=MjM5NzUyNzU2MA==&mid=2656475118&idx=1&sn=40fb821a3246bfdf7fb389ddb8bf39a8&chksm=bd7a9e898a0d179f265d4b36ce4d3601ec368b855a4f12544c3ead7a436b10c51775b388d9aa&scene=21#wechat_redirect)

(2) Qin, J., et al., A metagenome-wide association study of gutmicrobiota in type 2 diabetes. Nature. 2012.

(3) Arumugam, M. et al. Enterotypes of the human gut microbiome. Nature. 2011. 473, 174-180

(4) Dongen, v. Graph Clustering by Flow Simulation. PhD thesis (2000)
