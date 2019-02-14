<a name="content">目录</a>

[宏基因组shotgun分析套路](#title)
- [1. MGWAS](#mgwas)
	- [1.1. 构建基因集 (gene catalogue)](#gene-catalogue-construction)
	- [1.2. 宏基因组组成定量](#quantification-of-metagenome-content)
		- [1.2.1. Gene&Funtianl profiles](#funtional-profile)
		- [1.2.2. Taxonomic profiles](#taxonomic-profiles)
			- [1.2.2.1. RefMG.v1](#taxonomic-profiles-refmg)
			- [1.2.2.2. mOTU](#taxonomic-profiles-motu)





<h1 name="title">宏基因组shotgun分析套路</h1>

<p align="center"><img src=./picture/Strategies-metagenome-outline.png width=800 /></p>

主要包括三大部分：

- **基础分析**，包括数据预处理、denovo拼接和基因预测
- **构建gene catalog**（图中蓝色部分），包括基因去冗余，非冗余基因的功能注释
- **定量分析**（图中橙色部分），包括Gene&Funtianl profile（基因定量是为了后续功能分析作准备）与物种定量（taxonomic profile）
	- **Functional profiles**： 将 HQ reads mapping 到特定环境的参考基因集(human gut, human skin, mouse gut or the oceans)，也可以使用自己构建的基因集
	- **Taxonomic profiles**： 
		- 将 HQ reads mapping 到 **RefMG.v1** 数据库，该数据库中包含细菌基因组的单拷贝基因，因此可以得到 NCBI taxonomy 数据库水平的物种定量；
		- 也可以用 **mOTU** 数据库，基于已经鉴定出来的40个 protein-coding phylogenetic marker genes (MGs)——在大多数已知的细菌物种中，它们都是单拷贝的，经常被用作原核生物种水平的鉴定

<p align="center"><img src=./picture/Metagenome-flowchart-2.png width=800 /></p>

<a name="mgwas"><h2>1. MGWAS [<sup>目录</sup>](#content)</h2></a>

宏基因组关联分析（Metagenome-Wide Association Study，MWAS）

> 类似于全基因组关联分析（Genome-Wide Association Analysis, GWAS），宏基因组关联分析(MWAS)是在宏基因组范围内检测复杂疾病相关的微生物变化的一种群体关联分析方法。MWAS研究不仅能鉴定疾病相关的微生物种类差异、丰度差异，还能够发现相关微生物功能的增强或减弱，可以高深度、多角度研究微生物组与复杂疾病的关联，具有重要的科研和临床意义。

<a name="gene-catalogue-construction"><h3>1.1. 构建基因集 (gene catalogue) [<sup>目录</sup>](#content)</h3></a>

**（1）采用MetaHIT项目中相同的策略来构建肠道微生物基因集**

> - **序列拼接**：用SOAPdenovo进行序列拼接
> 
> - **基因预测与去冗余**：用GeneMark进行基因预测，然后对所有预测出来的基因用BLAST进行双序列比对，将那些90%序列长度能比对上其他序列的基因丢弃来去除冗余
> 
> - **合并MetaHIT中的基因集并去冗余**：

目前已经构建好的参考基因集 (pre-compiled reference gene catalogs)

> **IGC (human gut)**: Li,J. et al. (2014) An integrated catalog of reference genes in the human gut microbiome. Nat. Biotechnol., 32, 834–41.
> 
> **CRC-RGC (human gut)**: Zeller,G. et al. (2014) Potential of fecal microbiota for early-stage detection of colorectal cancer. Mol. Syst. Biol., 10, 766.
> 
> **skin-RGC (human skin)**: Oh,J. et al. (2014) Biogeography and individuality shape function in the human skin metagenome. Nature, 514, 59–64.
> 
> **mouse-RGC (human skin)**: Xiao,L. et al. (2015) A catalog of the mouse gut metagenome. Nat Biotech, 33, 1103–1108.
> 
> **OM-RGC (ocean)**: Sunagawa,S. et al. (2015) Structure and function of the global ocean microbiome. Science, 348 (6237), 1:10

**（2）确定基因的物种归属 (Taxonomic assignment of genes)**

> 将上一步得到的基因序列与 IMG 数据库中的微生物参考基因组进行比对，来确定基因的物种归属
> 
> 在MetaHIT发表的关于人的微生物肠型 (enterotype) 的文章中，对从系统发生角度进行物种鉴定的相似性参数进行了比较全面的探索，采样该文章推荐的参数：对于**种 (genus) 水平的鉴定**采用85%的 identity 阈值，以及80%的 alignment coverage；对于**门 (phylum) 水平的鉴定**采用65%的identity阈值

**（3）基因的功能注释**

> 将之前得到的基因序列进行翻译，得到其假定的氨基酸序列，然后用 BLASTP 与 eggNOG 和 KEGG 数据库的 proteins/domains 进行比对 (e-value ≤ 1e-5)
> 
> 依据 BLASTP 比对的 score 值最高的 hit(s)，将每个蛋白质归属到对应的 KEGG orthologue group (KO) 或 eggNOG orthologue group (OG) 中
> 
> 对于那些在eggNOG数据库中找不到注释信息的基因（孤儿基因），使用MCL进行新基因家族的鉴定 (inflation factor of 1.1, bit-score cutoff of 60)

<a name="quantification-of-metagenome-content"><h3>1.2. 宏基因组组成定量 [<sup>目录</sup>](#content)</h3></a>

<a name="funtional-profile"><h4>1.2.1. Gene&Funtianl profiles [<sup>目录</sup>](#content)</h4></a>

<a name="taxonomic-profiles"><h4>1.2.2. Taxonomic profiles [<sup>目录</sup>](#content)</h4></a>

<a name="taxonomic-profiles-refmg"><h5>1.2.2.1. RefMG.v1 [<sup>目录</sup>](#content)</h5></a>

<a name="taxonomic-profiles-motu"><h5>1.2.2.2. mOTU [<sup>目录</sup>](#content)</h5></a>

[mOTU.v1.padded](http://www.bork.embl.de/software/mOTU/download.html) 包含10个 MGs，代表了3,445个已知的原核参考基因组和一些未知的物种（基于可公开获取的263个metagenomes，数据来着 MetaHIT 和 HMP 计划）

如何构建mTOU MGs？



> 使用 [fetchMG](http://www.bork.embl.de/software/mOTU/download.html) 从已知细菌参考基因组和metagenomes中提取出这40个 MGs。它基于HMM进行序列搜索，从 [eggNOG](http://eggnogdb.embl.de/) 数据库得到的40个 MGs的多序列比对结果训练Hidden Markov Models([HMMER](http://hmmer.org/))

mOTU本来作为一个独立 (Stand-alone) 的分析工具被开发出来，后来被整合到MOCAT流程中

- **独立运行mOTU进行Taxonomic profiles**

	先下载`mOTUs.pl`，下载地址：`http://www.bork.embl.de/software/mOTU/download.html`，下载后要解压，解压后得到`mOTUs.pl`
	
	接着运行`mOTUs.pl`进行Taxonomic profiles
	
	```
	# 单样本single-end
	$ perl mOTUs.pl sample.fq.gz

	# 单样本paired-end
	$ perl motus.pl sample.1.fq.gz sample.2.fq.gz

	# 多样本对文件的组织形式有要求
	# 1）一个样本一个文件夹
	# 2）创建一个文本文件保存这些文件夹名（一个文件夹名一行）
	$ perl motus.pl --sample-file <SAMPLE_FILE>
	```

	注意：在首次运行该脚本时，它会在工作目录下创建一个`motus_data`文件夹，用来保存需要的dependencies

	mOTU profiles 结果会被保存在`motu.profiles/mOTU.v1.padded/`

	taxonomic profiles 结果会被保存在`taxonomic.profiles/Ref10.v1.padded/`



- **作为MOCAT的组件进行Taxonomic profiles**

	MOCAT提供了2中方法进行Taxonomic profiles：

	> - 使用`runMOCAT.sh`脚本
	> - 逐步执行MOCAT的命令

	**（1）方案一：使用`runMOCAT.sh`脚本**

	在安装MOCAT后，新建一个项目专用文件夹，例如`MOCAT_analysis`，然后将`MOCAT.cfg`拷贝到该文件夹下

	在`MOCAT_analysis`文件夹下，为每一个样本创建一个自文件夹，即`MOCAT_analysis/sample1`，然后将属于该样本的`the .fq(.gz)`文件保存到该文件夹下

	接着创建一个`sample file`，保存需要分析的样本的样本名，一个样本一行

	这样就可以执行`sh runMOCAT.sh`脚本进行分析了，该脚本以**交互形式**执行

	```
	$ runMOCAT.sh
	
	##############################################################################
	# WELCOME TO THE MOCAT EXECUTER v1.3 #
	##############################################################################
	
	This shell script is used to execute a number of MOCAT commands in a row.
	Typically this is used to process raw reads up to final taxonomic or mOTU
	profiles. Of course you can process each step individually using MOCAT.pl
	but we have created this software for your ease to execute these commands
	with ease without prior knowledge of how to run MOCAT. Enjoy!
	
	Usage: runMOCAT.sh [-sf SAMPLE_FILE -cfg CONFIG_FILE]
	
	SAMPLE_FILE not specified with option -sf SAMPLE_FILE
	Looking for valid sample files in the current folder:
	Getting files...
	Getting folders...
	Processing files................
	
	SELECT A SAMPLE FILE:
	- sample
	
	ENTER SAMPLE FILE:
	
	--- Type 'sample' and press enter ---
	```

	然后选择要执行的功能模块：

	```
	AVAILABLE SCRIPTS:
	1: assemble_revise_predict_genes_no_hg19_screen
	process raw reads, assemble, revise assembly and predict genes
	2: assemble_revise_predict_genes_with_hg19_screen
	process raw reads, remove human contaminants, assemble, revise assembly and predict genes
	3: taxonomic_and_motu_profiles_no_hg19_screen
	First process raw reads and then generate taxonomic and mOTU profiles
	4: taxonomic_and_motu_profiles_with_hg19_screen
	First process raw reads, remove humans reads and generate taxonomic and mOTU profiles
	
	STEP TO EXECUTE (enter number):
	--- Type '3' and press enter ---
	```

	**（2）逐步执行MOCAT的命令**

	```
	# 1. Initial sample processing
	$ MOCAT.pl -sf samples -rtf
	
	# 2. Generate mOTU profiles
	$ MOCAT.pl -sf samples -s mOTU.v1.padded -identity 97
	$ MOCAT.pl -sf samples -f mOTU.v1.padded -identity 97
	$ MOCAT.pl -sf samples -p mOTU.v1.padded -identity 97 -mode mOTU -o RESULTS
	
	# 3. Generate taxonomic profiles
	$ MOCAT.pl -sf samples -s RefMG.v1.padded -r mOTU.v1.padded -e -identity 97
	$ MOCAT.pl -sf samples -f RefMG.v1.padded -r mOTU.v1.padded -e -identity 97
	$ MOCAT.pl -sf samples -p RefMG.v1.padded -r mOTU.v1.padded -e -identity 97 -mode RefMG -previous_db_calc_tax_stats_file -o RESULTS
	```
	
---

参考资料：

(1) [MOCAT官方文档](http://vm-lux.embl.de/~kultima/MOCAT/about.html)

(2) [mOTU官网](http://www.bork.embl.de/software/mOTU/)

(3) Quince C, Walker A W, Simpson J T, et al. Shotgun metagenomics, from sampling to analysis[J]. Nature Biotechnology, 2017, 35(9):833.

(4) [【华大科技BGITech】Nature子刊：这些疾病都与肠道微生物相关（上）](https://mp.weixin.qq.com/s?__biz=MjM5NzUyNzU2MA==&mid=2656475118&idx=1&sn=40fb821a3246bfdf7fb389ddb8bf39a8&chksm=bd7a9e898a0d179f265d4b36ce4d3601ec368b855a4f12544c3ead7a436b10c51775b388d9aa&scene=21#wechat_redirect)

(5) Qin, J., et al., A metagenome-wide association study of gutmicrobiota in type 2 diabetes. Nature. 2012.

(6) Arumugam, M. et al. Enterotypes of the human gut microbiome. Nature. 2011. 473, 174-180

(7) Dongen, v. Graph Clustering by Flow Simulation. PhD thesis (2000)
