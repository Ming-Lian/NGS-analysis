<a name="content">目录</a>

[文献调研：RNA editing](#title)
- [背景知识](#background)
- [*Sci.Rep*: Symmetrical RNA Editing Events in the Mitochondria of Salvia miltiorrhiza](#sci-rep)
	- [REDItools](#reditools)





<h1 name="title">文献调研：RNA editing</h1>

<a name="background"><h2>背景知识 [<sup>目录</sup>](#content")</h2></a>

**RNA editing**: an important post-transcriptional mechanism that alters primary RNAs through the **insertion/deletion** or **modification of specific nucleotides**

脱氨基作用 (deamination)： 如 C->U，或腺苷脱氨酶 (adenosine deaminase (ADAR) family)作用于dsRNA导致 A->I

影响：

- **UTRs**：altered expression, preventing efficient ribosome binding or recognition by small regulatory RNAs
- **coding protein regions**：**amino acid replacements** with variable functional consequences
- In addition： influence the **activity of ncRNAs** such as siRNAs, miRNAs and potentially of piwiRNAs by affecting base-pairing interactions within RNA secondary structures

鉴定原理：

很简单，即将转录本与其对应的基因组序列进行比较，但是要在整个基因组范围内实现准确鉴定仍然充满挑战：**在存在测序错误与mapping不准确的干扰下，怎样从基因组范围内的SNPs中鉴定出真正的RNA editing位点？**

解决方法之一：use DNA-Seq data from single individuals, annotations in dbSNPs and several stringent filters

<a name="sci-rep"><h2>*Sci.Rep*: Symmetrical RNA Editing Events in the Mitochondria of Salvia miltiorrhiza [<sup>目录</sup>](#content")</h2></a>

<p align="center"><img src=./picture/RNA-editing-sci-rep-pipeline.png width=500 /></p>

**参数优化**：maping时bowtie选用的mismatch参数的大小？
- 小的mismatch：低估了RNA-editing位点数
> 1. 许多基因的RNA editing位点十分靠近
> 2. 平均100bp的reads有3.4个位点

- 大的mismatch：高假阳性

**解决策略一**："assembly-base"，用trinity进行拼接`‘-SS_lib_type FR’ `

**解决策略二**："mapping-based"，尝试不同的mismatch值，2~10，最终选择最佳值7

<a name="reditools"><h3>REDItools [<sup>目录</sup>](#content")</h3></a>

包含三个主要脚本，用于处理来自同一样本/个体的DNA-seq和RNA-seq数据

<li><strong>REDItoolDnaRNA.py</strong>：检测候选的RNA editing位点，通过比较pre-aligned RNA-Seq 和 DNA-Seq reads（BAM format）获得</li>

实现步骤：

> 1\. 对RNA-seq数据，逐一扫描基因组位点并返回一个表格，表格中包含
> 
> - coverage depth
> - the mean quality score
> - the observed base distribution
> - the strand if available
> - the list of observed substitutions as well as the frequency of variation
> 
> 2\. 若提供了DNA-seq数据，获得与步骤1中相似的表格数据，用于之后除去潜在的SNPs
> 
> 3\. 对一些位点按照一定规则进行过滤： read coverage, base quality score, mapping quality, bases supporting the variation, type of substitution and frequency
> 
> 4\. 去除一些位点位于：
>
> - homopolymeric regions of predefined length
> - intronic sequences surrounding known splice sites
> - invariant RNA-Seq positions
> - sites not supported by DNA-Seq
> - positions near read ends

<li><strong>REDItoolKnown.py</strong>：explore the RNA editing potential of RNA-Seq experiments by looking at known events only

<li><strong>REDItoolDenovo.py</strong>：不需要重测序数据，只利用RNA-seq数据进行RNA editiong的denovo检测


参考资料：

(1) Ernesto Picardi, Graziano Pesole; REDItools: high-throughput RNA editing detection made easy.[J]. Bioinformatics, 2013, 29:1813–1814

(2) Wu B, Chen H, Shao J, et al. Identification of Symmetrical RNA Editing Events in the Mitochondria of Salvia miltiorrhiza by Strand-specific RNA Sequencing.[J]. Scientific Reports, 2017, 7:42250.

