<a name="content">目录</a>

[文献调研：RNA editing](#title)
- [背景知识](#background)
- [*Sci.Rep*: Symmetrical RNA Editing Events in the Mitochondria of Salvia miltiorrhiza](#sci-rep)
	- [REDItools](#reditools)
- [*Nature*: Dynamic landscape and regulation of RNA editing in mammals](#nature)
	- [科普：mmPCR-seq](#mmpcr-seq)
- [*Nat.Meth*: Genome sequence–independent identification of RNA editing sites](#nat-meth)
	- [GIREMI](#giremi)





<h1 name="title">文献调研：RNA editing</h1>

<a name="background"><h2>背景知识 [<sup>目录</sup>](#content)</h2></a>

**RNA editing**: an important post-transcriptional mechanism that alters primary RNAs through the **insertion/deletion** or **modification of specific nucleotides**

脱氨基作用 (deamination)： 如 C->U，或腺苷脱氨酶 (adenosine deaminase (ADAR) family)作用于dsRNA导致 A->I

影响：

- **UTRs**：altered expression, preventing efficient ribosome binding or recognition by small regulatory RNAs
- **coding protein regions**：**amino acid replacements** with variable functional consequences
- In addition： influence the **activity of ncRNAs** such as siRNAs, miRNAs and potentially of piwiRNAs by affecting base-pairing interactions within RNA secondary structures

鉴定原理：

很简单，即将转录本与其对应的基因组序列进行比较，但是要在整个基因组范围内实现准确鉴定仍然充满挑战：**在存在测序错误与mapping不准确的干扰下，怎样从基因组范围内的SNPs中鉴定出真正的RNA editing位点？**

解决方法之一：use DNA-Seq data from single individuals, annotations in dbSNPs and several stringent filters

<a name="sci-rep"><h2>*Sci.Rep*: Symmetrical RNA Editing Events in the Mitochondria of Salvia miltiorrhiza [<sup>目录</sup>](#content)</h2></a>

<p align="center"><img src=./picture/RNA-editing-sci-rep-pipeline.png width=500 /></p>

**参数优化**：maping时bowtie选用的mismatch参数的大小？
- 小的mismatch：低估了RNA-editing位点数
> 1. 许多基因的RNA editing位点十分靠近
> 2. 平均100bp的reads有3.4个位点

- 大的mismatch：高假阳性

**解决策略一**："assembly-base"，用trinity进行拼接`‘-SS_lib_type FR’ `

**解决策略二**："mapping-based"，尝试不同的mismatch值，2~10，最终选择最佳值7

<a name="reditools"><h3>REDItools [<sup>目录</sup>](#content)</h3></a>

包含三个主要脚本，用于处理来自同一样本/个体的DNA-seq和RNA-seq数据

- **REDItoolDnaRNA.py**：检测候选的RNA editing位点，通过比较pre-aligned RNA-Seq 和 DNA-Seq reads（BAM format）获得

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

- **REDItoolKnown.py**：explore the RNA editing potential of RNA-Seq experiments by looking at known events only

- **REDItoolDenovo.py**：不需要重测序数据，只利用RNA-seq数据进行RNA editiong的denovo检测

<a name="nature"><h2>*Nature*: Dynamic landscape and regulation of RNA editing in mammals [<sup>目录</sup>](#content)</h2></a>

研究成果：
- dynamic spatiotemporal patterns
- new regulators of RNA editing
- discovered through an **extensive profiling of A-to-I RNA editing** from the Genotype-Tissue Expression (**GTEx**) project and in hundreds of other primate and mouse samples

结论：
- 与重复coding区域相比，非重复coding区域的RNA editing level在不同组织之间变化比较大
- ADAR1主要编辑repetitive coding sites，ADAR2主要编辑non-repetitive coding sites，而ADAR3能抑制RNA editing

<a name="mmpcr-seq"><h3>科普：mmPCR-seq [<sup>目录</sup>](#content)</h3></a>

a targeted RNA sequencing method that couples **microfluidics-based multiplex PCR** and **deep sequencing**

常规RNA-seq存在的问题：
> - the large dynamic range of RNA expression, which leads to inaccurate quantification of allelic ratios for genes with low-to-moderate expression levels
> - 即RNA丰度差异较大，对于中低丰度的RNA的定量不准

mmPCR-seq优点：
> - uniformly and simultaneously amplify up to 960 loci in 48 samples independently of their gene expression levels and to accurately and cost-effectively measure allelic ratios even for low-quantity or low-quality RNA samples
> 
> - 即成比例扩增RNA片段，而不影响基因表达水平的相对定量，同时能提高对低丰度RNA的灵敏度

<p align="center"><img src=./picture/RNA-editing-nature-mmPCR-seq.png width=600 /></p>

这个测序技术的关键在于进行类似454测序中用到的乳化PCR，即让每个RNA片段处于一个独立的PCR反应环境中

<a name="nat-meth"><h2>*Nat.Meth*: Genome sequence–independent identification of RNA editing sites [<sup>目录</sup>](#content)</h2></a>

当前RNA editing位点鉴别存在的挑战：
- 需要来自同一样本的genome sequence data 来过滤SNPs的影响
- 即使提供了genome sequence data，但是由于测序覆盖度（sequencing coverage）不一致等原因，使得仍然无法完全去除SNPs的干扰

> 其他不需要genome sequence data的方法：use multiple RNA-seq data sets to increase the confidence of finding individual sites
>
> 存在的问题：this precludes analysis of single data sets and may miss unique changes

开发的新工具：[**GIREMI**](https://www.ibp.ucla.edu/research/xiao/GIREMI.html)

优点：不需要genome sequence即可进行RNA editing位点的准确鉴定，即使RNA-seq dataset只有较低的测序深度

<a name="giremi"><h3>GIREMI [<sup>目录</sup>](#content)</h3></a>

鉴别RNA-editing/SNP的原理：

<p align="center"><img src=./picture/RNA-editing-nat-meth-GIREMI-principle.png width=500 /></p>

> A pair of SNPs in the same read (or read pair, in paired-end sequencing) maintains the same haplotype in the RNA as in reference genomic DNA
>
> In contrast, a SNP and an RNA editing site exhibit variable allelic linkage because RNA editing occurs post-transcriptionally to either allele randomly
>
> the allelic linkage for a pair of RNA editing sites may also appear random





参考资料：

(1) Ernesto Picardi, Graziano Pesole; REDItools: high-throughput RNA editing detection made easy.[J]. Bioinformatics, 2013, 29:1813–1814.

(2) Wu B, Chen H, Shao J, et al. Identification of Symmetrical RNA Editing Events in the Mitochondria of Salvia miltiorrhiza by Strand-specific RNA Sequencing.[J]. Scientific Reports, 2017, 7:42250.

(3) Tan M H, Li Q, Shanmugam R, et al. Dynamic landscape and regulation of RNA editing in mammals[J]. Nature, 2017, 550(7675):249-254.

(4) Rui Z, Xin L, Ramaswami G, et al. Quantifying RNA allelic ratios by microfluidic multiplex PCR and sequencing[J]. Nature Methods, 2014, 11(1):51.

(5) Zhang Q, Xiao X. Genome Sequence-Independent Identification of RNA Editing Sites[J]. Nature Methods, 2015, 12(4):347.

