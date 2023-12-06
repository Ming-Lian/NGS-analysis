
# 目录

- WGS & WES实验介绍与实验设计
  - 实验介绍
  - 实验设计
    - 样本量
    - 测序深度——用泊松分布评估测序覆盖率
- 数据分析：比对
  - 比对原理
  - 比较bowtie2 & bwa mem
  - realign, yes or no
  - 操作
- 数据分析：snp calling的前期处理
  - 去除PCR重复
    - duplicate产生原因
    - 用泊松分布解释duplicate问题
    - PCR bias的影响
    - 操作
  - 碱基质量校正
    - 质量校正原理
    - 操作
- 数据分析：变异检测
  - 场景一：Germline SNP & Indel
    - 变异位点基因型推断的数学原理
      - 单位点基因型推断
      - 单体型推断
    - 过滤可靠变异位点
      - Hard-filtering
      - Soft filtering
    - VCF格式简述
    - 操作
  - 场景二：Somatic SNP & Indel
  - 场景三：Somatic CNV
  - 分开call还是一起call（joint call）？
- snp位点注释
    - 变异产生的影响
    - 可用的注释信息
    - 操作
- 下游进阶分析
  - 癌症驱动基因的鉴定
- 拓展阅读
  - 比对的逻辑：从后缀数组到后缀树，再到BWT
  - BWT算法与动手实现
  - PairHMM
  - 强大的贝叶斯方法
- 附录

# WGS & WES实验介绍与实验设计

## 实验介绍

<p align='center'><img src=./WXS_doc_pic/Experiment-instrution-1.png width=600/></p>

<p align='center'><img src=./WXS_doc_pic/Experiment-instrution-2.png width=600/></p>

## 实验设计

### 样本量

### 测序深度——用泊松分布评估测序覆盖率

假设构建的基因组文库无区域偏好性，测序片段来自于基因组各个区域的概率均等，则我们可以估计特定建库方式下，一定的library size的序列所能覆盖的基因组区域

已知目标基因组的长度为$G$，测序片段长度（read size）为$S$，library size为$N$，则某一条read来自于一个长度为$L$的基因组区段的概率为$\frac{L}{G}$

此时，设随机变量：

$$D=起始于一个长度为L的区段的reads数$$

则D服从二项分布：

$$D \sim Binomial(N,\frac{L}{G}) $$

令$L=S$，则起始于该长度为S的区段的reads，均覆盖该区段的最后一个碱基，即此时D就是该碱基位置的测序深度

我们知道，对于二项分布，当其实验次数$N\to \infty$，概率$P\to 0$时，二项分布近似于泊松分布，在这里，因为$S<<G$，则$P=\frac{S}{G} \to 0$，且$N$非常大，可以用泊松分布近似，即：

$$D \sim Possion(\lambda), 其中\lambda=N\frac{S}{G}$$

因此，我们获得了全基因组各碱基位点的测序深度的概率分布（概率质量分布PMF）估计，如下图（以$\lambda=40$为例）：

<p align='center'><img src=./WXS_doc_pic/snp-calling-estimate-depth-distribution.png width=600/></p>

依据测序深度的概率分布，可以很容易推出测序深度大于指定阈值$d$的基因组区域比例：

$$P(D\ge d) = \sum_{i=d,...,\infty}P(D=i)$$

不过实际的基因组测序深度分布与泊松分布并不完全一致，由于GC偏好性等因素的影响，之前推导过程中依据的全基因组来源等概率的假设并不完全成立，从而使实际的分布相对于理论分布，存在明显的overdispersion（即$Var(D) > E(D)$）

<p align='center'><img src=./WXS_doc_pic/snp-calling-GCbias.png width=400/></p>

<p align='center'>Bentley et al, Nature, 2008</p>

<p align='center'><img src=./WXS_doc_pic/snp-calling-possion-overdispersion.png width=400/></p>

<p align='center'>Shen et al, Nature, 2008</p>

# 数据分析：比对

## 比对原理

精确匹配：

<p align='center'><img src=./WXS_doc_pic/BWT-exactMatch.png width=600/></p>

非精确匹配：

<p align='center'><img src=./WXS_doc_pic/BWT-inexactMatch.png width=600/></p>

## 比较bowtie2 & bwa mem

## realign

一般情况下，比对软件会对gap和mismatch赋予不同的分值，且mismatch的分值一般高于gap，这使得在比对过程中很可能将一个gap错认为mismatch，而**在read的两端这种情况更容易发生**：

> 因为如果在read的中部，将gap错认为mismatch，则会发生其某一侧序列大范围的错位，极大降低比对得分，因此一般不会被认是比对得分最高的最优比对形式，而如果这种情况发生在read的两端，其影响的范围很小，其带来的比对得分的下降也比较小，使其最终得分高于gap形式的比对，而被当作最优比对留下来

例如：

<p align='center'><img src=./WXS_doc_pic/realign-misalign-example.png width=600/></p>

比对结果的左侧，BWA和realign给出的比对结果分别为$\mathrm{C\color{red}{CA}}$和$\mathrm{CCA\color{red}{-}}$，从编辑距离（Editing distance）上来说，前者为2而后者为1，应该是后者的结果更优，但是在比对得分上：

$$\mathrm{C\color{red}{CA}} :1\times 0 + 2 \times (-1)=-2\\
\mathrm{CCA\color{red}{-}}: 3\times 0 + 1\times (=3)=-3$$

前者优于后者

下面realign操作前后的结果比较，可以看出这种差异：

<p align='center'><img src=./WXS_doc_pic/realign-before.png width=600/></p>

<p align='center'><img src=./WXS_doc_pic/realign-after.png width=600/></p>

然而，由于realign过程比较消耗计算资源，且有不少研究人员对有无realign操作下的变异检测结果作了比较，发现结果差异不大，因而提出了一些质疑：**realign步骤是否有必要？**

首先，是否有做realign对最终变异检测结果的影响比较小，这个是有原因的，因为GATK的变异检测不是仅仅基于单碱基位点的比对结果（单点基因型推断），而是将存在比对差异的变异候选区域作为一个整体进行分析，进行局部reads拼接和Haplotype的推断（局部单体型推断，通过HaplotypeCaller实现），这样可以解决由于局部碱基位点比对结果不准确而带来的不良影响，例如，以下形式的局部比对结果

<p align='center'><img src=./WXS_doc_pic/realign-necessary-1.png width=600/></p>

若按照单点基因型推断方法，在该区域会得到三个变异位点，从做到右分别是[+A]、[T->A/C]和[+C]

而通过对该局部进行reads拼接后，我们发现，其本质上只存在一个位点的变异，且为[T->AC]

<p align='center'><img src=./WXS_doc_pic/realign-necessary-2.png width=600/></p>

既然通过局部单体型推断，就可以解决由于局部碱基位点比对结果不准确而带来的不良影响，那为什么还要realign呢？

这是因为碱基质量值校正过程（BQSR）需要依赖准确的比对结果，realign之后得到的准确碱基质量评估有助于最终变异的成功检出，不过其实realign之后得到比对校正的区域/位点较少，因此对碱基质量值校正过程的影响也就比较小——付出这么大代价带来的效果提升很小，投入与产出不匹配，因此在GATK4中就果断放弃了这一步

## 操作

（1）建立参考序列索引

```
$ bwa index -a bwtsw ref.fa
```

参数`-a`用于指定建立索引的算法：
	
> - bwtsw 适用于>10M
> - is 适用于参考序列<2G (默认-a is)

可以不指定`-a`参数，bwa index会根据基因组大小来自动选择合适的索引方法

（2） 序列比对

```
$ bwa mem ref.fa sample_1.fq sample_2.fq -R '@RG\tID:sample\tLB:sample\tSM:sample\tPL:ILLUMINA' \
		2>sample_map.log | samtools sort -@ 20 -O bam -o sample.sorted.bam 1>sample_sort.log 2>&1
```

`-R` 选项为必须选项，用于定义头文件中的SAM/BAM文件中的read group和sample信息

> ```
> The file must have a proper bam header with read groups. Each read group must contain the platform (PL) and sample (SM) tags. 
> For the platform value, we currently support 454, LS454, Illumina, Solid, ABI_Solid, and CG (all case-insensitive)
> 
> The GATK requires several read group fields to be present in input files and will fail with errors if this requirement is not satisfied
> 
> ID：输入reads集的ID号； LB： reads集的文库名； SM：样本名称； PL：测序平台
> ```

查看SAM/BAM文件中的read group和sample信息：
	
```
$ samtools view -H /path/to/my.bam | grep '^@RG'
@RG ID:0    PL:solid    PU:Solid0044_20080829_1_Pilot1_Ceph_12414_B_lib_1_2Kb_MP_Pilot1_Ceph_12414_B_lib_1_2Kb_MP   LB:Lib1 PI:2750 DT:2008-08-28T20:00:00-0400 SM:NA12414  CN:bcm
@RG ID:1    PL:solid    PU:0083_BCM_20080719_1_Pilot1_Ceph_12414_B_lib_1_2Kb_MP_Pilot1_Ceph_12414_B_lib_1_2Kb_MP    LB:Lib1 PI:2750 DT:2008-07-18T20:00:00-0400 SM:NA12414  CN:bcm
@RG ID:2    PL:LS454    PU:R_2008_10_02_06_06_12_FLX01080312_retry  LB:HL#01_NA11881    PI:0    SM:NA11881  CN:454MSC
@RG ID:3    PL:LS454    PU:R_2008_10_02_06_07_08_rig19_retry    LB:HL#01_NA11881    PI:0    SM:NA11881  CN:454MSC
@RG ID:4    PL:LS454    PU:R_2008_10_02_17_50_32_FLX03080339_retry  LB:HL#01_NA11881    PI:0    SM:NA11881  CN:454MSC
...
```

read group信息不仅会加在头信息部分，也会在比对结果的每条记录里添加一个 `RG:Z:*` 标签

```
$ samtools view /path/to/my.bam | grep '^@RG'
EAS139_44:2:61:681:18781    35  1   1   0   51M =   9   59  TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA B<>;==?=?<==?=?=>>?>><=<?=?8<=?>?<:=?>?<==?=>:;<?:= RG:Z:4  MF:i:18 Aq:i:0  NM:i:0  UQ:i:0  H0:i:85 H1:i:31
EAS139_44:7:84:1300:7601    35  1   1   0   51M =   12  62  TAACCCTAAGCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA G<>;==?=?&=>?=?<==?>?<>>?=?<==?>?<==?>?1==@>?;<=><; RG:Z:3  MF:i:18 Aq:i:0  NM:i:1  UQ:i:5  H0:i:0  H1:i:85
EAS139_44:8:59:118:13881    35  1   1   0   51M =   2   52  TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA @<>;<=?=?==>?>?<==?=><=>?-?;=>?:><==?7?;<>?5?<<=>:; RG:Z:1  MF:i:18 Aq:i:0  NM:i:0  UQ:i:0  H0:i:85 H1:i:31
EAS139_46:3:75:1326:2391    35  1   1   0   51M =   12  62  TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA @<>==>?>@???B>A>?>A?A>??A?@>?@A?@;??A>@7>?>>@:>=@;@ RG:Z:0  MF:i:18 Aq:i:0  NM:i:0  UQ:i:0  H0:i:85 H1:i:31
...
```

若原始SAM/BAM文件没有read group和sample信息，可以通过AddOrReplaceReadGroups添加这部分信息

```
$ java -jar picard.jar AddOrReplaceReadGroups \
      I=input.bam \
      O=output.bam \
      RGID=4 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=20
```

# 数据分析：snp calling的前期处理

## 去除PCR重复

### duplicate产生原因

<p align="center"><img src=./WXS_doc_pic/GATK4-pipeline-remove-duplicates-reason-of-duplicates.jpg width=600/></p>

- **PCR duplicates（PCR重复）**

PCR扩增时，同一个DNA片段会产生多个相同的拷贝，第4步测序的时候，这些来源于同！一！个！拷贝的DNA片段会结合到Fellowcell的不同位置上，生成完全相同的测序cluster，然后被测序出来，这些相同的序列就是duplicate

- **Cluster duplicates**

生成测序cluster的时候，某一个cluster中的DNA序列可能搭到旁边的另一个cluster的生成位点上，又再重新长成一个相同的cluster，这也是序列duplicate的另一个来源，这个现象在Illumina HiSeq4000之后的Flowcell中会有这类Cluster duplicates

- **Optical duplicates（光学重复）**

某些cluster在测序的时候，捕获的荧光亮点由于光波的衍射，导致形状出现重影（如同近视散光一样），导致它可能会被当成两个荧光点来处理。这也会被读出为两条完全相同的reads

- **Sister duplicates**

它是文库分子的两条互补链同时都与Flowcell上的引物结合分别形成了各自的cluster被测序，最后产生的这对reads是完全反向互补的。比对到参考基因组时，也分别在正负链的相同位置上，在有些分析中也会被认为是一种duplicates。

### PCR bias的影响

1. DNA在打断的那一步会发生一些损失，主要表现是会引发一些碱基发生颠换变换（嘌呤-变嘧啶或者嘧啶变嘌呤），带来假的变异。PCR过程会扩大这个信号，导致最后的检测结果中混入了假的结果；

2. PCR反应过程中也会带来新的碱基错误。发生在前几轮的PCR扩增发生的错误会在后续的PCR过程中扩大，同样带来假的变异；

3. 对于真实的变异，PCR反应可能会对包含某一个碱基的DNA模版扩增更加剧烈（这个现象称为PCR Bias）。因此， 如果反应体系是对含有reference allele的模板扩增偏向强烈，那么变异碱基的信息会变小，从而会导致假阴。

<p align="center"><img src=./WXS_doc_pic/GATK4-pipeline-remove-duplicates-1.png width=900/></p>

### 用泊松分布解释duplicate问题

求解duplicate rate，相当于是在问这样一个问题：

> 对于已经建好的测序文库，其中有N种序列片段，每条片段长度均为l，每种片段的拷贝数为$k_i(i=1,...,N)$，文库大小（library size）为M，即：
>
> $$M=\sum_{i=1}^{N}k_i$$
>
> 现在，从这个文库M中随机抽取m条序列($m \ll M$)，进行测序
>
> 问：duplicate rate为多少？

先给大家一个结论，duplicate比例为：

$$d \approx 1-\frac{N}{m}\left(1-e^{-m/N}\right)$$



具体的推导过程，请阅读以下部分

解：

下面采用逆向思维来完成这个推导过程

假设，我们可以知道这$m$条序列中总共有$n$种片段，则我们可以很容易地求出目标duplicate rate $d$为：

$$d=1-\frac{n}{m} \tag{1}$$

那么，这个$n$为多少呢？

现在问题变成了：

> 从这个文库$M$中随机抽取$m$条序列($m \ll M$)，理论上我们能抽中多少种片段？

对于原始文库中的任意一种片段$i$，设以下随机变量：

$$
X_i = \left\{ 
  \begin{array}{ll}
    1 & 该片段被至少抽中一次 \\
    0 & 该片段未被抽中
  \end{array}
\right.
$$

注：当随机变量按照以上形式进行设定时，该随机变量的分布函数称为**示性函数**

则总共被抽中的片段种类为：

$$n=\sum_{i=1}^{N}X_i \tag{2}$$

我们需要求出$n$的期望，又

$$E(n)=E(\sum_{i=1}^{N}X_i)=\sum_{i=1}^{N}E(X_i) \tag{3}$$

则我们需要求出其中$E(X_i)$的通式

对于原始文库中的任意一种片段$i$，还可以设以下随机变量：

$$\theta_i=该种片段被抽中的次数$$

则$\theta_i$服从二项分布：$\theta_i \sim Binomial(m, \frac{k_i}{M})$

又由于$\frac{k_i}{M} \to 0$，且$m$也比较大，所以$\theta_i$服从泊松分布，即：$\theta_i \sim Possion(\lambda_i)$，其中$\lambda_i=\frac{m\cdot k_i}{M}$

则

$$
\begin{aligned}
&\quad E(X_i) \\
&= P(X_i=1) \\
&= P(\theta_i\ge 1) \\
&= 1-P(\theta_i=0) \\
&= 1-e^{-\lambda_i} \\
&= 1-e^{-m\cdot k_i/M}
\end{aligned}
\tag{4}
$$

这是我们在已知原始文库中该片段拷贝数$k_i$的情况下，能得出的结果，若我们不知道，则可以假设$k_i \sim Possion(\lambda)$，上式(4)就变成了

$$E(X_i)=E(1-e^{-m\cdot k_i/M})=\sum_{k=0}^{\infty} (1-e^{-m\cdot k_i/M})\cdot p_k \tag{5}$$

而且每种片段被抽中的可能性均满足(5)

所以

$$E(n)=\sum_{i=1}^{N}E(X_i)=N\cdot E(X_i)=N\cdot \sum_{k=0}^{\infty} (1-e^{-m\cdot k_i/M})\cdot p_k \tag{6}$$

还可推出

$$M=\sum_{i=1}^{N}k_i=N\cdot E(k_i)=N\lambda \tag{7}$$

所以

$$
\begin{aligned}
&\quad d \\
&=1-\frac{E(n)}{m} \\
&=1-\frac{N}{m}\sum_{k=0}^{\infty} (1-e^{-m\cdot k_i/M})\cdot p_k \\
&=1-\frac{N}{m}\sum_{k=0}^{\infty} (1-e^{-m\cdot k_i/M})\cdot \frac{\lambda^k}{k!}e^{-\lambda} \\
&=1-\frac{N}{m}+\frac{Ne^{-\lambda}}{m}\sum_{k=0}^{\infty} \frac{1}{k!}(\lambda e^{-m/M})^k
\end{aligned}
\tag{8}
$$

上式(7)中的$\sum_{k=0}^{\infty} \frac{1}{k!}(\lambda e^{-m/M})^k$（其中$\lambda e^{-m/M} \to 0$）近似指数函数$e^x$在$(0, f(0))$处的泰勒展开式：

$$e^x=\sum_{n=0}^{\infty} \frac{x^n}{n!} \tag{9}$$

由于$m \ll M$，则$-m/M \to 0^-$，且由于$e^x$在$x=0$处的一阶泰勒公式为：$e^x=1+x+o(x)$，则此时$e^{-m/M} \approx 1-m/M$

则(8)可以化简为：

$$
\begin{aligned}
&\quad d \\
&=1-\frac{N}{m}+\frac{Ne^{-\lambda}}{m}\sum_{k=0}^{\infty} \frac{1}{k!}(\lambda \underline{e^{-m/M}})^k \\
&\approx 1-\frac{N}{m}+\frac{Ne^{-\lambda}}{m}\sum_{k=0}^{\infty} \frac{1}{k!}[\lambda \underline{(1-m/M)}]^k \\
&\approx 1-\frac{N}{m}\left(1-e^{-\lambda}\cdot \underline{e^{\lambda(1-m/M)}}\right) \\
&=1-\frac{N}{m}\left(1-e^{-\lambda m/M}\right) \\
&=1-\frac{N}{m}\left(1-e^{-m/N}\right)
\end{aligned}
\tag{10}
$$

故，最终得到

$$d\approx 1-\frac{N}{m}\left(1-e^{-m/N}\right)$$


### 如何鉴定出这些duplicate？



### 操作

<p align="center"><img src=./WXS_doc_pic/GATK4-pipeline-remove-duplicates-2.png width=800/></p>

<p align="center"><img src=./WXS_doc_pic/GATK4-pipeline-remove-duplicates-3.png width=800/></p>

1、排序（SortSam）

- 对sam文件进行排序并生成bam文件，将sam文件中同一染色体对应的条目按照坐标顺序从小到大进行排序
- GATK4的排序功能是通过`picard SortSam`工具实现的。虽然`samtools sort`工具也可以实现该功能，但是在GATK流程中还是推荐用picard实现，因为SortSam会在输出文件的头信息部分添加一个SO标签用于说明文件已经被成功排序，且**这个标签是必须的**，GATK需要检查这个标签以保证后续分析可以正常进行

```
java -jar picard.jar SortSam \
  I=input.bam \
  O=sorted.bam \
  SORT_ORDER=coordinate
```

如何检查是否成功排序？

```
$ samtools view -H /path/to/my.bam
@HD     VN:1.0  GO:none SO:coordinate
@SQ     SN:1    LN:247249719
@SQ     SN:2    LN:242951149
@SQ     SN:3    LN:199501827
@SQ     SN:4    LN:191273063
@SQ     SN:5    LN:180857866
@SQ     SN:6    LN:170899992
@SQ     SN:7    LN:158821424
@SQ     SN:8    LN:146274826
@SQ     SN:9    LN:140273252
@SQ     SN:10   LN:135374737
@SQ     SN:11   LN:134452384
@SQ     SN:12   LN:132349534
@SQ     SN:13   LN:114142980
@SQ     SN:14   LN:106368585
@SQ     SN:15   LN:100338915
@SQ     SN:16   LN:88827254
@SQ     SN:17   LN:78774742
@SQ     SN:18   LN:76117153
@SQ     SN:19   LN:63811651
@SQ     SN:20   LN:62435964
@SQ     SN:21   LN:46944323
@SQ     SN:22   LN:49691432
@SQ     SN:X    LN:154913754
@SQ     SN:Y    LN:57772954
@SQ     SN:MT   LN:16571
@SQ     SN:NT_113887    LN:3994
...
```

若随后的比对记录中的contig那一列的顺序与头文件的顺序一致，且在头信息中包含`SO:coordinate`这个标签，则说明，该文件是排序过的

2、标记重复（Markduplicates）

标记文库中的重复

```
java -jar picard.jar MarkDuplicates \
  I=input.bam \
  O=marked_duplicates.bam \
  M=marked_dup_metrics.txt
```

## 碱基质量校正

### 质量校正原理

Phred碱基质量值是由测序仪内部自带的base-calling算法评估出来的，而这种base-calling算法由于受专利保护，掌握在测序仪生成商手中，研究人员并不能了解这个算法的细节，它对于人们来说就是一个黑盒子

而测序仪的base-calling算法给出的质量评估并不十分准确，它带有一定程度的系统误差（非随机误差），使得实际测序质量值要么被低估，要么被高估

<p align="center"><img src=./WXS_doc_pic/BQSR_quality_estimation_error-1.png width=700/></p>

BQSR试图利用机器学习的方法来对原始的测序质量值进行校正

例如：

> 对于一个给定的Run，我们发现，无论什么时候我在测序一个AA 的子序列时，改子序列后紧接着的一个任意碱基的测序错误率总是要比它的实际错误率高出1%，那么我就可以将这样的碱基找出来，将它的原始测序错误率减去1%来对它进行校正

会影响测序质量评估准确性的因素有很多，主要包括序列组成、碱基在read中的位置、测序反应的cycle等等，它们以类似于叠加的形式协同产生影响，这些可能的影响因素被称作协变量 (covariable)

注意：BQSR只校正碱基质量值而不改变碱基组成，特别是对于那些质量值偏低的碱基，我们只能说它被解析成当前碱基组成的准确性很低，但是我们又无法说明它实际更可能是哪种碱基，所以干脆不改

那么，BQSR的工作原理是怎样的？

BQSR本质上是一种回归模型

前提假设：影响质量评估的因素只有reads group来源，测序的cycle和当前测序碱基的序列组成背景（这里将它上游的若干个连续位点的碱基组成看作它的背景，一般为2~6，BQSR中默认为6）

则基于这个前提假设，我们可以得出以下结论：

> 相同reads group来源，同处于一个cycle，且序列背景相同的碱基，它们具有相同的测序错误率，这样的碱基组成一个bin

则可以建立这样的拟合模型：

$$X_i=(RG_i,Cyc_i,Context_i) \quad \begin{matrix} f \\ \to \end{matrix} \quad y_i$$

其中，i表示当前碱基，$RG_i$表示碱基所属的Reads Group来源，$Cyc_i$表示该碱基所在的测序cycle，$Context_i$表示该碱基的序列组成背景，$y_i$表示该碱基的实际测序质量(emprical quality)

这三个分量可以直接通过输入的BAM文件的记录获得，那如何获得实际的实际测序质量呢？

可以通过BAM文件中的比对结果推出

用给定的大型基因组测序计划得到的人群变异位点作为输入，将样本中潜在变异位点与人群注释位点overlap的部分过滤掉，则剩下的那些位点，我们假设它们都是“假”的变异位点，是测序错误导致的误检

则实际测序质量为：

$$EQ=-10\log \frac{\# mismatch + 1}{\#bases + 2}$$

注意：emprical quality是以bin为单位计算出来的

这样，有了X和Y，就可以进行拟合模型的训练了，训练好的模型就可以用于碱基质量值的校正

上述只是BQSR的基本逻辑框架，在实际的实现细节上会稍有一些差别

那质量校正后的效果如何呢？

将校正后的碱基质量和实际碱基质量进行比较，若越靠近实际质量，则能说明得到了比较理想的校正效果

<p align="center"><img src=./WXS_doc_pic/BQSR_quality_estimation_error-2.png width=800/></p>

### 操作

1、建立较正模型

质量值校正，这一步需要用到variants的known-sites，所以需要先准备好已知的snp，indel的VCF文件：

```
# 下载known-site的VCF文件，到Ensembl上下载
$ wget -c -P Ref/mouse/mm10/vcf ftp://ftp.ensembl.org/pub/release-93/variation/vcf/mus_musculus/mus_musculus.vcf.gz >download.log &
$ cd Ref/mouse/mm10/vcf && gunzip mus_musculus.vcf.gz && mv mus_musculus.vcf dbsnp_150.mm10.vcf
# 建好vcf文件的索引，需要用到GATK工具集中的IndexFeatureFile，该命令会在指定的vcf文件的相同路径下生成一个以".idx"为后缀的文件
$ gatk IndexFeatureFile -F dbsnp_150.mm10.vcf

# 建立较正模型
$ gatk BaseRecalibrator -R Ref/mouse/mm10/bwa/mm10.fa -I PharmacogenomicsDB/mouse/SAM/ERR118300.enriched.markdup.bam -O \
PharmacogenomicsDB/mouse/SAM/ERR118300.recal.table --known-sites Ref/mouse/mm10/vcf/dbsnp_150.mm10.vcf
```


2、质量值校准

```
# 质量校正
$ gatk ApplyBQSR -R Ref/mouse/mm10/bwa/mm10.fa -I PharmacogenomicsDB/mouse/SAM/ERR118300.enriched.markdup.bam -bqsr \
PharmacogenomicsDB/mouse/SAM/ERR118300.recal.table -O PharmacogenomicsDB/mouse/SAM/ERR118300.recal.bam
```

# 数据分析：变异检测

## 场景一：Germline SNP & Indel

### 变异位点基因型推断的数学原理

#### 单点基因型推断：Locus based Caller

问题描述

> 某一区域的比对结果如下（目标位点由大写字母标出）：
> 
> ```
> REFERENCE: atcatgacggcaGtagcatat
> --------------------------------
> READ1:     atcatgacggcaGtagcatat
> READ2:         tgacggcaGtagcatat
> READ3:     atcatgacggcaAtagca
> READ4:            cggcaGtagcatat
> READ5:     atcatgacggcaGtagc
> ```
> 
> 该区域存在一个候选变异位点，且在上面的比对结果中用大写字母标出：该碱基位置共用5条reads成功比对上，其中有4条reads在该位置的碱基组成与参考基因组一致，为G，而只有一条read的碱基组成为A，与参考基因组不同，即G:A=4:1
> 
> 该个体该位点似乎是一个G/A杂合，因为按照最简单的杂合位点的定义标准，若一个位点中检出的碱基组成中，与参考不一致的碱基组成比例在20%~80%的范围内，就可是算作是一个杂合位点。显然，在本示例中，是满足这个条件的
> 
> 然而，仅依据碱基组成比例来作这个推断，显然是不够可靠的，因为支持其为G/A杂合的read只有一条，而该read可能来自：
> 
> - 实际存在的变异
> - 建库过程中引入的错误（PCR）
> - 测序过程中错误的碱基检测（base calling）
> - 比对错误
> 
> 因此，为了得到更为可靠的变异检测结果，我们应该综合考虑更多的信息，在给出给出对于位点genotype推断时，同时给出该推断的置信度
> 
> 那么，如何实现呢？

简单来说，就是基于各genotype（$G_i$，其中$i$表示该基因型中与参考基因组一致的allele数）的人群先验概率$P(G_i)$以及特定样本对应位点的测序结果$S$，推断其可能为各个$G_i$的后验概率$P(G | S)$，根据贝叶斯推断的方法，将后验概率最高的那种$G_i$作为其最终推断出的genotype，即

$$G=arg \max_{G_i} P(G_i|S)$$

而根据贝叶斯公式：

$$P(G_i|S)=\frac{P(S|G_i)\cdot P(G_i)}{P(S)}$$

由于测序数据确定，$P(S)$固定，可以省略，则$P(G_i|S) \propto P(S|G_i)\cdot P(G_i)$

另外$P(G_i)$是$G_i$基因型的人群频率，已经给出，所以要想指定$P(G_i|S)$，只需要额外算出$P(S|G_i)$即可

> 注：
> 
> 当只有少数样本时，可以从对应的大型人群重测序项目中获取基因型的人群频率，例如dbSNP
> 
> 当有多个样本构成一个人群时，可以从当前研究的人群测序数据中，根据allele频率推算出基因型频率（Hardy-Weinberg equilibrium (HWE)）。
> 
> 例如，对于某个以A/T形式存在的位点，A的频率为1%（$p=0.01, q=1-p=0.99$），则
> | AA($p^2$) | AT($2pq$) | TT($q^2$) |
> |:--- |:---|:---|
> |0.0001|0.0198|0.9801|

我们假设各个read相互独立，则$P(S|G_i)$就是各个read的$P(S_j | G_i)$的乘积，即

$$P(S|G_i)=\prod_j P(S_j | G_i)$$

则，下面只要分别算出各个$P(S_j | G_i)$，最后就能得到$P(G_i|S)$

那么，怎么算这个$P(S_j | G_i)$呢？

我们令

$$
G_i = \left\{
   \begin{array}{ll}
    G_0 & \textrm{基因型的两个allele与参考基因组都不同} \\
    G_1 & \textrm{基因型的两个allele中有一个与参考基因组都相同}\\
    G_2 & \textrm{基因型的两个allele均与参考基因组相同}
  \end{array} 
  \right.
$$

则

$$
\left\{
  \begin{array}{ll}
  P(S_j | G_0) = \epsilon_j^{I_j}\cdot(1-\epsilon_j)^{1-I_j} \\
  P(S_j | G_1) = \frac{1}{2}(1-\epsilon_j) + \frac{1}{2}\epsilon_j=\frac{1}{2} \\
  P(S_j | G_2) = (1-\epsilon_j)^{I_j} \cdot \epsilon_j^{1-I_j}
  \end{array} 
\right.
$$

其中，$I_j$表示该read当前位点碱基组成是否与参考位点一致，若一致$I_j=1$，否则$I_j=0$

> 对上面的$P(S_j | G_i)$的公式，作一个简单的说明：
> 
> 以上的3个式子来源于下面的同一形式
> 
> $$P(S_j | G_i)=P(I_j=1)^{I_j}\cdot P(I_j=0)^{1-I_j}$$
> 
> 而不同genotype下，$P(I_j=1)$或$P(I_j=0)$因为表示测对和测错对应的事件不同，而得到最终不同的公式

为了加深理解，下面来动手做一个计算：

> 某一基因组区域的比对结果如下：
> 
> <p align='center'><img src=./WXS_doc_pic/snp-calling-mathmatical-theory-individual-site-1.png width=600/></p>
> 
> 其中有一潜在的变异位点：A/C
>
> 假设每个碱基的错误率为0.01，试推段断该位点的基因型？
> 
> 先看第一条read：
> 
> <p align='center'><img src=./WXS_doc_pic/snp-calling-mathmatical-theory-individual-site-2.png width=600/></p>
> 
> $$
> \left\{
  \begin{array}{ll}
> P(S_1|G_0)= P(C|A/A)=\epsilon_j^{I_j}\cdot(1-\epsilon_j)^{1-I_j}=0.01^1\cdot (1-0.01)^{1-1}=0.01\\
> P(S_1|G_1)=P(C|A/C)=0.5\\
> P(S_1|G_2)=P(C|C/C)=(1-\epsilon_j)^{I_j} \cdot \epsilon_j^{1-I_j}=(1-0.01)^1\cdot 0.01^{1-1}=0.99
> \end{array} 
> \right.
> $$
> 
> 再来看下一条read：
> 
> <p align='center'><img src=./WXS_doc_pic/snp-calling-mathmatical-theory-individual-site-3.png width=600/></p>
> 
> $$
> \left\{
  \begin{array}{ll}
> P(S_2|G_0)= P(A|A/A)=\epsilon_j^{I_j}\cdot(1-\epsilon_j)^{1-I_j}=0.01^0\cdot(1-0.01)^{1-0} =0.99\\
> P(S_2|G_1)=P(A|A/C)=0.5\\
> P(S_2|G_2)=P(A|C/C)=(1-\epsilon_j)^{I_j} \cdot \epsilon_j^{1-I_j}=(1-0.01)^0\cdot 0.01^{1-0}=0.01
> \end{array} 
> \right.
> $$
> 
> 后面以此类推，最后得到
> 
> $$
> \left\{
  \begin{array}{ll}
> P(S|G_0)= \prod_i P(S_i|G_0) = 0.00000098 \\
> P(S|G_1)= \prod_i P(S_i|G_1) = 0.03125 \\
> P(S|G_2)= \prod_i P(S_i|G_2) = 0.000097
> \end{array} 
> \right.
> $$
> 
> 另外，给定各基因型的人群先验概率：
> 
> $$
> \left\{
  \begin{array}{ll}
> P(G_0)= P(A/A)=0.04 \\
> P(G_1)= P(A/C)=0.32 \\
> P(G_2)= P(C/C)=0.64
> \end{array} 
> \right.
> $$
> 
> 则可以算出
> 
> $$
> \left\{
  \begin{array}{ll}
> P(A/A|S)=P(G_0|S)=P(G_0)\cdot P(S|G_0)<0.001 \\
> P(A/C|S)=P(G_1|S)=P(G_1)\cdot P(S|G_1)=0.999 \\
> P(C/C|S)=P(G_2|S)=P(G_2)\cdot P(S|G_2)<0.001
> \end{array} 
> \right.
> $$
> 
> 最后，根据贝叶斯推断
> 
> $$G=arg \max_{G_i} P(G_i|S)$$
> 
> 得到最可能的基因型应该是A/C

上面是将三种可能基因型的$P(S_j | G_i)$分别表示出来，为了将它们合并在一个公式中得到更为简洁的表达方式，可采用下面的形式：

$$
\begin{aligned}
&\quad P(S_j|G_i) \\
&=P(I_j=1)^{I_j}\cdot P(I_j=0)^{1-I_j} \\
&=\left[ \frac{i}{2}(1-\epsilon_j)+\frac{2-i}{2}\epsilon_j \right]^{I_j}\cdot \left[ \frac{2-i}{2}(1-\epsilon_j)+\frac{i}{2}\epsilon_j\right]^{1-I_j} \\
&=\frac{1}{4}[i(1-\epsilon_j)+(2-i)\epsilon_j]^{I_j}\cdot [(2-i)(1-\epsilon_j)+i\epsilon_j]^{1-I_j}
\end{aligned}
$$

当$i$分别取0，1和2时，就能得到前面的结果，即：

$$
P(S_j | G_i) = \left\{
   \begin{array}{ll}
    \epsilon_j^{I_j}\cdot (1-\epsilon_j)^{1-I_j}, & i=0 \\
    0.5, & i=1 \\
    (1-\epsilon_j)^{I_j}\cdot \epsilon_j^{1-I_j}, & i=2
  \end{array} 
  \right.
$$


#### 单体型推断：Haplotype based Caller

GATK进行SNP calling的核心算法为HaplotypeCaller，这个也是GATK中最核心的算法，理解了这个算法基本上就明白了GATK变异检测的原理

HaplotypeCaller它本质上是对贝叶斯原理的应用，只是相于同类算法它有点不同之处

算法思想概述：

> HaplotypeCaller首先是根据所测的数据，先构建这个群体中的单倍体组合$H_i$（这也是其名字的由来），由于群体中的单倍体是有多个的，所以最好是多个人一起进行HaplotypeCaller这样构建出来的单倍体组合会越接近真实情况
>
> 构建出单倍体的组合之后（每一个单倍体都有一个依据数据得出的后验概率值），再用每个样本的实际数据$S$去反算它们自己属于各个单倍体组合的后验概率$P(H_i/H_j \mid S)$，这个组合一旦计算出来了，对应位点上的碱基型（或者说是基因型，genotype）也就跟着计算出来了：
>
> 计算每一个后候选变异位置上的基因型（Genotype）后验概率，最后留下基因型（Genotype）中后验概率最高的哪一个

下面进行详细地说明：

在HaplotypeCaller中变异检测过程被分为以下四个大的步骤

<p align=center><img src=./WXS_doc_pic/Algorithms-Bioinf-variants-calling-algorithmn-GATK-1.png width=600/></p>

（1）确定候选变异区域（ActiveRegion）

通过read在参考基因组上的比对情况，筛选出潜在的变异区域，这些区域在GATK中被称为ActiveRegion

<p align=center><img src=./WXS_doc_pic/Algorithms-Bioinf-variants-calling-algorithmn-GATK-5.png width=400/></p>

根据每个碱基位点的序列比对情况，计算每个位点的碱基组成多样性，用信息熵（entropy）来定量描述，得到初始active score，表现为沿基因组坐标延伸的波浪线，称为active profile。随后采用一些数据平滑技术，对active profile进行修正，以突出active score较高的连续的局部。最后通过设置阈值来筛选出所谓的候选变异区域

（2）通过对候选变异区域进行重新组装来确定单倍型

对于每个ActiveRegion，GATK会利用**比对到该区域上的所有read**（如果有多个样本那么是所有这些样本的reads而不是单样本进行）构建一个类似于de Bruijn的图对ActiveRegion进行局部重新组装，构建出该区域中可能的单倍型序列。然后，使用Smith-Waterman算法将每个单倍型序列和参考基因组进行重新比对，重新检测出潜在的变异位点

（3）依据所给定的read比对数据计算各个单倍型的似然值

在步骤2的基础上，我们就得到了在ActiveRegion 中所有可能的单倍型序列，接下来需要评估现有数据中对这些单倍型的支持情况

GATK使用PairHMM算法把原本比对于该区域中的每一条read依次和这些单倍型序列进行两两比对，这样我们就可以得出一个read-单倍型序列成对的似然值矩阵，例如以单倍型为列，以read为行，将矩阵记作$(a_{i,j})$

$$
\begin{array}{l|c|c|c|c}
\hline
0 & H_1 & H_2 & .. & H_m \\
\hline
S_1 & a_{11} & a_{12} & .. & a_{1m} \\
S_2 & a_{21} & a_{22} & .. & a_{2m} \\
.. & .. & .. & .. & .. \\
S_n & a_{n1} & a_{n2} & .. & a_{nm} \\
\hline
\end{array}
$$

则矩阵中的某一个元素$a_{ij}$表示在read $S_i$支持单体型为$H_j$的似然，或者说由$H_j$产生$S_i$的概率，即$P(S_i\mid H_j)=a_{ij}$

例如：

候选单体型：

<p align='center'><img src=./WXS_doc_pic/Algorithms-Bioinf-variants-calling-algorithmn-GATK-2.png width=500/></p>

read-单倍型序列成对的似然值矩阵：

<p align='center'><img src=./WXS_doc_pic/Algorithms-Bioinf-variants-calling-algorithmn-GATK-3.png width=300/></p>

这个似然值矩阵很重要，因为在获得这个矩阵之后，GATK会在每一个潜在的变异位点计算等位基因的边缘概率，这个边缘概率实际上是每一个read在该位点上支持其为变异的似然值

那是如何计算的呢？

以第一个变异位点为例进行说明，该位点为“T/-”二等位型，其可能的基因型组合为T/T、T/-和-/-，以计算T/-的后验概率为例

$$P(T/-\mid S) \propto P(T/-)P(S \mid T/-) \tag{1}$$

要使该位点的基因型为T/-，其可能的二倍体单体型组合有四种形式：$R/H_1$，$R/H_3$，$H_1/H_2$和$H_2/H_3$，则

$$
\begin{aligned}
  &\quad P(S \mid T/-)\\
  &=P(S \mid R/H_1) + P(S \mid R/H_3) + P(S \mid H_1/H_2) + P(S \mid H_2/H_3)
\end{aligned} \tag{2}
$$

各种二倍体单体型组合可以统一记作$H_i/H_j$，则

$$
\begin{aligned}
  &\quad P(S \mid T/-)\\
  &=\sum_{H_i/H_j} P(S \mid H_i/H_j)\\
  &=\sum_{H_i/H_j} \prod_k P(S_k \mid H_i/H_j) \\
  &=\sum_{H_i/H_j} \prod_k \frac{P(S_k\mid H_i)+P(S_k\mid H_j)}{2}
\end{aligned} \tag{3}
$$

而(3)中的$P(S_k\mid H_i)$可以从之前得到的read-单倍型序列成对的似然值矩阵直接获得

然后不知道为什么，GATK并没有采取这种计算思路

该区域有两个候选变异位点，下面以第一个位点为例进行说明

对于该候选变异位点，其可能的allele形式为“T/-”，而在候选的4种haplotype中，Reference与Haplotype2支持“T”，Haplotype1与Haplotype3支持“-”

然后来逐一分析每条read对各个allele的支持概率，从它与Reference和Haplotype2的配对似然值中选择最大值，作为该条read支持“T”的似然值，另外，从它与Haplotype1和Haplotype3的配对似然值中选择最大值，作为该条read支持“-”的似然值

从而得到read-allele的成对似然矩阵

<p align='center'><img src=./WXS_doc_pic/Algorithms-Bioinf-variants-calling-algorithmn-GATK-4.png width=300/></p>

（4）计算每一个样本在最佳单倍型组合下的基因型（Genotype）

在完成了步骤3之后，我们就知道了**每一条read在每个候选变体位点上支持每一种等位基因（Allele）的概率**了。那么，最后要做的就是通过这些似然值，计算出候选变异位点上最可能的样本基因型，也就是Genotype——这也是发现真正变异的过程。这就需要应用贝叶斯原理来完成这个计算了—— GATK这也是到这一步才使用了该原理，通过计算就可以得到每一种Genotype的可能性，最后选择后验概率最高的那一个Genotype作为结果输出至VCF中

后面的分析中，对于每一个变异位点假设只有二等位形式——注意：这和一个ActiveRegion中存在多种单体型不矛盾，若一个ActiveRegion在群体中存在n个变异位点，在只考虑二等位形式的前提下，该区域具有的单体型总共有$2^n$种

下面来推导某个样本中的某一个变异位点最可能的SNP形式

该样本在该位点的genotype为G的后验概率为：

$$P(G \mid D) = \frac{P(G)P(D \mid G)}{\sum_i P(G_i)P(D \mid G_i)} \tag{1}$$

由于分母部分对于任何形式genotype都一样，即它是个定值，所以可以忽略，因此上面的公式可以简化成：

$$P(G \mid D) = P(G)P(D \mid G) \tag{2}$$

其中，$P(G)$为genotype为G的先验概率，理论上为样本来源的群体中allele为G的频率，这个一般需要前期给定，若不给定的话，GATK会默认每种G的频率均等

$P(D \mid G)$表示在已知样本genotype为G的前提下，对样本进行测序得到的测序数据为D（仅考虑该ActiveRegion范围内的）的条件概率，我们假设每条reads之间是相互独立的，所以

$$P(D \mid G)=\prod_j P(D_j \mid G)\tag{3}$$

其中，$D_j$表示该样本测序数据D中的第j条read

由于我们正常人都是二倍体，则对于某一条reads，它既可能来自于同源染色体1，记作$H_1$，也可能开自于同源染色体2，记作$H_2$，所以

$$
\begin{aligned}
&\quad P(D_j \mid G) \\
&= P(D_j,H_1 \mid G) + P(D_j,H_2 \mid G) \\
&= P(H_1 \mid G)P(D_j \mid H_1) + P(H_2 \mid G)P(D_j \mid H_2)
\end{aligned} \tag{4}
$$

由于理论上一条read来源于$H_1$还是$H_2$的概率是均等的，都为1/2，即$P(H_1 \mid G)=P(H_2 \mid G)=1/2$，所以

$$P(D_j \mid G)=\frac{P(D_j \mid H_1)}{2} + \frac{P(D_j \mid H_2)}{2} \tag{5}$$

因此(3)可以改写成

$$P(D \mid G)=\prod_j \left( \frac{P(D_j \mid H_1)}{2} + \frac{P(D_j \mid H_2)}{2}\right) \tag{6}$$

现在如果想算出$P(G \mid D)$，就差$P(D_j \mid H_n)$了，那么，如何算$P(D_j \mid H_n)$呢？

上面已经提到，$P(D_j \mid H_n)$表示的是由同源染色体$H_n$产生read $D_j$的条件概率，而每条同源染色体有它各自的单体型，所以这里可以把$H_n$理解为它对应的单体型，则$P(D_j \mid H_n)$可以理解为在特定单体型$H_n$的前提下，产生read $D_j$的条件概率

### 过滤可靠变异位点

通过上面的一系列操作后，我们一般会得到大量的变异位点，其中不可避免的会存在一些假阳性的检出位点，此时需要采用一些方法和策略来从中过滤提取出高可信的位点，而将其中比较可疑的位点丢弃

此时，有两种过滤策略可供选择，且它们基于这样一个共同的假设——真变异（truth callset）与假变异（false callset）在一些关键指标（metric）上存在差异：

- **Hard filter**：基于对一些关键指标的分析观察，找出truth callset与false callset之间的差异，得到这些关键指标各自的分割点（或称之为阈值），用这些阈值去实现变异位点的筛选

    然而该方法存在一些缺点：对各个指标分别进行阈值比较筛选，本质上是在样本空间上切割出一个边缘平整的多维空间，将truth callset与false callset分别框在该空间的内外两侧，然而这样并不能准确反映这两个集合的真实分布特征；不存在放之四海而皆准的通用筛选阈值，具体阈值的选择，依赖于研究人员对数据特征的了解和个人经验，这是一个比较主观的过程，而且为了简化、统一和可重复，科学共同体往往倾向于使用一套统一的筛选阈值

- **Soft filter**：采用一些统计学习方法，从数据集中学习出所谓的truth callset与false callset的一些关键指标的分布特征，从而得到一个判别模型，那么对于一个待评估的变异位点，将它放到这个模型中，让模型去评判，它到底是更像truth callset还是更像false callset，从而给出判断：若它更像truth callset，那它应该就是真的变异，否则，它就是假的变异

    这种方法完全依赖于从数据集中学习到的判别特征，自动化程度高，减少人工的主观选择引入的偏见，相对更加客观，但是要想得到比较可靠的判别模型，需要足够大量的训练数据集，而这个前提要求常常无法满足，尤其是WES，其所覆盖的变异位点远少于全基因组测序，达不到训练出一个可靠的判别模型的要求，此时，只能退而求其次，使用Hard filter

那么，在此之前需要解决以下这两个问题：

- **truth callset与false callset分别是什么，怎么得到的？**

    常用的truth sets有：
    - **dbSNP**：包含目前研究报道的所有变异位点，但是许多变异位点没有通过严格的验证和评估
    - **对应样本的基因芯片检测结果（Sample-matched genotyping chip）**：检测结果可靠，但是成本高，而且只能对有限的一些位点进行检测
    - **HapMap**：国际人类基因组单体型图计划（International HapMap Project）得到变异位点，其有可靠的证据支持，但是只包含常见的基因组变异（common variantion）
    - **1000 Genomes**：国际千人基因组计划（1000 Genomes Project），其中提供的变异位点信息也比较可靠
    - **Omni 2.5M SNP chip array**：Illumina公司开发的基因型检测芯片（Genotyping BeadChip）中所覆盖的变异位点

    GATK中一般使用HapMap3和Omni中的变异位点作为truth callset

- **采用哪些指标？**

  例如，候选变异位点周围的序列组成背景（sequence context），候选位点的read覆盖深度（Depth of Coverage，DP）和链偏倚性（Strand Bias，SB，即覆盖该位点的reads，分别有多少为正链方向，多少为负链方向）等，均可以作为变异位点的一些描述性指标，这些信息在VCF文件的INFO列

  例如，有一个位点的信息如下：

  <p align=center><img src=./WXS_doc_pic/VQSR-1.png width=600/></p>

#### Hard-filtering

通过一些方法得到一系列所谓的good variants和bad variants，比较分析它们在一些关键指标上具有的一定区分度的差异

<p align=center><img src=./WXS_doc_pic/Hard-filter-01.png width=800/></p>

<p align=center><img src=./WXS_doc_pic/Hard-filter-02.png width=800/></p>

<p align=center><img src=./WXS_doc_pic/Hard-filter-03.png width=800/></p>

基于以上的分布情况，GATK官方给出的Hard-filtering阈值（SNP和Indel分别进行）如下：

> SNP：
> 
> ```
> QD<2.0
> FS>60.0
> SOR>3.0
> MQ<40.0
> MQRankSum < -12.5
> ReadPosRankSum < -8.0
> ```
> 
> Indel：
> 
> ```
> QD < 2.0
> ReadPosRankSum < -20.0
> InbreedingCoeff < -0.8
> FS > 200.0
> SOR > 10.0
> ```

#### Soft-filtering

简单来说，就是将具有一定可信度的已知变异位点作为训练集（training sets，一般选择dbSNP），这些变异位点来自前人的研究，其中既包含ture-positive variants（good variants），也包含false-positive variants（bad variants），此时我们想利用变异位点构成的训练集来训练出一个判别模型，以判别good vs. bad variants——这是很自然的想法

但是现在有一个问题：我们不知道这里面到底哪些是good variants（构成truth sets），哪些是bad variants（构成false sets），怎么办？

此时，如果手上有一批高度可信的变异位点，例如HapMap或Omni，我们可以把它们当作truth sets，这个应该没什么问题，而且truth sets一般来说包含在训练集中，是训练集的子集。那false sets呢？简单，将训练集扣除truth sets，剩下的就是false sets了——不能这么干，对于此时的truth sets，我们只能说它其中的位点基本都是good variants，但是并不能保证它囊括了训练集中所有的good variants，因此直接将训练集扣除truth sets，得到的不是完全由bad variants构成的false sets，其中既有bad variants，也有一部分未被truth sets囊括的good variants。

那么，到底应该怎么获取truth sets和false sets呢？

将HapMap或Omni中的位点作为truth sets，用一些方法来定量评估训练集中每个位点与truth sets特征的相似度，将其中最不像的那些位点当作false sets，那么现在的问题变成了如何定量评估这种相似度？

一个最简单的思路就是**对于任意一个位点，计算它在特征空间中与truth sets中心的距离，距离越近，则说明它越像truth sets**。不过这种做法的局限性很大，它要求比较的目标集合以一个近似于圆形（在三维空间下为球形）的形式分布，而且只能聚集成一个中心团块（可以简单概括为：**圆形分布，一个中心**），如果有多个中心团块，则很难比较，而现实情况下，不仅很难做到圆形或球形的分布，也不能保证只有一个分布中心

例如，有一个truth sets的特征空间分布如下图：

<p align=center><img src=./WXS_doc_pic/VQSR-2.jpg width=600/></p>

它就明显不满足“圆形分布，一个中心”的要求

而且可以看出，truth sets中的各个位点在特征空间中聚集成椭圆形的两团，若划分得再细致一点，右侧的一团可以看作是由左上角的稀疏团块和右下角的致密团块重叠组成，则总共有三个团块，分别记作$\theta_1$，$\theta_2$和$\theta_3$，这便是truth sets特征的直观展示

那么对于其中的任意一点，它与truth sets有多像呢？

可以将这些散点在特征空间的分布看作是truth sets的概率云（类似于电子云，它是电子在电子轨道出现的概率描述），散点分布越密集的位置，其越可能是good variant出现的位置，可以用下面的椭圆形的概率等高线来刻画good variant在各个位置出现的概率：越靠近分布的中心，概率越高

<p align=center><img src=./WXS_doc_pic/VQSR-3.png width=600/></p>

正好可以用二元高斯分布（Gaussian distribution，或者称正态分布，Normal distribution）来对其进行拟合，有三团概率云，需要用三个高斯分布去拟合，这便是所谓的高斯混合模型（Gaussian mixture model）

其对应高斯分布为$p(x \mid \theta_k) \sim Norm(\mu_k, \Sigma_k), k=1,2,3$，至于如何得到高斯混合模型的参数估计，不在本文的讨论范围

> 二元高斯分布的例子
> 
> 二元高斯分布的三维图和它的二维等高线图：
> 
> <p align=center><img src=./WXS_doc_pic/Gaussian-distribution-1.png width=250/><img src=./WXS_doc_pic/Gaussian-distribution-2.png width=250/></p>

对于特征空间的任意一点$x$，只要它离其中一个概率云团块足够近（这里概率的大小可以等价描述距离的远近），甚至直接落入团块的核心区域，它就可能是good variant，因此它可能是good variant，或者说它与truth set的相似度，可以定量表示为所有概率云生成它的概率：

$$p(x)=\sum_{k=1}^3 p(x, \theta_k)=\sum_{k=1}^3 p(\theta_k)p(x \mid \theta_k)$$

$p(x)$越大，则说明它与truth sets越相似

> 这是对高斯混合模型的一种创造性的应用
> 
> 因为传统上，高斯混合模型被看作是一种无监督聚类方法，是一种对k-means算法的概率版本的改进，当得到高斯混合模型（$P(x \mid \theta_k) \sim Norm(\mu_k, \Sigma_k), k=1,2,...,n$）后，对于每一个样本点，通过计算后验概率来确定其应该归属于哪个高斯分布所代表的样本簇，即：
> 
> $$\begin{aligned}
> &\quad\theta^{\ast} \\
> &= \arg \max_{\theta_k} P(\theta_k\mid x)\\
> &=\arg \max_{\theta_k} \frac{P(\theta_k)P(x \mid \theta_k)}{\sum_{\theta_i}P(\theta_i)P(x \mid \theta_i)}
> \end{aligned} $$
> 
> 而在这个应用中，并不关心样本点到底应该落在哪个团块中，或者说应该归属于哪个团块，而是考虑所有团块对它的叠加影响，也就是说，算的是上式中的分母部分，这是一个全概率

得到$P(x)$这个相似度的指标后，对所有候选位点的$P(x)$从大到小进行排序，将其中$P(x)$最大的一部分位点作为最终进行模型训练的truth sets，将其中$P(x)$最小的一部分位点作为false sets，至于$P(x)$阈值到底设多大，GATK VQSR式这么做的：选择合适的$P(x)$阈值，使得初始truth sets中有99.9%（这个比例值可以通过参数`--truth-sensitivity-tranche`调整）位点的$P(x)$值在阈值之上；对于false sets，一般默认选择$P(x)$值最小的1000个位点，当然也可以通过设置参数`--minimum-bad-variants`来调整

下图是利用初始truth sets训练得到的模型，其中左图是各位点x的$P(x)$的概率分布热图，右图是由该模型得到的最终用于模型训练的good variants和bad variants的划分

<p align=center><img src=./WXS_doc_pic/VQSR-4.jpg width=300/><img src=./WXS_doc_pic/VQSR-5.jpg width=300/></p>

以上过程可以概况为下图：

<p align=center><img src=./WXS_doc_pic/VQSR-9.png width=600/></p>

得到正式用于判别模型构建的truth sets和false sets后，对truth sets和false sets采用前面提到的高斯混合模型，分别训练出代表good variants的positive model，和代表bad variants的negative model

对于任意位点，positive model和negative model会分别给出其由对应模型生成的概率，称为似然（likelihood），分别记作$p(x\mid \mathrm{positive})$和$p(x\mid \mathrm{negative})$

按照最大似然估计（maximum likelihood estimation，MLE）的思想，哪个似然值越大，则将其判作哪一类

但是仅考虑大小，而不考虑具体的相对大小的差异程度，往往无法给出可靠的判断：如果两个似然值差异很小，此时我们还能拍着胸脯说实际结果就应该是这个吗？

为了简单起见，假设依据一维的x训练得到的positive model和negative model的似然函数为以下分布形式：

<p align=center><img src=./WXS_doc_pic/VQSR-6.png width=300/></p>

横坐标为位点一维特征的取值，纵坐标为其在对应取值下的两个模型给出的似然值，此时若按照最大似然方法，则相当于将分割线选在两条曲线的交点处

<p align=center><img src=./WXS_doc_pic/VQSR-7.png width=300/></p>

将分割线左侧的位点都预测为good variants，此时上图中红色区域所代表的位点就会被误判，因此需要在灵敏度（sensitivity）特异性（specificity）之间作权衡：分割线往左移动，灵敏度提高但特异性下降；分割线往左移动，灵敏度下降而特异性提高

> 简单解释一下灵敏度和特异性
> 
> 灵敏度：实际为positive的个体中被正确判断为positive的比例，相当于图中分割线左侧True positive部分的面积
> 
> 特异性：实际为negative的个体中被正确判断为negative的比例，相当于图中分割线左侧True negative部分的面积
> 
> <p align=center><img src=./WXS_doc_pic/VQSR-8.png width=300/></p>
> 
> 因此，通过分析分割线左移或右移后，Sensitivity部分和Specificity部分面积的变化，即可以推断出灵敏度和特异性的变化情况

此时设：

$$\mathrm{VQSLOD}(x)=\log \frac{p(x|\mathrm{positive})}{p(x|\mathrm{negative})}$$

将它作为判别的指标并设定相应阈值$\mathrm{VQSLOD}_{\alpha}(\mathrm{VQSLOD}_{\alpha}\ge 0)$，取值大于阈值判定为good variant，否则判定为bad variant。若调大阈值，相当于左移分割线；若调小阈值，相当于右移分割线

至于如何权衡灵敏度和特异性，阈值该设多大，就由用户自己决定了

### VCF格式简述

通过上面的一通分析，我们得到了一系列变异位点，包括SNP和Indel，那么它们在我们的结果中到底长什么样呢？

要回答这个问题，我们需要对VCF格式有一个大致的了解，下面作一个简单的介绍

变异结果最终会保存在VCF文件中，VCF文件一般包括两大部分：头信息（Header）和变异信息（Variant records）

<p align=center><img src=./WXS_doc_pic/VCF-format-1.png width=600/></p>

头信息部分以`##`开头，其中包含了大量的文件注释信息，具体每条注释信息都包含哪些内容，这里不再赘述

变异信息部分，每一行记录一条变异位点的信息，每一列记录该变异位点的各个具体指标，其表头位于头信息的最后一行，第一条变异信息的前一行，以`#`开头

例如，有以下形式的变异检测结果：

| 示例 | 序列 | 说明 |
|:---|:---|:---|
| Ref | $\text{a t C g a}$ | C为参考序列中的allele |
| 1 | $\text{a t G g a}$ | 发生了C->G的替换 |
| 2 | $\text{a t  - g a}$ | 发生C碱基的删除 |
| 3 | $\text{a t CA g a}$ | 发生A碱基的插入 |

它们在VCF文件中表示为以下形式：

```
#CHROM POS ID REF ALT QUAL FILTER INFO
20  3 . C G . PASS DP=100
20  2 . TC  T . PASS DP=100
20  2 . TC  TCA . PASS  DP=100
```

### 操作

1. 变异检测

```
gatk HaplotypeCaller  \
  -R Homo_sapiens_assembly38.fasta \
  -I input.bam \
  -O output.g.vcf.gz \
  -ERC GVCF
```

```
gatk GenotypeGVCFs \
  -R Homo_sapiens_assembly38.fasta \
  -V input.g.vcf.gz \
  -O output.vcf.gz
```


2. 过滤可靠变异

（1）Hard-filtering

选择其中的SNP位点

```
gatk SelectVariants \
  -V cohort.vcf.gz \
  -select-type SNP \
  -O snps.vcf.gz
```

选择其中的Indel位点

```
gatk SelectVariants \
  -V cohort.vcf.gz \
  -select-type INDEL \
  -O indels.vcf.gz
```

对SNP位点执行Hard-filtering

```
gatk VariantFiltration \
  -V snps.vcf.gz \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "SOR > 3.0" --filter-name "SOR3" \
  -filter "FS > 60.0" --filter-name "FS60" \
  -filter "MQ < 40.0" --filter-name "MQ40" \
  -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  -O snps_filtered.vcf.gz
```

对Indel位点执行Hard-filtering

```
gatk VariantFiltration \ 
  -V indels.vcf.gz \ 
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "FS > 200.0" --filter-name "FS200" \
  -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \ 
  -O indels_filtered.vcf.gz
```

（2）Soft-filtering

训练模型

```
gatk VariantRecalibrator \
   -R Homo_sapiens_assembly38.fasta \
   -V input.vcf.gz \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.sites.vcf.gz \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.sites.vcf.gz \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 Homo_sapiens_assembly38.dbsnp138.vcf.gz \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -mode SNP \
   -O output.recal \
   --tranches-file output.tranches \
   --rscript-file output.plots.R
```

分别对SNP和Indel执行Soft-filtering

```
# SNP
gatk ApplyVQSR \
  -V indel.recalibrated.vcf.gz \
  --recal-file ${snps_recalibration} \
  --tranches-file ${snps_tranches} \
  --truth-sensitivity-filter-level 99.7 \
  --create-output-variant-index true \
  -mode SNP \
  -O snp.recalibrated.vcf.gz
# Indel
gatk ApplyVQSR \
  -V cohort_excesshet.vcf.gz \
  --recal-file cohort_indels.recal \
  --tranches-file cohort_indels.tranches \
  --truth-sensitivity-filter-level 99.7 \
  --create-output-variant-index true \
  -mode INDEL \
  -O indel.recalibrated.vcf.gz
```

## 场景二：Somatic SNP & Indel

## 场景三：Somatic CNV

### snp位点注释

#### 变异产生的影响

#### 可用的注释信息

#### 操作

# 下游进阶分析

# 拓展阅读

## 比对的逻辑：从后缀数组到后缀树，再到BWT

对于参考序列T=gggtaaagctataactattgatcaggcgtt，其所有后缀子序列为：

```
gggtaaagctataactattgatcaggcgtt
ggtaaagctataactattgatcaggcgtt
gtaaagctataactattgatcaggcgtt
taaagctataactattgatcaggcgtt
aaagctataactattgatcaggcgtt
aagctataactattgatcaggcgtt
agctataactattgatcaggcgtt
gctataactattgatcaggcgtt
ctataactattgatcaggcgtt
tataactattgatcaggcgtt
ataactattgatcaggcgtt
taactattgatcaggcgtt
aactattgatcaggcgtt
actattgatcaggcgtt
ctattgatcaggcgtt
tattgatcaggcgtt
attgatcaggcgtt
ttgatcaggcgtt
tgatcaggcgtt
gatcaggcgtt
atcaggcgtt
tcaggcgtt
caggcgtt
aggcgtt
ggcgtt
gcgtt
cgtt
gtt
tt
t
```

将其按照字母顺序进行排序后，可以得到：

```
aaagctataactattgatcaggcgtt
aactattgatcaggcgtt
aagctataactattgatcaggcgtt
actattgatcaggcgtt
agctataactattgatcaggcgtt
aggcgtt
ataactattgatcaggcgtt
atcaggcgtt
attgatcaggcgtt
caggcgtt
cgtt
ctataactattgatcaggcgtt
ctattgatcaggcgtt
gatcaggcgtt
gcgtt
gctataactattgatcaggcgtt
ggcgtt
gggtaaagctataactattgatcaggcgtt
ggtaaagctataactattgatcaggcgtt
gtaaagctataactattgatcaggcgtt
gtt
t
taaagctataactattgatcaggcgtt
taactattgatcaggcgtt
tataactattgatcaggcgtt
tattgatcaggcgtt
tcaggcgtt
tgatcaggcgtt
tt
ttgatcaggcgtt
```

这便是后缀数组（suffix array）。此时，若要查询其中的子串，如aagctataacta，则可以利用排序后的后缀子序列，通过二分法进行高效检索

<p align='center'><img src=./WXS_doc_pic/Alignment-algorithmn-1.png width=600/></p>

这个其实就是我们小时候就懂的查字典技巧：比如说，我想用拼音在字典里查一下“牛”字在哪里，怎么办？“牛”字的拼音是niu，那我先在字典的拼音索引里查拼音首字母是“n”的那一部分，再在其中找到拼音前缀为“ni”的部分，最后再锁定到“niu”，一步步缩小检索的范围，最终确定我们要找的少数的几个候选目标，或者就是唯一一个

其搜索长度正比与子串长度n，因此其算法时间复杂度为线性复杂度，即为$O(n)$，但是这是一种以空间换时间的策略，需要消耗$\frac{n(n+1)}{2}$个单位存储空间，因此其算法的空间复杂度为平方复杂度，即为$O(n^2)$

但是相对于直接以排序后的所有后缀子串进行存储，其实还有极大的优化空间：局部紧邻的几个后缀子串存在公共前缀，可以让它们共享一个公共前缀子串，以实现数据压缩：

<p align='center'><img src=./WXS_doc_pic/Alignment-algorithmn-2.png width=600/></p>

<p align='center'><img src=./WXS_doc_pic/Alignment-algorithmn-3.png width=600/></p>

最终得到的压缩结果如下：

<p align='center'><img src=./WXS_doc_pic/Alignment-algorithmn-4.png width=600/></p>

这便是后缀树（suffix trie tree），因此可以说，后缀树来自于后缀数组，是后缀数组数据压缩的一种表示方式。所以，在后缀树中进行子串的查找与后缀数组类似，都是二分查找法：

<p align='center'><img src=./WXS_doc_pic/Alignment-algorithmn-5.png width=600/></p>

同理，我们也可以从后往前倒着查，称为回溯（backtrack），此时需要构建的就是前缀数组（preffix array）和其简化压缩的前缀树（preffix trie tree）

> 这里举一个比前面简单一点的例子，以参考序列为T=AGGAGC为例
> 
> 前缀子串为：
> 
> ```
> AGGAGC
> AGGAG
> AGGA
> AGG
> AG
> A
> ```
> 
> 排序后（从右到左按照字母表顺序进行排序）得到的前缀数组为：
> 
> ```
>      A
>   AGGA
> AGGAGC
>     AG
>  AGGAG
>    AGG
> ```
> 那么，按照从右到左的顺序，对应前缀树从根到叶，建立前缀树，如下：
> 
> <p align='center'><img src=./WXS_doc_pic/Alignment-algorithmn-6.png width=600/></p>
> 
> 如果要查找子串GGAG的话，查找路径为：
> 
> <p align='center'><img src=./WXS_doc_pic/Alignment-algorithmn-7.png width=600/></p>

虽然利用树结构可以有效减少算法的空间复杂度，提高空间利用效率，但是其存储空间的开销仍然不小，而且其数据压缩依赖于前缀子串的共享，例如存在一些极端情况

当$n=m$（m为字符组成类型数，当为碱基序列时，m=4）时，让参考序列T=ACGT，其后缀数组为：

```
ACGT
CGT
GT
T
```

此时相邻两行之间不存在公共前缀子串，无法做任何压缩

当$n=m^2$时，让参考序列T=ACGTTGCAATCTAGGA，其后缀数组为：

```
A
AATCTAGGA
ACGTTGCAATCTAGGA
AGGA
ATCTAGGA
CAATCTAGGA
CGTTGCAATCTAGGA
CTAGGA
GA
GCAATCTAGGA
GGA
GTTGCAATCTAGGA
TAGGA
TCTAGGA
TGCAATCTAGGA
TTGCAATCTAGGA
```

此时相邻两行之间存在最长公共前缀子串的长度为1，只能对第一列做有限的压缩

……

凡此种种不再一一列举，可以从上面推出一个结论：对于长度n=$m^i$参考序列，其最小的最长公共前缀子串的长度为$i-1$，此时只能对后缀数组的初始前$i-1$进行有限程度的共享压缩

因此前缀树/后缀树仍然不是一个足够理想的数据结构以实现快速比对/子串检索，有没有一种更高效的数据存储组织方式呢？

其实不论是前缀数组还是后缀数组，其**任意相邻两列其实已经包含了参考序列所有的排列信息**，一旦选择好相邻两列后，只要充分利用其中排列信息，完全可以做到与前缀数组/后缀数组相近的二分查找效率

以前缀数组的回溯查找为例，令参考序列为T=AGGAGC，为了方便识别参考序列的起始端，在序列的起始位置加上“\^”，则T=\^AGGAGC

其前缀数组为：

```
      ^
     ^A
  ^AGGA
^AGGAGC
    ^AG
 ^AGGAG
   ^AGG
```

此时要查询的目标子序列为S=GGAG

通过从右到左的回溯查找，可以快速地找到含有“AG”后缀的前缀子串

<p align='center'><img src=./WXS_doc_pic/Alignment-algorithmn-8.png width=100/></p>

此时，接着回溯找含有“GAG”后缀的前缀子串，若按照前缀数组的常规回溯查找法，需要看含有“AG”后缀的前缀子串的前一列

<p align='center'><img src=./WXS_doc_pic/Alignment-algorithmn-9.png width=100/></p>

除此之外，其实还有另一个思路：若能将当前碱基对应到前缀数组的最后一列上，例如

<p align='center'><img src=./WXS_doc_pic/Alignment-algorithmn-10.png width=300/></p>

那么回溯查找下一位时，需要查看的就是紧挨着最后一列的前一列，即倒数第二列，查找成功后再对应到最后一列并向前查找，如此循环往复向前回溯，直到目标子序列的第一个碱基，则能查到所有包含目标子序列的位置，否则未到目标子序列的第一个碱基就无法再向前延伸，则查找失败——整个过程只需要查看最后两列

那么如何将当前碱基对应到前缀数组的最后一列上呢？

基于前缀数组/后缀数组推广得到BWT矩阵，BWT矩阵第一列与最后一列所具有的一个特殊性质：第一列当前碱基c在当前列（倒数第二列）中的索引与最后一列的一致，即第一列的第i个c对应最后一列的第i个c

<p align='center'><img src=./WXS_doc_pic/Alignment-algorithmn-12.png width=700/></p>

而依据参考序列中的碱基相对排序，BWT矩阵的最后一行其实是第一列的前一列

因此，若是想利用该方法实现类似与于**前缀数组的回溯查找**，则选择由前缀数组推广得到的BWT矩阵，选择其中的第一列和最后一列，组成保留相邻排列信息的前后两列：

<p align='center'><img src=./WXS_doc_pic/Alignment-algorithmn-13.png width=700/></p>

然后按照从后往前的方向进行回溯查找：

<p align='center'><img src=./WXS_doc_pic/Alignment-algorithmn-14.png width=700/></p>

若是想利用该方法实现类似与于**后缀数组的前向查找**，则选择由前缀数组推广得到的BWT矩阵，选择其中的第一列和最后一列，组成保留相邻排列信息的前后两列，不过这里还要额外对这两列进行排序：

<p align='center'><img src=./WXS_doc_pic/Alignment-algorithmn-15.png width=700/></p>

然后按照从前往后的方向进行前向查找：

<p align='center'><img src=./WXS_doc_pic/Alignment-algorithmn-16.png width=700/></p>

## BWT算法与动手实现

本小节节选自github项目：BioLeetCode，https://github.com/Ming-Lian/Bioinfo_LeetCode

### Burrows-Wheeler Transformation

> 给定：一条核酸序列，例如：TCATC
>
> 任务：将该序列经过BW转换（Burrows-Wheeler Transformation），得到BWT输出，或者称为BWT索引（BWT index）
>
> <p align='center'><img src=./WXS_doc_pic/BioLeetCode_issue_Hard_2.png /></p>
> <p align='center'>Burrows-Wheeler Transformation</p>

其实，这一小题里的任务，以及算法的基本操作过程已经很清楚了，关键是采用什么样的数据结构去实现它

最简单的选择就是数组（以下以Perl语言的逻辑进行说明）：

> 例如：给定的序列为 `TCATC`
> 
> 首先，将其打散成一个字符一个元素的数组，变成`@seq=('T','C','A','T','C')`，并在数组最后追加上一个终止标识符`$`，从而变成`@seq=('T','C','A','T','C','\$')`（注意这个终止标识符属于特殊字符，需要进行转义）
> 
> 随后，将这个数组当作一个队列使用，对队尾元素先出队，然后再将其在队首入队，例如：
> 
> ```perl
> $tail = pop @seq;
> unshift @seq, $tail;
> ```
> 
> 重复执行这样的操作直到从队尾取出一个元素后，新的队尾是终止标识符时截至
> 
> 每执行一次，就将新队形的数组合并成字符串，追加保存到一个新数组中，例如`@BWT`：
> 
> ```perl
>  $queue = join "",@seq;
>  push @BWT, $queue;
> ```
> 最后，将`@BWT`中的字符串元素按照字母顺序逐一进行处理：提取最后一个字母

示例代码（Perl版本），见附录：`BWT_Trans.pl`

### BWT reverse transformation

> 给定：BWT output，以[2.1. Burrows-Wheeler Transformation](#for-veterans-2-1)的输出结果`CCTTA$`为例
>
> 任务：根据BWT output，还原出其原始的输入序列

实现BWT reverse transformation的方法有以下两种：

**（1）方法一**：

根据First Column(FC)与Last Column(LC)之间的存在的以下两种关系：

> - FC与LC在同一行i的碱基，在参考序列中是一后一前紧挨着的关系
>
>    利用这种关系，我们就可以某一个碱基的在参考序列的前一个碱基或者后一个碱基的组成
>
> - 按照从上到下的计数方式，FC中第i次遇到的碱基a与LC中第i次遇到的碱基a，是参考序列中的同一个碱基，即它们在参考序列上对应于同一个碱基位置
>
>    利用这种关系，我们可以FC（或LC）第i行碱基在LC（或FC）所对应的碱基的位置

那么，利用上面的两条性质，我们可以按照下面的方法实现BWT逆变换：

> 初始条件为：当前碱基位置index=0，即为FC的第一行，且当前碱基组成为base=\$，当前参考序列组成为T=\$
> 
> 循环执行下面的操作，直至当前碱基组成再一次为\$：
> 
> （1）回溯：获取LC同处于index行的碱基组成c，更新当前参考序列T=cT
> 
> （2）LC到FC的定位：获取LC处于index行的碱基c在FC对应的碱基的位置i，更新当前碱基位置index=i

<p align='center'><img src=./WXS_doc_pic/BioLeetCode_issue_Hard_2-2-1.png height=300/></p>

那么现在的问题是：该如何分别实现上面的**回溯**和**LC到FC的定位**呢？

对于**回溯**问题，很好解决，我只需要用到基于BWT得到的Last Column（下面为了方便讨论，将其称为BWT数组），则根据已知条件里给定的当前碱基位置index，则可以直接回溯得到LC同处于index行的碱基组成c=BWT[index]

对于**LC到FC的定位**问题，由于FC列中的碱基组成是严格按照碱基优先级递增排列的，因此，只需要知道在参考序列T中碱基优先级小于当前碱基c的所有可能的碱基类型的计数总和Pre(c)，再加上当前碱基是同类碱基的计数次数（Observed Character Count，OCC）OCC(c)，如下图：

<p align='center'><img src=./picture/BioLeetCode_how2deal_Hard_2-2-1.png height=300/></p>

则LC中索引位置为index的碱基c，它再FC中对应的同一碱基的索引位置为：

$$\mathrm{LF(c,index)=Pre(c) + OCC(c)}$$

因此，为了能够实现以上操作，我们需要额外增加以下两个索引信息：

> - **Pre散列**：键为各种可能的碱基类型，键值为参考序列中碱基优先级小于对应碱基的碱基总数
> - **OCC数组**：LC对应行碱基，在LC该行及该行之前出现的次数（0-base）

<p align='center'><img src=./picture/BioLeetCode_how2deal_Hard_2-2-2.png height=400/></p>

其中BWT数组和OCC数组是FM-index（Ferragina & Manzini, 2000）的组成部分，除此之外还包含SA数组（在BWT逆变换中还暂时用不到它，所以在该小节未设置该变量，它在下一小节BWT search中会用到），详细的FM-index，见下方说明：

> 后缀数组（suffix array）的FM-index：
>
> <p align='center'><img src=./picture/BioLeetCode_how2deal_Hard_2-2-3.png height=400/></p>
>
> 以上图的BWT matrix为例，后缀数组由下面三个数组组成：
>
> - **SA数组**：LC对应行碱基在原始序列T中的位置（1-base）
> - **OCC数组**：LC对应行碱基，在LC该行及该行之前出现的次数（0-base），例如对于LC的第5行（0-base）的碱基为a，在LC的第5行及第5行之前，a碱基只出现2次，故OCC[5]=1
> - **BWT数组**：即BWT output，上图的给的例子中即为`gc$aaac`

有了上面定义的三个参考序列的索引信息（BWT数组、OCC数组和Pre散列）后，BWT逆变换算法过程，就可以表示为以下形式：

> 给定：BWT数组、OCC数组和Pre散列
> 
> 初始条件为：当前碱基位置index=0，即为FC的第一行，且当前碱基组成为base=\$，当前参考序列组成为T=\$
> 
> 循环执行下面的操作，直至当前碱基组成再一次为\$：
> 
> （1）回溯：获取LC同处于index行的碱基组成c，即c=BWT[index]，更新当前参考序列T=cT
> 
> （2）LC到FC的定位：获取LC处于index行的碱基c在FC对应的碱基的位置i，即i=Pre[c] + OCC[index]，更新当前碱基位置index=i

示例代码（Perl版本），见附录：`BWT_revTrans_V1.pl`

**（2）方法二**：

根据BWT output逐列推出

<p align='center'><img src=./WXS_doc_pic/BioLeetCode_issue_Hard_2-2-2.png height=300/></p>

通过构建出完整是BWT矩阵，来得到原始的输入序列

可以通过设置两个数组来实现：

> - BWT数组：保存BWT转换后的BWT索引，即@BWT=('C', 'C', 'T', 'T', 'A', '$')
> - Array数组：保存已经得到的部分BWT矩阵，每一个元素为部分BWT矩阵的一行

按照以下流程来完成完整的BWT矩阵的构建：

> - 初始状态：
>
>    @BWT=('C', 'C', 'T', 'T', 'A', '$')
>
>    @Array为空
>
> - 循环执行以下操作，直到Array的列数达到@BWT的长度：
> 
>   （1）@Array**之前**追加一列，为@BWT；
> 
>   （2）对@Array的行按照字母表顺序进行排序，得到重排后的BWT矩阵，以此来更新@Array；

示例代码（Perl版本），见附录：`BWT_revTrans_V2.pl`

# 附录

## BWT_Trans.pl

```perl
#/usr/bin/perl
use strict;
use warnings;

sub BWT_Trans(){
	my ($infile) = @_;
	open IN, "<$infile" or die "$!\n";
	# 解析fasta文件
	my $fa_seq;
	while(<IN>){
		chomp;
		next if(/^>/);
		$fa_seq .= $_;
	}
	close IN;
	# 将以字符串形式组织的序列，打散成数组，并在末尾追加终止标识符
	my @seq = split //,$fa_seq;
	push @seq, "\$";
	print STDERR "Input:\n";
	print STDERR "\t".join("",@seq)."\n\n";
	# 循环出队入队，得到BWT matrix
	my @BWT;
	push @BWT, join("", @seq);
	my $n=0;
	print STDERR "Cyclic rotations:\n";
	print STDERR "\t$n:$BWT[$n]\n";
	my $tail = pop @seq;
	while($seq[$#seq] ne "\$"){
		$n++;
		unshift @seq, $tail;
		push @BWT, join("",@seq);
		print STDERR "\t$n:$BWT[$n]\n";
		$tail = pop @seq;
	}
	print STDERR "\n";
	# 对@BWT中的各个字符串元素按照字母顺序逐一进行遍历，截取最后一个字母
	print STDERR "BWT matrix:\n";
	my @LastColumn;
	$n=0;
	foreach my $string (sort {$a cmp $b} @BWT){
		print STDERR "\t$n:$string\n";
		push @LastColumn, substr($string, length($string)-1);
		$n++;
	}
	print STDERR "\n";
	print STDERR "BWT output:\n";
	print STDERR "\t".join("", @LastColumn)."\n\n";
	return @LastColumn;
}

my $infile = $ARGV[0];
# 构建BWT的LastColumn
my @BWT_LC = &BWT_Trans($infile);
print "BWT index:".join("", @BWT_LC)."\n";
```

## BWT_revTrans_V1.pl

```perl
#/usr/bin/perl
use strict;
use warnings;

sub FM_index(){
	my ($infile, $SA, $OCC, $LastColumn, $Pre) = @_;
	open IN, "<$infile" or die "$!\n";
	# 解析fasta文件
	my $fa_seq;
	while(<IN>){
		chomp;
		next if(/^>/);
		$fa_seq .= $_;
	}
	close IN;
	# 将以字符串形式组织的序列，打散成数组，并在末尾追加终止标识符
	my @seq = split //,$fa_seq;
	push @seq, "\$";
	print STDERR "Input:\n";
	print STDERR "\t".join("",@seq)."\n\n";
	# 循环出队入队，得到BWT matrix，并在每行末尾加上末尾字符在原始序列的位置索引（1-base）
	my @BWT;
	my @index;
	my $N = $#seq + 1;
	push @BWT, join("", @seq);
	push @index, $N;
	my $n=0;
	print STDERR "Cyclic rotations:\n";
	print STDERR "\t$n:$BWT[$n]\t$index[$n]\n";
	my $tail = pop @seq;
	while($seq[$#seq] ne "\$"){
		$n++;
		$N--;
		unshift @seq, $tail;
		push @BWT, join("",@seq);
		push @index, $N;
		print STDERR "\t$n:$BWT[$n]\t$index[$n]\n";
		$tail = pop @seq;
	}
	print STDERR "\n";
	# 对@BWT中的各个字符串元素按照字母顺序逐一进行遍历，得到FM-index的3个数组和1个散列
	print STDERR "Sort and construct FM-index:\n";
	print STDERR "\tBWTmatrix\tSA\tOCC\n";
	$n=0;
	my %Count;
	## 获得3个数组：SA, OCC, BWT
	foreach my $i (sort {$BWT[$a] cmp $BWT[$b]} 0..$#BWT){
		my $LC = substr($BWT[$i], length($BWT[$i])-1);
		push @$LastColumn, $LC;
		push @$SA, $index[$i];
		$Count{$LC} = defined $Count{$LC} ? $Count{$LC} + 1 : 0;
		push @$OCC, $Count{$LC};
		print STDERR "\t$n:$BWT[$i]\t$$SA[$n]\t$$OCC[$n]\n";
		$n++;
	}
	print STDERR "\n";
	print STDERR "FM-index:\n";
	## 获得1个散列：Pre
	print STDERR "\tPre:";
	$n = 0;
	foreach my $base (sort {$a cmp $b} keys %Count){
		$$Pre{$base} = $n;
		$n += $Count{$base}+1;
		print STDERR "$base=$$Pre{$base};";
	}
	print STDERR "\n";
	print STDERR "\tSA :".join("", @$SA)."\n";
	print STDERR "\tOCC:".join("", @$OCC)."\n";
	print STDERR "\tBWT:".join("", @$LastColumn)."\n\n";
}

sub BWT_reverse_transformation(){
	my ($SA, $OCC, $BWT, $Pre) = @_;
	my @seq;
	my $currentBase = "\$";
	my $currentIndex = 0;
	my $index;
	do{
		unshift @seq, $currentBase;
		# Back_trace，得到当前碱基的前一位碱基的组成
		$currentBase = &Back_trace($currentIndex,$BWT);
		# LC2FC，得到LC中的当前碱基在FC对应的碱基位置
		$currentIndex = &LC2FC($currentBase, $currentIndex, $Pre, $OCC);
	}while($currentBase ne "\$");
	print STDERR "Reverse transformation:".join("", @seq)."\n\n";
}

sub Back_trace(){
	my ($index, $BWT) = @_;
	return $$BWT[$index];
}

sub LC2FC(){
	my ($base, $index, $Pre, $OCC) = @_;
	$index = $$Pre{$base} + $$OCC[$index];
	return $index;
}

my $infile = $ARGV[0];

# 构建FM-index
my (@SA, @OCC, @BWT, %Pre);
&FM_index($infile, \@SA, \@OCC, \@BWT, \%Pre);
print "FM-index:\n";
print "\tPre:";
foreach my $base (sort {$a cmp $b} keys %Pre){
	print STDERR "$base=$Pre{$base};";
}
print "\n";
print "\tSA :".join("", @SA)."\n";
print "\tOCC:".join("", @OCC)."\n";
print "\tBWT:".join("", @BWT)."\n\n";

# 从FM-index还原原始参考序列
&BWT_reverse_transformation(\@SA, \@OCC, \@BWT, \%Pre);
```

## BWT_revTrans_V2.pl

```perl
#/usr/bin/perl
use strict;
use warnings;

sub BWT_Trans(){
	my ($infile) = @_;
	open IN, "<$infile" or die "$!\n";
	# 解析fasta文件
	my $fa_seq;
	while(<IN>){
		chomp;
		next if(/^>/);
		$fa_seq .= $_;
	}
	close IN;
	# 将以字符串形式组织的序列，打散成数组，并在末尾追加终止标识符
	my @seq = split //,$fa_seq;
	push @seq, "\$";
	print STDERR "Input:\n";
	print STDERR "\t".join("",@seq)."\n\n";
	# 循环出队入队，得到BWT matrix
	my @BWT;
	push @BWT, join("", @seq);
	my $n=0;
	print STDERR "Cyclic rotations:\n";
	print STDERR "\t$n:$BWT[$n]\n";
	my $tail = pop @seq;
	while($seq[$#seq] ne "\$"){
		$n++;
		unshift @seq, $tail;
		push @BWT, join("",@seq);
		print STDERR "\t$n:$BWT[$n]\n";
		$tail = pop @seq;
	}
	print STDERR "\n";
	# 对@BWT中的各个字符串元素按照字母顺序逐一进行遍历，截取最后一个字母
	print STDERR "BWT matrix:\n";
	my @LastColumn;
	$n=0;
	foreach my $string (sort {$a cmp $b} @BWT){
		print STDERR "\t$n:$string\n";
		push @LastColumn, substr($string, length($string)-1);
		$n++;
	}
	print STDERR "\n";
	print STDERR "BWT output:\n";
	print STDERR "\t".join("", @LastColumn)."\n\n";
	return @LastColumn;
}

sub BWT_revTrans(){
	my ($BWT) = @_;
	my @Array;
	foreach my $i (0..$#$BWT){
		# 在@Array之前追加新的一列，为@BWT
		foreach my $j (0..$#$BWT){
			$Array[$j] = defined $Array[$j] ? $$BWT[$j].$Array[$j] : $$BWT[$j];
		}
		print STDERR "Roud-$i\n\t";
		print STDERR join("\n\t", @Array)."\n\n\t";
		# 对@Array进行排序
		my @sArray;
		foreach my $string (sort {$a cmp $b} @Array){
			push @sArray, $string;
		}
		@Array = @sArray;
		print STDERR join("\n\t", @Array)."\n";
	}
	return $Array[0];
}

my $infile = $ARGV[0];
# 构建BWT的LastColumn
my @BWT = &BWT_Trans($infile);
my $seq = &BWT_revTrans(\@BWT);

print STDERR "The reconstructed sequence: $seq\n";
```

---

参考资料：

（1）DePristo MA, Banks E, Poplin R, et al. A framework for variation discovery and genotyping using next-generation DNA sequencing data. Nat Genet. 2011 May;43(5):491-8.

（2）GATK官方文档：https://gatk.broadinstitute.org/hc/en-us

（3）GATK Workshop

（4）Heng Li. Mathematical Notes on SAMtools Algorithms. 2010.

（5）VCF格式官方文档：https://samtools.github.io/hts-specs/VCFv4.2.pdf

（6）龙星计划2019

（7）李航. 统计学方法. 

（8）周志华. 机器学习.

（9）Burrows M, Wheeler DJ. Technical report 124. Palo Alto, CA: Digital Equipment Corporation; 1994. A block-sorting lossless data compression algorithm.

（10）Ferragina P, Manzini G. Proceedings of the 41st Symposium on Foundations of Computer Science (FOCS 2000) IEEE Computer Society; 2000. Opportunistic data structures with applications; pp. 390–398.

（11）Li H, Durbin R. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics. 2009 Jul 15;25(14):1754-60.

（12）Langmead B, Trapnell C, Pop M, Salzberg SL. Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biol. 2009;10(3):R25.

（13）Stuart M. Brown等. 第二代测序信息处理. 科学出版社, 2014.

（14）邓俊辉. 数据结构（C++语言版）（第三版）. 清华大学出版社, 2013.
