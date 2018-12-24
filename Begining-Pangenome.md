<a name="content">目录</a>

[泛基因组入门](#title)
- [概念](#introduction)
- [分析流程](#analysis-workflow)
	- [组装方法](#assembly-approaches)
- [应用](#applyication)


<h1 name="title">泛基因组入门</h1>

<a name="introduction"><h2>概念 [<sup>目录</sup>](#content)</h2></a>

泛基因组即某一物种全部基因的总称

<p align="center"><img src=./picture/Pangenome-classification.jpg width=800 /></p>

核心基因：在所有动植物品系或者菌株中都存在的基因；

可变基因组：在1个以及1个以上的动植物品系或者菌株中存在的基因

品系或者菌株特有基因：某个基因，仅存在某一个动植物品系或者菌株中

另外，结构变异中的存在/缺失变化(presnece/absence variation)是泛基因组的重点研究对象，因为可变基因组可能就是使个体产生不同性状（抗病性，抗寒性等）的原因。

开展泛基因组测序的意义：

> 在漫长的进化过程中，由于地域因素，环境因素等的影响，每个个体都形成了极其特别的遗传性状，单一个体的基因组已经不能涵盖这个物种的所有遗传信息，另外一个原因，由于基因测序变得更加廉价，为近年来火爆的泛基因组的研究提供了可能性

<a name="analysis-workflow"><h2>分析流程 [<sup>目录</sup>](#content)</h2></a>

<p align="center"><img src=./picture/Pangenome-analysis-workflow.jpg width=500 /></p>

<a name="assembly-approaches"><h3>组装方法 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/Pangenome-analysis-workflow-assembly-approaches.jpg width=800 /></p>

1. 重头组装 (De novo assembly) 的方法

	<p align="center"><img src=./picture/Pangenome-analysis-workflow-assembly-approaches-1.jpg width=800 /></p>

	分别对多个个体进行，分别的De novo assembly，然后将所得的每个个体的新组装的序列与参考序列reference基因组进行比对，找出比对不上的区域，再进行进一步的assembly，然后注释。

	需要更多的电脑资源，因为需要对每一个个体进行分别进行重头组装，然后还需要全基因组比对。该方法比较适合基因组相对较小的植物。

2. 迭代组装的方法

	<p align="center"><img src=./picture/Pangenome-analysis-workflow-assembly-approaches-2.jpg width=400 /></p>

	相当于一种迭代的方式，分别将每一个材料的reads比对到参考基因组中，然后找出没有比对上的部分进行组装，得到新的基因序列进而扩展原有的参考序列。一步一步这样迭代，直到所有的种系都处理完。最后建立起的泛基因组，再进行注释

	相对需要更少电脑资源，比较适合构建基因size相对较大的植物泛基因组，但是可能会产生更多的小片段

<a name="applyication"><h2>应用 [<sup>目录</sup>](#content)</h2></a>

- 选择**不同亚种**材料进行泛基因组测序，可以研究物种的**起源及演化**等重要生物学问题；

- 选择**野生种和栽培种**等不同特性的种质资源进行泛基因组测序，可以**发掘重要性状相关的基因资源**，为科学育种提供指导；

- 选择**不同生态地理类型的种质资源**进行泛基因组测序，可以开展物种的**适应性进化**，外来物种入侵性等热门科学问题的研究，为分子生态学等学科提供新的研究手段；




---

参考资料：

(1) [【生信菜鸟团】生信杂谈：认识泛基因组测序](https://mp.weixin.qq.com/s/4rUetrV9V1CZGFmBGqkftg)
