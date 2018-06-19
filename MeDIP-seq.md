<a name="content">目录</a>

[Analysis pipeline for MeDIP-seq](#title)
- [原理](#principle)
	- [Bisulfite-Sequencing原理](#principle-sequencing)
	- [Bismark原理](#principle-bismark)





<h1 name="title">Analysis pipeline for MeDIP-seq</h1>

<a name="principle"><h3>原理 [<sup>目录</sup>](#content)</h3></a>

<a name="principle-sequencing"><h4>Bisulfite-Sequencing原理 [<sup>目录</sup>](#content)</h4></a>

<p align="center"><img src=./picture/MeDIP-seq-principle-1.png width=800 /></a>

<p align="center"><img src=./picture/MeDIP-seq-principle-2.png width=800 /></a>

<p align="center"><img src=./picture/MeDIP-seq-principle-3.png width=800 /></a>

优点：可以获得单碱基分辨率的甲基化位点

缺点：

- 重硫酸盐处理步骤成本高且费时；
- 苛刻但又是必须的反应条件会导致DNA降解；
- 转换后的基因组的复杂性降低，会给随后的PCR扩增过程的引物设计带来一些限制；
- 增加了与参考基因组比对的难度；
- 无法区别 cytosine，mC 和 hmC；

<a name="principle-bismark"><h4>Bismark原理 [<sup>目录</sup>](#content)</h4></a>

Bisulfite将序列正负链的C全部转换为T，所以也要将基因组序列进行转换。

> - 基因组的正链C->T，才能匹配原正链的reads
> - 基因组的负链C->T相当于正链G->A，才能匹配原负链的reads









参考资料：

(1) [甲基化RRBS流程：原理+bismark使用](https://blog.csdn.net/qq_29300341/article/details/53302961)

(2) Flusberg, B.A. et al. Direct detection of DNA methylation during single-molecule, real-time sequencing. Nat. Methods 7, 461–465 (2010).
