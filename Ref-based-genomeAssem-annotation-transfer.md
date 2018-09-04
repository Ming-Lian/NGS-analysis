<a name="content">目录</a>

[基于重测序得到Population/Strain 的consensus序列和对应的注释](#title)
- [方法一：基于variant calling](#based-on-variant-calling)
- [方法二：基于genome comparison](#based-on-genome-comparison)


<h1 name="title">基于重测序得到Population/Strain 的consensus序列和对应的注释</h1>

<a name="based-on-variant-calling"><h2>方法一：基于variant calling [<sup>目录</sup>](#content)</h2></a>

<p align="center"><img src=./picture/Ref-based-genomeAssem-annotation-transfer-strategy-1-workflow.png width=600/></p>

<a name="based-on-genome-comparsion"><h2>方法二：基于genome comparison [<sup>目录</sup>](#content)</h2></a>

该方法需要先进行基因组拼接，然后将拼接完的基因组与标准的参考基因组进行比较，得到相对于参考基因组的SNP/Indel，从而得到对于参考基因组的坐标转换文件（genome coordinate conversion file）：

| 列序号 | 说明 |
|:---:|:---|
| 1 | TargetChrom |
| 2 | TargetStart |
| 3 | TargetEnd |
| 4 | QueryStart |
| 5 | QueryEnd |
| 6 | QueryChrom |
| 7 | QueryStrand |


然后基于坐标对应关系，对参考基因组的注释文件进行迁移修改，就得到了新拼接基因组的注释文件

<p align="center"><img src=./picture/Ref-based-genomeAssem-annotation-transfer-strategy-2-workflow.png width=600/></p>

---

参考资料：

Wanfei Liu et al. RGAAT - Reference Based Genome Assembly and Annotation Tool for New Genome and Known Genome Upgrade. On publish.
