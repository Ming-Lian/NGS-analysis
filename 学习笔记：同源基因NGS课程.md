<a name="content">目录</a>

[学习笔记：同源基因NGS课程](#title)
- [测序原理](#principle)
	- [FAQ](#faq)
- [测序饱和度评估](#saturation)
- [数据质控与预处理](#qc-and-processing)
	- [数据质控](#qc)
	- [数据预处理](#processing)
- [短序列比对](#short-reads-align)
	- [比对算法与工具](#algorithm-and-tools)
	- [SAM格式](#sam-format)
	- [估计intertsize](#estimate-intertsize)
- [基因差异表达计算](#diff-expression-analysis)
	- [表达定量](#estimate-expression-abundance)
	- [差异统计分析](#diff-statistic)
- [变异检测](#variant-calling)
- [物种组成与丰度分析](#species-composition-and-abundance)
	- [16s高变区测序](#16s-seq)
	- [16s测序分析流程](#16s-workflow)




<h1 name="title">学习笔记：同源基因NGS课程</h1>

<a name="principle"><h2>测序原理 [<sup>目录</sup>](#content)</h2></a>

<p align="center"><img src=./picture/NGS-course-note-1.jpg width=800 /></p>

<p align="center"><img src=./picture/NGS-course-note-2.jpg width=800 /></p>

<a name="faq"><h3>FAQ [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/NGS-course-note-3.jpg width=800 /></p>

<a name="saturation"><h2>测序饱和度评估 [<sup>目录</sup>](#content)</h2></a>

<p align="center">测序乘数与覆盖率之间的关系</p>

<p align="center"><img src=./picture/NGS-course-depth-coverage.png width=700 /></p>

<p align="center"><img src=./picture/NGS-course-note-4-1.jpg width=700 /></p>

<a name="qc-and-processing"><h2>数据质控与预处理 [<sup>目录</sup>](#content)</h2></a>

<a name="qc"><h3>数据质控 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/NGS-course-note-4-2.jpg width=700 /></p>

<p align="center"><img src=./picture/NGS-course-note-5-1.jpg width=700 /></p>

<a name="processing"><h3>数据预处理 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/NGS-course-note-5-2.jpg width=700 /></p>

<a name="short-reads-align"><h2>短序列比对 [<sup>目录</sup>](#content)</h2></a>

<p align="center"><img src=./picture/NGS-course-note-6.jpg width=700 /></p>

<p align="center"><img src=./picture/NGS-course-2-type-coverage.png width=500 /></p>

<a name="algorithm-and-tools"><h3>比对算法与工具 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/NGS-course-note-7.jpg width=700 /></p>

<p align="center"><img src=./picture/NGS-course-note-8-1.jpg width=700 /></p>

<a name="sam-format"><h3>SAM格式 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/NGS-course-note-8-2.jpg width=700 /></p>

第二列（FLAG）标记信息：

<p align="center"><img src=./picture/NGS-course-sam-format-2nd-col-flag.png width=600 /></p>

最后一列比对详细信息：

<p align="center"><img src=./picture/NGS-course-sam-format-last-col.png width=500 /></p>

<a name="estimate-intertsize"><h3>估计intertsize [<sup>目录</sup>](#content)</h3></a>

将pair-end reads 比对到参考基因组上，根据得到的SAM文件中两端reads定位的位置，计算得到intertsize分布图

<table>
<tr>
	<td><img src=./picture/NGS-course-calculate-intertsize.png width=500 /></td>
	<td><img src=./picture/NGS-course-distribution-of-insertsize.png width=400 /></td>
</tr>
</table>

<a name="diff-expression-analysis"><h2>基因差异表达计算 [<sup>目录</sup>](#content)</h2></a>

<a name="estimate-expression-abundance"><h3>表达定量 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/NGS-course-note-9-1.jpg width=700 /></p>

<a name="diff-statistic"><h3>差异统计分析 [<sup>目录</sup>](#content)</h3></a>

p-value的计算：认为P(x)服从泊松分布

<p align="center"><img src=./picture/NGS-course-different-expression-analysis-calculate-pvalue.png width=700 /></p>

p(y|x) 为基因A在两个样本中表达量相等的概率

<p align="center"><img src=./picture/NGS-course-note-10-1.jpg width=700 /></p>

<a name="variant-calling"><h2>变异检测 [<sup>目录</sup>](#content)</h2></a>

<p align="center"><img src=./picture/NGS-course-note-10-2.jpg width=700 /></p>

<p align="center"><img src=./picture/NGS-course-note-11-1.jpg width=700 /></p>

检测插入缺失：

<p align="center"><img src=./picture/NGS-course-detect-indel.png width=500 /></p>

<a name="species-composition-and-abundance"><h2>物种组成与丰度分析 [<sup>目录</sup>](#content)</h2></a>

物种组成检测的案例：

<p align="center"><img src=./picture/NGS-course-species-abundance-analysis-case.png width=700 /></p>

<a name="16s-seq"><h3>16s高变区测序 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/NGS-course-species-abundance-analysis-16s.png width=700 /></p>

<p align="center"><img src=./picture/NGS-course-note-11-2.jpg width=700 /></p>

<a name="16s-workflow"><h3>16s测序分析流程 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/NGS-course-species-abundance-analysis-16s-workflow.png width=700 /></p>


