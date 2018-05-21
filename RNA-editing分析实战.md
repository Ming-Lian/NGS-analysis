<a name="content">目录</a>

[RNA-editing分析实战](#title)
- [获取合适的测试数据](#acquire-test-dataset)








<h1 name="title">RNA-editing分析实战</h1>

<a name="acquire-test-dataset"><h2>获取合适的测试数据 [<sup>目录</sup>](#content)</h2></a>

我们要找到的数据集要满足以下几个条件：
- ribominus建库方法获得的 total RNA-seq
- 数据量大小为双端至少16G，即单端至少8G
- 人的样品，最好能与心脏病相关

到ENCODE数据库中寻找，Data->Matrix，在打开的Experiment Matrix页面中设置过滤条件
> - Assay: total RNA-seq
> - Organism: Homo sapiens
> - Biosample type: tissue
> - Organ: heart

<p align="center"><img src=./picture/InAction-RNA-editing-analysis-acquire-datasets-ENCODE-matched-list.png width=500 /></p>

Experiment list
- heart
	- Homo sapiens heart male adult (34 years) —— [GSE87943](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87943)，一个sample
		`SRP013565` -> `SRR4421689` 19.7Gb
	- Homo sapiens heart female embryo (19 weeks) and female embryo (28 weeks) —— [GSE78567](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78567)，两个sample
		`SRP013565` -> `SRR3192433` 10.9Gb
		`SRP013565` -> `SRR3192434` 13.6Gb
