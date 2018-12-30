<a name="content">目录</a>

[群体遗传学知识点](#title)
- [Fst及计算方法](#fst)






<h1 name="title">群体遗传学知识点</h1>

<a name="fst"><h2>Fst及计算方法 [<sup>目录</sup>](#content)</h2></a>

Fst是什么？

> Fst：群体间遗传分化指数，是种群分化和遗传距离的一种衡量方法，分化指数越大，差异越大。适用于亚群体间多样性的比较。
> 
> 用于衡量种群分化程度，取值从0到1，为0则认为两个种群间是随机交配的，基因型完全相似；为1则表示是完全隔离的，完全不相似。它往往从基因的多样性来估计，比如SNP或者microsatellites(串联重复序列一种，长度小于等于10bp)。是一种以哈温平衡为前提的种群遗传学统计方法。

Fst的计算公式如下：

<p align="center"><img src=./picture/PopulationGenetics-Fst-formula.png width=200 /></p>

Hs：亚群体中的平均杂合度
Ht：复合群体中的平均杂合度

在遗传学中，F一词通常代表“近亲繁殖”，它倾向于减少群体中的遗传变异。遗传变异可以用杂合度来衡量，所以F一般表示群体中杂合性的减少。 FST是与它们所属的总群体相比，亚群体中杂合性的减少量。

**如何计算Fst？**

以一个例子进行说明：

> 基因SLC24A5是黑色素表达途径的关键部分，其导致皮肤和毛发色素沉着。与欧洲较轻的皮肤色素密切相关的SNP是rs1426654。 SNP有两个等位基因A和G，其中G与轻度皮肤相关，在犹他州的欧裔美国人中，频率为100％。
> 
> 美洲印第安人与美国印第安人混血儿的SNP在频率上有所不同。
> - 墨西哥的样本有38％A和62％G;
> - 在波多黎各，频率分别为59％A和41％G；
> - 查尔斯顿的非裔美国人样本中有19％A和81％G；
> 
> 这个例子中的FST是什么？

<p align="center"><img src=./picture/PopulationGenetics-Fst-example.png width=800 /></p>

Fst值的计算可以使用VCFtools实现

```
# 对每一个SNP变异位点进行计算

vcftools --vcf test.vcf --weir-fst-pop 1_population.txt --weir-fst-pop 2_population.txt  --out p_1_2—single

# 按照区域来计算

vcftools --vcf test.vcf --weir-fst-pop 1_population.txt --weir-fst-pop 2_population.txt  --out p_1_2_bin --fst-window-size 500000 --fst-window-step 50000
```
参数说明：

> - `--vcf`：SNP calling 过滤后生成的vcf 文件；
> - `--weir-fst-pop`：群体组成文件，包含同一个群体中所有个体，一般每行一个个体。个体名字要和vcf的名字对应；
> - `--out`：生成结果的prefix；
> - `--fst-window-size`与`--fst-window-step`：在按照区域计算Fst时，指定计算窗口和步长

<p align="center"><img src=./picture/PopulationGenetics-Fst-plot-1.png width=800 /></p>

<p align="center">对每个SNP位点计算Fst的散点分布图</p>

<p align="center"><img src=./picture/PopulationGenetics-Fst-plot-2.png width=800 /></p>

<p align="center">对每个区块计算Fst的散点分布图</p>

---

参考资料：

(1) [【简书】Fst的计算原理与实战](https://www.jianshu.com/p/b73a8d6233be)
