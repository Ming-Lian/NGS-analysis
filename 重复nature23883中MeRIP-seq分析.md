<a name="content">目录</a>

[重复nature23883中MeRIP-seq分析](#title)
- [研究思路与成果](#ideas-and-result)
- [科普：MeRIP–seq和miCLIP–seq技术](#introduct-technology)







<h1 name="title">重复nature23883中MeRIP-seq分析</h1>

<a name="ideas"><h3>研究思路 [<sup>目录</sup>]()</h3></a>

<p align="left"> <li>基因敲低技术发现了m6A对斑马鱼胚胎的早期发育有重要作用，并选择了该实验室擅长的斑马鱼造血系统作为研究对象。</li></p>

<p align="left"> <li>通过降低m6A的甲基化酶 （mettl3）在细胞中的水平，然后发现了4000多个mRNA的m6A 水平降低。这些target基因中，表达升高的mRNA富集在了Notch信号通路中</li></p>

**mettl3 blocked** => **4000多个mRNA的m6A ↓** => **表达升高部分富集于Notch信号通路**

即： **mettl3 blocked** => **Notch信号激活**

说明： **mettl3** => **m6A特异性修饰通路的关键基因** => **Notch信号通路低水平表达** => **EMT**

<p align="left"> <li>从这些上调基因中，选择了Notch通路中重要受体蛋白notch1a作为研究对象进行研究</li></p>

通过一系列实验证明mettl3 KO的胚胎中：

**mettl3 KO** => **notch1a m6A ↓** => **notch1a表达水平升高** => **内皮无法正常向造血细胞转化**

<p align="left"> <li>接着搞清楚notch1a与EMT之间的机制</li></p>

结合之前发表的研究，m6A结合蛋白YTHDF2参与RNA的降解，通过RIP-seq技术，发现ythdf2能够结合notch1a，

<a name="introduct-technology"><h3>科普：MeRIP–seq和miCLIP–seq技术 [<sup>目录</sup>]()</h3></a>

**MeRIP–seq**: 通过识别m6A修饰的抗体来富集含有m6A的RNA片段，然后对这些片段进行建库和高通量测序，并通过和input比较，可以得到m6A peak的修饰强度。缺点：分辨率有限，不能够准确定位到底是哪个碱基发生化学修饰

<img src=./picture/MeRIP-seq.jpg width=500 />

**miCLIP**: 也是用抗体富集RNA片段，通过鉴定紫外照射产生的突变位点以及RT-PCR阻断位点，可实现单碱基精度

<table>
<tr>
<td>
<p>1. RNA extraction from HEK293 cell lines (human) or liver nuclei (mouse) using Trizol.</p>
<br>
<br>
<p>2. Fragmentation of RNA to 30–130 nucleotide lengths, and incubation with anti-m6A antibody.</p>
<br>
<br>
<br>
<p>3. UV cross-linking of RNA to the bound antibody.</p>
<br>
<br>
<p>4. Recovery of antibody-RNA complexes with protein A/G affinity purification, SDS-PAGE and nitrocellulose membrane transfer.</p>
<br>
<br>
<p>5. Adapter ligation, and release of RNA with proteinase K.</p>
<br>
<p>6. Reverse transcription of RNA to cDNA, PCR amplification and sequencing.</p>
<br>
<p>7. Identification of C-T transitions or truncations and alignment against known genomic sequences. Mapping and annotation of these binding sites, identified as m6A/m6Am residues, to the transcriptome.</p>
</td>
<td>
<img src=./picture/MeRIP-seq-inAction-miCLIP.jpg width=600 />
</td>
</tr>
</table>




参考资料：

(1) [探索m6A修饰在早期胚胎造血过程中的作用 | 专访基因组所Nature 共同第一作者陈宇晟](http://mp.weixin.qq.com/s/vIZXp60JM838-qE4vmszcw)
