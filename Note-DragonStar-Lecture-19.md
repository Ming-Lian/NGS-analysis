<a name='content'>目录</a>

[2019年龙星计划学习笔记](#title)
- [Day-1](#day-1)
    - [1. Data-format & Variant-call](#day-1-part1)
        - [1.1. 不同文件格式采用的坐标系统](#day-1-1)
        - [1.2. Coverage的两个含义](#day-1-2)
        - [1.3. 测序量与覆盖率的关系：泊松分布建模](#day-1-3)
        - [1.4. Variant calling的基本思路](#day-1-4)
    - [2. Genetic technology](#day-1-part2)
        - [2.1. 基因组变异与影响](#day-1-5)
        - [2.2. 测序技术](#day-1-6)
            - [2.2.1. 华大测序仪及其测序原理](#day-1-6-1)
            - [2.2.2. SMRT测序](#day-1-6-2)
- [Day-2](#day-2)
    - [1. Alignment](#day-2-alignment)
        - [1.1. Blast](#day-2-alignment-1)
        - [1.2. BWT](#day-2-alignment-2)
        - [1.3. Minimap2](#day-2-alignment-3)
    - [2. Genome Assembly](#day-2-assembly)
        - [2.1. 拼接算法](#day-2-assembly-1)
            - [2.1.1. OLC](#day-2-assembly-1-1)
            - [2.1.2. de Bruijn Graph](#day-2-assembly-1-2)
        - [2.2. Lander-Waterman statistcs：建模测序覆盖率和gaps数](#day-2-assembly-2)

<h1 name="title">2019年龙星计划学习笔记</h1>

<a name="day-1"><h2>Day-1 [<sup>目录</sup>](#content)</h2></a>

<a name="day-1-part1"><h2>1. Data-format & Variant-call [<sup>目录</sup>](#content)</h2></a>

<a name="day-1-1"><h3>1.1. 不同文件格式采用的坐标系统 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/Note-DragonStarLecture-01_1.jpg width=800 /></p>

<a name="day-1-2"><h3>1.2. Coverage的两个含义 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/Note-DragonStarLecture-01_2.jpg width=800 /></p>

<a name="day-1-3"><h3>1.3. 测序量与覆盖率的关系：泊松分布建模 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/Note-DragonStarLecture-02.jpg width=800 /></p>

<p align="center"><img src=./picture/Note-DragonStarLecture-03.jpg width=800 /></p>

<a name="day-1-4"><h3>1.4. Variant calling的基本思路 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/Note-DragonStarLecture-04.jpg width=800 /></p>

<p align="center"><img src=./picture/Note-DragonStarLecture-05.jpg width=800 /></p>

<p align="center"><img src=./picture/Note-DragonStarLecture-06.jpg width=800 /></p>

<p align="center"><img src=./picture/Note-DragonStarLecture-07_1.jpg width=800 /></p>

<a name="day-1-part2"><h2>2. Genetic technology [<sup>目录</sup>](#content)</h2></a>

<a name="day-1-5"><h3>2.1. 基因组变异与影响 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/Note-DragonStarLecture-07_2.jpg width=800 /></p>

<p align="center"><img src=./picture/Note-DragonStarLecture-08.jpg width=800 /></p>

<a name="day-1-6"><h3>2.2. 测序技术 [<sup>目录</sup>](#content)</h3></a>

<a name="day-1-6-1"><h4>2.2.1. 华大测序仪及其测序原理 [<sup>目录</sup>](#content)</h4></a>

<p align="center"><img src=./picture/Note-DragonStarLecture-09.jpg width=800 /></p>

<a name="day-1-6-2"><h4>2.2.2. SMRT测序 [<sup>目录</sup>](#content)</h4></a>

<p align="center"><img src=./picture/Note-DragonStarLecture-10.jpg width=800 /></p>

<a name="day-2"><h2>Day-2 [<sup>目录</sup>](#content)</h2></a>

<a name="day-2-alignment"><h2>1. Alignment [<sup>目录</sup>](#content)</h2></a>

<a name="day-2-alignment-1"><h3>1.1. Blast [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/Note-DragonStarLecture-11.jpg width=800 /></p>

<p align="center"><img src=./picture/Note-DragonStarLecture-12.jpg width=800 /></p>

<p align="center"><img src=./picture/Note-DragonStarLecture-13_1.jpg width=800 /></p>

<a name="day-2-alignment-2"><h3>1.2. BWT [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/Note-DragonStarLecture-13_2.jpg width=800 /></p>

<p align="center"><img src=./picture/Note-DragonStarLecture-14.jpg width=800 /></p>

<p align="center"><img src=./picture/Note-DragonStarLecture-15_1.jpg width=800 /></p>

<a name="day-2-alignment-3"><h3>1.3. Minimap2 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/Note-DragonStarLecture-15_2.jpg width=800 /></p>

<p align="center"><img src=./picture/Note-DragonStarLecture-16.jpg width=800 /></p>

<p align="center"><img src=./picture/Note-DragonStarLecture-17_1.jpg width=800 /></p>

<a name="day-2-assembly"><h2>2. Genome Assembly [<sup>目录</sup>](#content)</h2></a>

<p align="center"><img src=./picture/Note-DragonStarLecture-17_2.jpg width=800 /></p>

<a name="day-2-assembly-1"><h3>2.1. 拼接算法 [<sup>目录</sup>](#content)</h3></a>

<a name="day-2-assembly-1-1"><h4>2.1.1. OLC [<sup>目录</sup>](#content)</h4></a>

<p align="center"><img src=./picture/Note-DragonStarLecture-18.jpg width=800 /></p>

<p align="center"><img src=./picture/Note-DragonStarLecture-19_1.jpg width=800 /></p>

<a name="day-2-assembly-1-2"><h4>2.1.2. de Bruijn Graph [<sup>目录</sup>](#content)</h4></a>

<p align="center"><img src=./picture/Note-DragonStarLecture-19_2.jpg width=800 /></p>

<p align="center"><img src=./picture/Note-DragonStarLecture-20.jpg width=800 /></p>

<a name="day-2-assembly-2"><h3>2.2. Lander-Waterman statistcs：建模测序覆盖率和gaps数 [<sup>目录</sup>](#content)</h3></a>

<p align="center"><img src=./picture/Note-DragonStarLecture-21.jpg width=800 /></p>

<p align="center"><img src=./picture/Note-DragonStarLecture-22.jpg width=800 /></p>

---

参考资料：

[2019-DragonStar Bioinformatics ppt](https://github.com/WGLab/dragonstar2019)
