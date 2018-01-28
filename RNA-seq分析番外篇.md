<p name="content">目录</p>

[RNA-seq分析番外篇](#title)
- [过滤细菌库](#filt-baclib)
	- [与细菌库比对](#filt-map)
	- [过滤出mapping上细菌库的序列](#filt-bac)
	- [过滤掉mapping上细菌库的序列](#remove-bac)
- [饱和度分析](#saturation)
	- [按照比例](#percent)
	- [按照测序深度](#depth)


<h1 name="title">RNA-seq分析番外篇</h1>

<a name="filt-baclib"><h3>过滤细菌库 [<sup>目录</sup>](#content)</h3></a>

<a name="filt-map"><h4>与细菌库比对（进行单端mapping)</h4></a>

```
bwa aln -t 10 -f $wd/$map_dir/${1}_BacLib_1.sai $wd/Ref/NCBI_bacteria/referance.fasta $wd/$cleandata_dir/${1}_Clean_Data1.${mate_str}fq 1>$wd/$map_dir/${1}_BacLib_1.sai.log 2>&1
echo "[filt] Finish searching SA coordinate for ${1}_Clean_Data1.${mate_str}fq"
bwa aln -t 10 -f $wd/$map_dir/${1}_BacLib_2.sai $wd/Ref/NCBI_bacteria/referance.fasta \
	$wd/$cleandata_dir/${1}_Clean_Data2.${mate_str}fq 1>$wd/$map_dir/${1}_BacLib_2.sai.log 2>&1
echo "[filt] Finish searching SA coordinate for ${1}_Clean_Data2.${mate_str}fq"
bwa samse -f $wd/$map_dir/${1}_BacLib_1.sam $wd/Ref/NCBI_bacteria/referance.fasta $wd/$map_dir/${1}_BacLib_1.sai \
	$wd/$cleandata_dir/${1}_Clean_Data1.${mate_str}fq 1>$wd/$map_dir/${1}_BacLib_1.sam.log 2>&1
echo "[filt] Finish single-end mapping for ${1}_Clean_Data1.${mate_str}fq"
bwa samse -f $wd/$map_dir/${1}_BacLib_2.sam $wd/Ref/NCBI_bacteria/referance.fasta $wd/$map_dir/${1}_BacLib_2.sai \
	$wd/$cleandata_dir/${1}_Clean_Data2.${mate_str}fq 1>$wd/$map_dir/${1}_BacLib_2.sam.log 2>&1
echo "[filt] Finish single-end mapping for ${1}_Clean_Data2.${mate_str}fq"
```

<a name="filt-bac"><h4>过滤<font color="red">出</font>mapping上细菌库的序列</h4></a>

```
perl -ane 'chomp;next if (/^\@/);if ($F[2] ne "*"){print "$_\n"}' $wd/$map_dir/${1}_BacLib_1.sam >$wd/$map_dir/${1}_BacLib_1.sam.filt
perl -ane 'chomp;next if (/^\@/);if ($F[2] ne "*"){print "$_\n"}' $wd/$map_dir/${1}_BacLib_2.sam >$wd/$map_dir/${1}_BacLib_2.sam.filt
```
<a name="remove-bac"><h4>过滤<font color="red">掉</font>mapping上细菌库的序列</h4></a>
```
perl $wd/perlscript/sel_seq_for_Hiseq3000.pl $wd/$map_dir/${1}_BacLib_1.sam.filt $wd/$map_dir/${1}_BacLib_2.sam.filt \ 			$wd/$cleandata_dir/${1}_Clean_Data1.${mate_str}fq $wd/$cleandata_dir/${1}_Clean_Data2.${mate_str}fq \
	$wd/$cleandata_dir/${1}_Clean_Data1.filtBacLib.fq $wd/$cleandata_dir/${1}_Clean_Data2.filtBacLib.fq 
echo "[filt] Finish filt bac contaminate"
```

<a name="saturation"><h3>饱和度分析 [<sup>目录</sup>](#content)</h3></a>

```
bwa aln -t 10 -f $wd/$map_dir/cds_${sample}_1.sai $wd/$ref_cds $wd/$cleandata_dir/${sample}_Clean_Data1.filtBacLib.fq 1>$wd/$map_dir/cds_${sample}_1.sai.log 2>&1
echo "[saturation] Finish search SA coordinate for ${sample}_Clean_Data1.filtBacLib.fq"
bwa aln -t 10 -f $wd/$map_dir/cds_${sample}_2.sai $wd/$ref_cds $wd/$cleandata_dir/${sample}_Clean_Data2.filtBacLib.fq 1>$wd/$map_dir/cds_${sample}_2.sai.log 2>&1
echo "[saturation] Finish search SA coordinate for ${sample}_Clean_Data2.filtBacLib.fq"
bwa sampe -f $wd/$map_dir/cds_${sample}.sam $wd/$ref_cds $wd/$map_dir/cds_${sample}_1.sai \
	$wd/$map_dir/cds_${sample}_2.sai $wd/$cleandata_dir/${sample}_Clean_Data1.filtBacLib.fq $wd/$cleandata_dir/${sample}_Clean_Data2.filtBacLib.fq 1>$wd/$map_dir/cds_${sample}.sam.log 2>&1
echo "[saturation] Finish pair-end mapping for $sample"
```

<a name="percent"><h4>按照比例</h4></a>

```
for i in 5 10 15 20 30 40 50 60 70 80 90
do
	export i
	perl -ne 'chomp;next if (/^\@/);if (rand()<0.01*$ENV{"i"}){print "$_\n";}' $wd/$map_dir/cds_${sample}.sam | cut -f 3 | sort | uniq | wc -l >>$wd/$satu_dir/stat/${sample}.txt
	echo "[saturation] Analyse saturation for $sample: ${i}%"
done 
rm $wd/$map_dir/cds_${sample}_[12].sa[im]
rm $wd/$map_dir/cds_${sample}.sam
```

<a name="depth"><h4>按照测序深度（一个外显子组大小为50M）</h4></a>

```
total_reads=`perl -ne 'chomp;next if (/^\@/);print "$_\n"' $wd/$map_dir/cds_${sample}.sam | wc -l | perl -ane 'chomp;print "$F[0]"'`
export total_reads
exome_length=50000000
reads_length=150
depth1_reads=$[$exome_length/$reads_length]
export depth1_reads
max_depth=$[$total_reads*$reads_length/$exome_length]
if [ -f $wd/$satu_dir/stat/${sample}_depth.txt ]
then
	rm $wd/$satu_dir/stat/${sample}_depth.txt
fi
for ((i=5;i<=max_depth;i+=5))
do
	current_depth=$i
	export current_depth
	echo -ne "$current_depth\t" >>$wd/$satu_dir/stat/${sample}_depth.txt
	perl -ne 'chomp;next if (/^\@/);if (rand()<($ENV{"current_depth"}*$ENV{"depth1_reads"}/$ENV{"total_reads"})){print "$_\n";}' $wd/$map_dir/cds_${sample}.sam | cut -f 3 | sort | uniq | wc -l >>$wd/$satu_dir/stat/${sample}_depth.txt
	echo "[saturation] Analyse saturation for $sample: ${current_depth}X"
done
```
