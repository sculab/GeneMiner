# GeneMiner

​	感谢您选择GeneMiner(GM)。本文档将帮助您学会如何使用GeneMiner。请注意，GeneMiner仍处于测试阶段，因此本文档可能会在未来的版本升级中更新。

​	如果您有任何问题，可以联系**xiepulin@stu.scu.edu.cn**  或者 15162176893@163.com

[TOC]

## 1. 简介

​	GeneMiner(Gene-Miner，基因矿工)是一款新开发的用于从下一代测序数据中挖掘系统发育标记的软件。GeneMiner 应用场景广阔。可以从浅层基因组中提取叶绿体全基因组，全部或者部分线粒体基因组，以及核基因组中高度重复区 (如rDNA)等；也可以从转录组中提取大量的系统发育标记，如被子植物353单拷贝基因。	GeneMiner运行速度很快。因为软件选择了数据量大小合适的原始数据，避免算力浪费；采用了后缀树算法,加快reads过滤速度;选择了适合于短reads拼接的软件（minia）；支持多线程并行，充分调用计算机资源等一系列方法，大幅度增加软件运行速度，让用户能在较短的时间挖掘大量系统发育标记。GeneMiner的结果十分准确。我们通过实验和自展检测的方法对结果进行校验，准确率均超过99%。GeneMiner操作简单，使用方便。用户可以在Windows、Mac和Linux平台上使用，在未来的版本中，我们也将推出对用户更加友好的图形界面版本。最后，GeneMiner在降低测序成本，扩大基因选择方面也有显著优势。

## 2.下载

​		https://github.com/sculab

​        后续补充

## 3. 安装

For Linux users:

```shell
#(1)解压
tar -zxvf  GeneMiner.tar.gz

#(2)校验环境与安装依赖
cd GeneMiner
python setup.py 

如果这上一步安装依赖失败，可以自己安装
pip3 install -r requirements --user #非管理员
sudo pip3 install -r requirements   #管理员

#(3)将GeneMiner写入环境变量
#打开系统环境变量配置文件
vim ~/.bashrc
#打开后在文件最后面添加以下语句
export PATH="your/GeneMiner/path:$PATH" 
#保存退出，然后在bash中执行source命令使其生效
source ~/.bashrc

#(4)检测是否配置成功
GeneMiner -h
```

**For Windows users:** 

**For Mac users:** 



## 4.快速入门

​		在对你的测序数据进行挖掘之前，我们建议你先了解自己的测序数据。了解测序方式、测序深度、质量、数据量大小等等。我们的软件主要适用于Illumina平台所返回的二代测序数据。经大量的数据测试验证，GeneMiner能从10-20x的二代测序数据中挖掘到nuclear_gene、cp-gene以及mito-gene等，也可以从转录组数据中挖掘出部分单拷贝基因。

（1）挖掘单个目标基因

```shell
 GeneMiner -1 data.1.fq.gz -2 data.2.fq.gz -rt ref.fa
```

值得注意的是：每一个fasta格式中只能存放同一种基因。如以下样例：

```
>species1 GeneA
ATCGATCG
>species2 GeneA
ATCGATCC
>species3 GeneA
ATTGATCC
```

（2）挖掘多个目标基因

可以将多个fasta格式的文件放在一个文件夹中，便能够批量挖掘基因

```shell
mkdir ref

mv ref1.fa ref2.fa ref3.fa ref

GeneMiner -1 data.1.fq.gz -2 data.2.fq.gz -rt ref
```

（3）批量挖掘不同类型的基因

```shell
GeneMiner -1 data.1.fq.gz -2 data.2.fq.gz -rn nuclear_ref -rmito mito_ref -rcp cp_re
```

（4)检测结果正确性

```shell
GeneMiner -1 data.1.fq.gz -2 data.2.fq.gz -rn ref --check True --bootstrap_number 10
```

值得注意的是：采用自展检测的方法对结果正确性进行评估，同时计算耗时将大幅度增加



## 5.详细介绍

### 5.1参数解读

使用GeneMiner -h ，计算机将显示所有选项，后面我将详细介绍它们的使用方法和一些小窍门

```shell
GeneMiner -h
usage: GeneMiner <-1 -2|-s|-12>  <-rn|rcp|-rmito|rt>  [options]

GeneMiner: a software for extracting phylogenetic markers from skimming genome
Version: 1.0
Copyright (C) 2021 Pu-lin Xie
Please contact <xiepulin@stu.edu.scu.cn>, if you have any bugs or questions

optional arguments:
  -h, --help            show this help message and exit

Basic option:
  -1                    forward paired-end reads,support fastq/fastq.gz/fastq.bz2
  -2                    reverse paired-end reads,support fastq/fastq.gz/fastq.bz2
  -12                   interlaced forward and reverse paired-end reads,support fastq/fastq.gz/fastq.bz2
  -s , --single         single-read, support fastq/fastq.gz/fastq.bz2
  -o , --out            Specify the result folder [default='auto']
  -rcp <file|dir>       reference of chloroplast genome,only support GenBank-format
  -rmito <file|dir>     reference of mitochondrial genome,only support GenBank-format
  -rn <file|dir>        reference of nuclear_genes,only support fasta-format
  -rt <file|dir>        reference of target_genes,only support fasta-format

Advanced option:
  -n , --number         The number of rows of raw data from skimming genomes,default=1000000
  -k , --kmer           size of a kmer  [default =31]
  -max                  gene maximum length  [default =5000]
  -min                  gene minimum length  [default =300]
  -t , --thread         Specify the number of threads you want to run [default='auto']
  -b , --boundary       extend the length to both sides of the gene while extracting  genes from  Genebank file [default=75]

Bootstrap option:
  --check [False,True]  Evaluate the accuracy of the results by using bootstrap method.
                        If you use this parameter, the computation will be greatly increased.default=False
  --bootstrap_number    Specify the bootstrap number . The number ranges from 1 to 1024

```



#### 5.1.1 基础参数

```shell
-1				
双末端测序数据的正向数据  支持fastq/fastq.gz/fastq.bz2格式，务必保留正确的文件名后缀
example: GeneMiner -1 data.1.fq.gz -2 data.2.fq.gz -rt ref
    
-2   			
双末端测序数据的反向数据 支持fastq/fastq.gz/fastq.bz2格式，务必保留正确的文件名后缀
 example: GeneMiner -1 data.1.fq.gz -2 data.2.fq.gz -rt ref

-12  			
交错合并的双末端测序数据 支持fastq/fastq.gz/fastq.bz2格式，务必保留正确的文件名后缀
example: GeneMiner -12 data.fq.gz -rt ref

-s   			
单端测序数据 支持fastq/fastq.gz/fastq.bz2格式，务必保留正确的文件名后缀
example: GeneMiner -s data.fq.gz -rt ref

-rn  <file|dir>       	
核基因参考序列，仅支持fasta格式。
可以输入一个fasta文件，该文件里既可以包含一个种，也可以包含多个种；还可以将多个fasta文件放入一个文件夹下。
值得注意的是，同一个fasta格式的文件里，只能存放同一个基因的数据，例子如下：
>species_a ITS
ACGT
>species_b ITS 
AATT
>species_c ITS
ATAT
t4 species_d ITS
CCGT

-rcp 	<file|dir>	
叶绿体参考基因组，只支持GenBank的格式。
可以输入一个GenBank文件，该文件里既可以包含一个种，也可以包含多个种；还可以将多个Genebank文件放入一个文件夹下。 
值得注意的是，用户可以修改GeneBank文件，仅保留自己感兴趣的部分基因。同时在选择参考基因组的时候尽量选择近缘参考基因组
example: GeneMiner -1 data.1.fq.gz -2 data.2.fq.gz -rcp ref

-rmito    <file|dir>	
线粒体参考基因组，只支持GenBank的格式。具体用法同rcp

-rt     <file|dir>  		
目标基因参考序列，仅支持fasta格式。可用于寻找自己感兴趣的基因，自然包括了线粒体，叶绿体，核基因。
值得注意的是，-rt的功能囊括了-rn -rcp -rmito。之所以将其区分开来，是为了方便用户使用以及后期功能模块添加
example: GeneMiner -1 data.1.fq.gz -2 data.2.fq.gz -rt ref


-o , --out            
指定输出文件夹。如果不指定输出文件夹，GeneMiner将自动使用 ''GM+时间戳'' 作为输出文件夹名字
example: GeneMiner -1 data.1.fq.gz -2 data.2.fq.gz -rt ref -o GeneMiner_out

     
```



#### **5.1.2 高级参数**

```shell
-n , --number     
输入原始数据量行数，默认为100w行。
经过测试，仅选择原始数据的一部分就能得到很好的结果，却能够大大减少软件运行时间
对于读长为150bp的二代测序数据，100w行的数据量大概在80~100M之间。对于linux用户，可以使用cat,zcat,wc等命令查看原始数据行数
如:cat data.fq|wc -l  or zcat data.fq.gz|wc -l
值得注意的是:并非选择的数据量越大越好-n在100w~1000w之间效果较佳
example:
GeneMiner -1 data.1.fq.gz -2 data.2.fq.gz -rt ref -n 2000000

-k , --kmer   
指定k-mer长度，k-mer是de Bruijn图中节点的长度。kmer的取值严重依赖于数据集。
默认大小为31，kmer的取值范围在15~127之间
example:
GeneMiner -1 data.1.fq.gz -2 data.2.fq.gz -rt ref -n 2000000 -k 43

-max 
指定挖掘基因的最大长度，默认为5000bp.可结合用于-rcp/-rmito参数
example:
GeneMiner -1 data.1.fq.gz -2 data.2.fq.gz -rcp ref -n 2000000 -k 43 -max 3000

-min    
指定挖掘基因的最大长度，默认为300bp.可结合用于-rcp/-rmito参数
example:
GeneMiner -1 data.1.fq.gz -2 data.2.fq.gz -rcp ref -n 2000000 -k 43 -max 3000 -min 500

-b , --boundary   
指定软边界的长度。
当沿着挖掘出的基因的两侧延申的时候，延申长度越长，准确率越低。但这种下降不是断崖式的，而是在一定长度的缓冲区内逐渐下滑。我们将保留一段长度的缓冲区，并将这段缓冲区称为软边界。默认为0.5*reads长度，取值范围在0~200之间
example:
GeneMiner -1 data.1.fq.gz -2 data.2.fq.gz -rcp ref -n 2000000 -k 43 -max 3000 -min 500 -b 75

-t , --thread  
example:
指定线程数量 ，软件默认根据计算机性能选择合适的线程数量
GeneMiner -1 data.1.fq.gz -2 data.2.fq.gz -rcp ref -n 2000000 -k 43 -max 3000 -min 500 -b 75 -t 8

```

#### 5.1.3自展检测参数

```shell
--check [False,True] 
通过自展检测的方法评估结果的正确性。
example:
GeneMiner -1 data.1.fq.gz -2 data.2.fq.gz -rt ref --check True

--bootstrap_number  
指定自展检测次数，默认为10次，取值范围为1~1024.取值越大，计算量越大，耗时越久
GeneMiner -1 data.1.fq.gz -2 data.2.fq.gz -rt ref --check True --bootstrap_number 20
```

### 5.2 结果解读

#### 5.2.1目录结构

```shell
GeneMiner将产生大量的文件，主要分为三大部分：
第一部分
合适大小的原始数据 ，forward.fq 和 reverse.fq
第二部分
结果信息，包含各种统计信息的excel文件(results_information.xlsx)和日志(log.txt)，后文将会对统计信息进行详细说明
第三部分
结果汇总文件，根据基因类型可分为四种。nuclear_gene, cp_gene, mito_gene, target_gene分别对应nuclear_genes,cp_genes,mito_genes,target_genes.  
结果汇总文件还可以细分为：
a.target/cp/mito/nuclear_genes_fasta
    储存参考基因组的文件夹。根据参考序列，将每个基因单独写为一个fasta文件，这些文件会作为后续过滤reads的参考序列
b.filtered_out
	储存经filter脚本过滤后的的reads，格式为fastq格式
c.assembled_out
    储存经minia拼接后的contigs和untigs
d.GM_results:
	将所有挖掘的基因汇总在一起。一类是原始结果（gene.fasta）,一类是经过各种条件筛选后根据参考序列比对切齐的结果 （gene.trimmed.fasta）.值得注意的是，当两种结果都存在的时候，结果往往更准确。当然，用户也可以根据自己的需要将序列切齐
e.bootstrap_out(可选)
	存储自展检测的结果，包括参考序列库（bootstrap_db），过滤结果（filter_out）,组装结果（assembled_out）,最终挖掘结果（GM_results）和自展结果信息（GM_gene_trimmed_bootstrap.xlsx）。
	值得注意的是，当GM_results中不存在gene.trimmed.fasta，将不会参与自展检测


example:
GeneMiner -1 forward.fq  -2 reverse.fq  -rt ref --check True --bootstrap_number 10 -o GeneMiner_out
tree -L 4 GeneMiner_result/
GeneMiner_result/
|-- forward.fq
|-- log.txt
|-- results_information.xlsx
|-- reverse.fq
`-- target_genes
    |-- GM_results
    |   |-- GM_matK.fasta
    |   |-- GM_matK_trimmed.fasta
    |   |-- GM_psbA.fasta
    |   `-- GM_psbA_trimmed.fasta
    |-- assembled_out
    |   |-- target_matK
    |   |   |-- assembled_out.contigs.fa
    |   |   `-- assembled_out.unitigs.fa
    |   `-- target_psbA
    |       |-- assembled_out.contigs.fa
    |       `-- assembled_out.unitigs.fa
    |-- bootstrap_out
    |   |-- target_GM_matK_trimmed
    |   |   |-- GM_matK_trimmed_bootstrap.xlsx
    |   |   |-- GM_results
    |   |   |-- assembled_out
    |   |   |-- bootstrap_db
    |   |   `-- filtered_out
    |   `-- target_GM_psbA_trimmed
    |       |-- GM_psbA_trimmed_bootstrap.xlsx
    |       |-- GM_results
    |       |-- assembled_out
    |       |-- bootstrap_db
    |       `-- filtered_out
    |-- filtered_out
    |   |-- target_matK
    |   |   |-- Filtered_reads__R1.fastq
    |   |   |-- Filtered_reads__R2.fastq
    |   |   `-- filtered.fq
    |   `-- target_psbA
    |       |-- Filtered_reads__R1.fastq
    |       |-- Filtered_reads__R2.fastq
    |       `-- filtered.fq
    `-- target_genes_fasta
        |-- matK.fasta
        `-- psbA.fasta


```



### 5.2.2 结果统计信息说明

**results_information.xlsx详细说明**

| 列名                      | 内容       | 注释                                                         |
| ------------------------- | ---------- | ------------------------------------------------------------ |
| nuclear_gene_name         | ITS        | 基因名                                                       |
| filtered_reads_number     | 300        | 经filter过滤后获得的reads                                    |
| richness                  | 75.88      | 丰度，过滤后的reads比对到参考序列上的平均深度，richness=(n*L1)/L2  <br /> n: reads数量 ，L1 reads平均长度  L2: 基因长度 |
| graph_construction        | 30.766     | graph construction step(minia)                               |
| assembled_percentage      | 0.186      | 可以用于拼接的reads百分比                                    |
| assembled_max_length      | 1181       | 拼接后获得的最长序列的长度                                   |
| result_max_length         | 1181       | 经筛选后的结果中最长序列的长度                               |
| identity_trimmed_sequence | 85.83%     | 经剪切后的序列与参考序列的一致度                             |
| coverage_trimmed_sequence | 98.52%     | 参考序列对于经剪切序列的覆盖度                               |
| Gene_extraction           | successful | 基因是否提取成功                                             |
| Gene_aligned_cut          | successful | 基因是否能够对齐剪切成功                                     |



### 5.3 例子

#### 5.3.1 提取叶绿体基因：

```shell
GeneMiner -1 /data/Daucus_carota.1.fq.gz -2 /data/Daucus_carota.2.fq.gz  -rmito ref_cp.gb -b 0 -max  3000    -min 300 -o Daucus_carota_cp_genes 
```

当有较为充足的数据量，合适的参考序列时，本软件几乎能从浅层基因组数据中提取所有的叶绿体基因。这为叶绿体组装不成环提供了第二种解决方法。

#### 5.3.2 提取线粒体基因：

```shell
GeneMiner -1 /data/Daucus_carota.1.fq.gz -2 /data/Daucus_carota.2.fq.gz -rmito ref_mito.gb -o Daucus_carota_mito_genes
```

我们挖掘到的线粒体主要在（mito_genes/GeneMiner_results)文件夹内，其中有后缀为trimmed.fasta的文件，说明该结果较好，而其他没有trimmed后缀的文件可能是因为序列较短而不能进一步的剪切。

对于植物来说，线粒体基因所发生的突变较多，因此，如果你想挖掘出更多线粒体基因，可能需要输入更多的数据，或寻找更加近缘的参考序列。

由于本软件所取数据默认为原始数据的前1000000行reads，因此，你可以在现有数据的基础上适当提高-n的输入：

```shell
#使用前3000000行的数据输入
GeneMiner -1 /data/Daucus_carota.1.fq.gz -2 /data/Daucus_carota.2.fq.gz -rmito ref_mito.gb -o Daucus_carota_mito_genes -n 3000000
```

#### 5.3.3 提取核基因：

```shell
#ITS
GeneMiner -1 /data/Daucus_carota.1.fq.gz -2 /data/Daucus_carota.2.fq.gz -rn ITS_ref.fasta -o Daucus_carota_ITS
#18S
GeneMiner -1 /data/Daucus_carota.1.fq.gz -2 /data/Daucus_carota.2.fq.gz -rn 18S_ref.fasta -o Daucus_carota_18S
#ETS
GeneMiner -1 /data/Daucus_carota.1.fq.gz -2 /data/Daucus_carota.2.fq.gz -rn ETS_ref.fasta -o Daucus_carota_ETS
#26S
GeneMiner -1 /data/Daucus_carota.1.fq.gz -2 /data/Daucus_carota.2.fq.gz -rn 26S_ref.fasta -o Daucus_carota_26S
```

#### 5.3.4 提取单拷贝基因

我们推荐从转录组数据中提取核单拷贝基因，可以根据目前已被开发的适用于系统发育构建的被子植物353单拷贝基因作为参考序列。

```shell
GeneMiner -1 /data/Daucus_carota.1.fq.gz -2 /data/Daucus_carota.2.fq.gz -rt 353_fastafile -o Daucus_carota_35

```

值得注意的是

第一：由于原始数据测序深度较浅，而单拷贝基因仅在单倍体基因组中出现一次，所有有一定概率并不能得到结果。

第二：挖掘某个基因特别是单拷贝基因时，若提供多个参考序列且参考序列之间差异度太大,GeneMiner将产生多个结果. 因此，输入的参考序列建议使用单条序列，或者输入的参考序列之间的差异度应当保证在15%之内。

第三：由于输入的原始数据是转录组，参考序列应该为对应的mRNA等转录本信息



## 6.结果校验(可选)

采用bootstrap方法，对GeneMiner产生的结果进行检测。

```shell
#自展检验10次
GeneMiner -1 /data/Daucus_carota.1.fq.gz -2 /data/Daucus_carota.2.fq.gz -rn ITS_ref.fasta -o Daucus_carota_ITS --check True --bootstrap_number 10
```



## 7.原理详细介绍

等小郭的图





## 8. 常见问题

### 8.1 如何验证结果的可靠性？

​		我们的软件设置了可以以自展检验的方式对结果进行检验，如果你对你的结果产生了质疑，可以使用--check选项，并设置自展次数，但自展检验需要反复的调用GeneMiner，因此，建议设置一个合适的值，以免耗费过多时间。

​		同时也可以利用NCBI的二代测序数据或者真实的一代数据作为补充验证



### 8.2 没有得到结果怎么办？

​		影响基因挖掘成功率的的主要原因大概有如下方面：第一，原始数据质量的好坏。相较于raw_data，clean_data的效果更好。第二，数据量的大小，数据量太小可能导致基因丰度不够无法拼接；数据量太大可能导致多种拼接情况的出现。我们推荐将-n设置为100w~1000w。第三，参考序列的选择。选择近缘属或者同属不同种的序列作为参考序列，往往能够得到较好的结果。同时，多条参考序列的时候，注意参考序列之间的差异度不要太大。



### 8.3 有GM_gene.fasta而没有GM_gene_trimmed.fasta怎么办

​		GM_gene_fasta是GeneMiner挖掘出原始contigs，后续对挖掘出来的contigs做了进一步的过滤。比如根据参考序列的长度，删除了太长或太短的序列；根据参考序列的方向，判断contigs是否需要反向互补从而保证contigs和参考序列方向一致；使用双序列局部比对方法确定出参考序列和contigs中高度吻合的部分。经过上述等方法最终获得GM_gene_trimmed.fasta。如果仅有GM_gene_trimmed.fasta而没得到GM_gene.fasta，说明原始contigs不满足所有过滤要求。值得注意的是，GM_gene.fasta依然值得使用，用户可以根据自己的需求对序列进行比对切齐。



### 8.4 结果有多条序列怎么办？

​		出现这种情况的原因可能是由于用户提供了多条参考序列，且参考序列之间的变异度大于20%，对于这种情况，我们推荐使用单条参考序列。一般来说，多条序列的结果中应当有一条序列是正确的，用户可以尝试用blast进行手动筛选并检验。在未来的版本中，我们将对用户输入的多条参考序列进行处理，去除相对差异较大的序列。



### 8.5 如何寻找合适的参考序列？



（到时候会补充一个随着变异度变化  结果准确度变化的图）

​	

### 8.6  -n 设置多少合适？

​	提取较为保守的序列，用默认的1000000行的数据，即250000条reads就足够了，例如18s、26s、5.8s，如果对这些序列用过多的reads反而可能提取不到结果，但对于ITS、ETS这样在nrDNA中相对不太保守的序列来说，则需要用更多的数据。



（到时候会补充一个改变数据量 结果准确度变化的图  ）



## 9.引用

当你使用GeneMiner的时候请引用：

GeneMiner :  a software for extracting phylogenetic markers from skimming genome



请同时引用：

K. Salikhov, G. Sacomoto and G. Kucherov. [*Using cascading Bloom filters to improve the memory usage for de Brujin graphs*](http://minia.genouest.org/files/cascading-wabi13.pdf), WABI 2013

Chikhi R, Rizk G. Space-efficient and exact de Bruijn graph representation based on a Bloom filter[J]. Algorithms for Molecular Biology, 2013, 8(1): 1-9.

