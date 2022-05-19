

# GeneMiner

​	感谢您选择GeneMiner(GM)。本文档将帮助您学会如何使用GeneMiner。请注意，GeneMiner仍处于测试阶段，因此本文档可能会在未来的版本升级中更新。

​	如果您有任何问题，可以联系**xiepulin@stu.scu.edu.cn**  或者 15162176893@163.com

[TOC]

# 1. 简介

​		GeneMiner(Gene-Miner，基因矿工)是一款用于从二代测序(NGS)数据中挖掘基因的软件，能够从极低质量和深度的源数据中获取高质量的特定目标基因。例如从浅层基因组测序数据中准确提取叶绿体/线粒体的全部或者部分基因、核基因组中高度重复区 (如*nr*DNA)等；从转录组测序数据中提取大量的单拷贝系统发育标记；从宏基因组数据中获取特定微生物的环境响应基因等。可广泛应用于系统发育与进化研究、海关检验检疫、特定功能基因的挖掘等，在降低测序成本，扩大基因选择方面具有显著优势。GeneMiner的结果十分准确，在真实的实验验证中达到或接近了一代测序的结果，即便参考序列与目标序列相似度低于90%，依然能够依靠梯度逼近算法获得100%准确的结果。软件创新性的提出了基于自展检测的校验方法，可以在不依赖参考序列的情况下对目标序列进行重复验证，输出更加可靠的一致性序列。基于算法层面的大量优化，GeneMiner的运算速度和内存消耗都非常优越，支持多线程并行，充分调用计算机资源，既能够部署在高性能运算集群上，也可以在普通个人电脑上运行。GeneMiner对用户非常友好，支持Windows、Mac和Linux各种主流的操作系统，用户可以选择命令行界面或图形界面进行使用。

# 2.下载

GeneMiner 在MIT License下是开源的。它通过github存储库分发:https://github.com/sculab/GeneMiner，你可以随时下载最新的的版本。请务必关注github，以保持最新的代码更改。我们对以前版本的代码不提供任何支持!版本号遵循符号x.y.z，其中x随着主要的代码重组而改变，y当添加新特性时发生更改，并且伴随着bug修复而发生更改

# 3. 安装

## 3.1 For Linux users
自动安装（推荐）
```shell
tar -zxvf geneminer.tar.gz # 解压缩
cd geneminer
python install.py   #根据脚本提示完成，自动安装依赖并将GeneMiner写入环境变量
```
手动安装（自动安装遇到问题时使用）
```shell
tar -zxvf geneminer.tar.gz # 解压缩
cd geneminer
# 手动安装所需库
pip3 install  biopython --user
#也可以根据软件提供的依赖文件安装
pip3 install -r requirements.txt --user

# 配置环境变量
echo "export PATH=\$PATH:$(pwd)" >> ~/.bashrc
source ~/.bashrc
```
关闭并重启终端，检测是否配置成功
```shell
geneminer.py -h # 命令行界面
```

## 3.2 For Mac users
下载对应版本的打包好的图形界面app，双击即可运行。（推荐）
要手动配置命令行和图形界面版本，使用如下命令进行自动安装：

```shell
tar -zxvf geneminer.tar.gz # 解压缩
cd geneminer 
python3 setup.py #脚本会自动安装依赖，并将GeneMiner写入环境变量
```
或者使用如下命令手动安装：
```shell
tar -zxvf geneminer.tar.gz # 解压缩
cd geneminer
# 手动安装所需库
pip3 install  biopython --user
#也可以根据软件提供的依赖文件安装
pip3 install -r requirements.txt --user

# 配置环境变量
# 对于macOS Catalina (10.15) 及其之后的系统：
echo "export PATH=\$PATH:$(pwd)" >> ~/.zshrc
source ~/.zshrc
# 对于macOS Catalina (10.15) 之前的系统：
echo "export PATH=\$PATH:$(pwd)" >> ~/.bash_profile
source ~/.bash_profile
```
关闭并重启终端，检测是否配置成功
```shell
geneminer.py -h # 命令行界面
```
## 3.3 For Windows users
### 3.3.1 安装打包好的图形界面应用程序（推荐）
下载GeneMiner对应版本的图形界面应用程序，双击即可运行。GeneMiner windows版maual详见：xxxxxxxxxxxxxxxxxx
### 3.3.2 配置命令行界面应用程序
配置命令行界面的GeneMiner需要您在系统中安装python 3.6以上的版本，您也可以安装Anaconda或者miniconda，具体安装方式请参阅python及其不同发行版的官方技术文档。您可以参考如下步骤配置GeneMinier：
自动安装

- 下载GeneMiner的Windows安装包，双击进行自动安装。
  如果自动安装遇到问题，可以执行手动安装：

手动安装

- 解压缩：将下载的windows安装包解压到指定的文件夹，例如geneminer。
- 打开命令行，手动安装所需库：
```shell
pip3 install  biopython --user
pip3 install  pysimplegui --user
#也可以根据软件提供的依赖文件批量安装
pip3 install -r requirements.txt --user
```
- 将geneminer文件夹加入用户环境变量path中。
- 关闭并重启终端，检测是否配置成功
```shell
geneminer.py -h # 命令行界面
```

# 4.快速入门

​		在对您的测序数据进行挖掘之前，我们建议您先了解自己的测序数据，包括测序方式、深度、质量、数据量大小等等。我们的软件主要适用于Illumina平台所返回的二代测序数据。经大量的测试数据验证，即便对于较低的测序深度(10x以下)，GeneMiner也能从中挖掘到单拷贝核基因(nuclear gene)、叶绿体基因(cp gene)以及线粒体基因(mito gene)等，也可以用于从转录组数据中挖掘单拷贝基因。

​		用于练手的数据存放在 GeneMiner/example/

（1）挖掘单个目标基因

```shell
 geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa matK.fasta 
```
其中data1.fq.gz和data2.fq.gz是二代测序返回的测序数据，matK.fasta是近缘物种（同属或者同科）的同源基因。
**注意：每一个fasta文件中只能存放同一种基因（可以是不同物种的）**。如以下样例：

```
>GeneA_species1
ATCGATCG
>GeneA_species2
ATCGATCC
>GeneA_species3
ATTGATCC
```

（2）挖掘多个目标基因

可以将多个fasta格式的文件放在同一个文件夹中，便能够批量挖掘基因

```shell
 geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref_fasta
```

其中cp.gb是GenBank格式的叶绿体参考基因组，mito.gb是GenBank格式的线粒体参考基因组，-rcp和-rmito用于指定参考基因组类型


（3）评估可靠性

```shell
geneminer.py -1 data1.fq.gz -2 data2.fq.gz  -rtfa  matK.fasta  -bn 10
```

-bn 自展检测的次数。采用基于Bootstrap思想的方法评估结果的准确性。在不依赖参考序列的情况下对目标序列进行重复验证，输出更加可靠的一致性序列。



# 5.详细使用指南

## 5.1参数解读

使用geneminer -h ，计算机将显示所有选项，后面我将详细介绍它们的使用方法和一些小窍门

```shell
GeneMiner: a software for extracting phylogenetic markers from next generation sequencing data
Version: 1.0.0
Copyright (C) 2022 Pulin Xie
Please contact <xiepulin@stu.edu.scu.cn> if you have any questions

optional arguments:
  -h, --help          show this help message and exit

Basic option:
  -1                  One end of the paired-end reads, support fastq format
  -2                  Another end of the paired-end reads, support fastq format
  -s , --single       Single reads, support fastq format
  -o , --out          Specify the result folder
  -rtfa <file|dir>    References of target genes, only support fasta format
  -rtgb <file|dir>    References of target genes, only support GenBank format

Advanced option:
  -k1 , --kmer1       Specify the size of the wordsize to filter reads  [default = 29]
  -k2 , --kmer2       Specify the size of the kmer to assemble reads  [default = 31]
  -d , --data_size    Specifies the number of reads to reduce raw data. If you want to                         use all the data, you can set as 'all' [default = 'all']
  -step_length        the length of the sliding window on the reads [default = 4]
  -limit_count        limit of kmer count [default=auto]
  -limit_min_length   limit of contig length
  -limit_max_length   limit of contig length
  -change_seed        times of changing seed [default = 32]
  -scaffold           make scaffold
  -max                The maximum length of contigs [default = 5000]
  -min                The minimum length of contigs [default = 300]
  -t , --thread       The number of threads [default = 'auto']
  -b , --boundary     Extend the length to both sides of the gene while extracting genes 					   from Genbank file [default = 75]
  -bn , --bootstrap   Specify the bootstrap number. Evaluate the results based on the 			              bootstrap method

```



### 5.1.1 基础参数

```shell
-1				
双末端测序数据其中一端的数据，支持fastq/fastq.gz格式。务必保留正确的文件扩展名。
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa  ref.fasta

-2   			
双末端测序数据中另一端的数据，支持fastq/fastq.gz/fastq.bz2格式。务必保留正确的文件扩展名。
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta

-s   			
单端测序数据，支持fastq/fastq.gz格式。务必保留正确的文件扩展名。
example: geneminer.py -s data1.fq.gz -rtfa  ref.fasta

-o , --out            
指定输出文件夹。
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -o geneminer_out

-rtfa     <file|dir>  
目标基因参考序列，仅支持fasta格式。可用于寻找自己感兴趣的基因。
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta
同一个fasta格式的文件里，只能存放同一个基因的数据，例子如下：
>species_a ITS
AGCTAGCT
>species_b ITS 
AGCTAGCC
>species_c ITS
AGCTAGCA
>species_d ITS
AGCTAGAA

-rtgb     <file|dir>  		
目标基因参考序列，仅支持GenBank格式。可用于寻找自己感兴趣的基因。
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtgb gb_folder
```



### **5.1.2 高级参数**

```shell
-k1 , --kmer1  
指定wordsize的大小。geneminer将长度为L的read拆分为(L-k1+1)条长度为k1的word,如果其中某一个word能与参考序列匹配，则将这条read保留。k1的大小。取决于亲缘关系，亲缘关系越近，待挖掘基因与参考序列相似度越大，k1取值越大。该参数用于过滤reads,默认大小为17，kmer的取值范围在17~127之间

-k2 , --kmer2   
指定k-mer长度，k-mer是de Bruijn图中节点的长度。kmer的取值严重依赖于数据集。
该参数用于拼接reads，默认大小为31，kmer的取值范围在17~127之间
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -n 2000000 -k1 17 -k2 31

-d , --data     
指定reads的数目。对于测序深度较高的仅选择原始数据的一部分（100w~1000w）就已经能得到很好的结果，同时又大大减少了软件运行时间
对于读长为150bp的二代测序数据，1000w行的数据量大概在800~1000MB之间。

如果您需要将全部的原始数据输入，可以将-n 设置为all

如果您对原始数据量的行数感兴趣，您可以使用以下命令查看您的数据
zcat your_data.fq.gz|wc -l 
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -n 2000000
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -n all

-step_length
指定reads过滤滑动窗口的长度 默认值为4
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -step_length 4

-limit_count
kmer数量的最低阈值
在reads过滤的过程中，我们会将reads拆分为长度为k1的子序列，并统计这些kmer出现的次数
我们自然的认为，低频的kmer置信度较低。所以，geneminer将移除kmer数量（kmercount）低于limit_count的kmer
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -limit_count auto

-limit_min_length
目标序列占参考序列平均长度的最小比率
limit_min_length=目标序列长度/参考序列平均长度 默认值=1

-limit_max_length
目标序列占参考序列平均长度的最大比率
limit_max_length=目标序列长度/参考序列平均长度 默认值=2
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -limit_min_length 0.5 -limit_max_length 1.5


-max 
指定挖掘基因的最大长度，默认为5000bp.
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rcp cp.gb -n 2000000 -k 43 -max 3000

-min    
指定挖掘基因的最小长度，默认为300bp.
example:
geneminer.py  -1 data1.fq.gz  -2 data2.fq.gz -rcp mito.gb -n 2000000 -k 43 -max 3000     -min 200

-t , --thread  
example:
指定线程数量 ，如果不指定的话，软件会根据计算机性能自动选择合适的线程数量
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rcp cp.gb -n 200000 -k 43 -max 3000 -min 200 -t 8

-b , --boundary   
指定软边界的长度。
当沿着挖掘出的基因的两侧延申的时候，随着延申长度的增加，碱基的准确率越来越低。但这种下降不是断崖式的，而是在某一定长度的缓冲区内逐渐下滑。我们将保留一段缓冲区，并将这段缓冲区称为软边界。推荐大小为0.5*reads的长度，软边界取值范围在0~200之间。该参数与-rcp参数同时使用时才生效
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rcp cp.gb -n 2000000 -k 43 -max 3000 -min 500 -b 75

-bn , --bootstrap
指定自展检测的次数
基于自展检测的校验方法，可以在不依赖参考序列的情况下，对结果进行准确性评估，并对目标序列进行重复验证
[默认值=10]
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta  -bn 100
```





## 5.2 结果解读

```shell

```

## 5.3 例子

### 5.3.1 提取叶绿体基因：

​		当有较为充足的数据量，近缘的参考序列时，本软件几乎能从浅层基因组数据中提取所有的叶绿体基因同时在系统发育研究中.这为叶绿体组装不成环提供了另一种解题思路。

```shell
#example
geneminer.py -1 data1.fq.gz -2 data2.fq.gz  -rmito ref_cp.gb -b 0 -max  3000    -min 300 -o results
```

### 5.3.2 提取线粒体基因：

​		线粒体基因挖掘难度往往会比叶绿体基因大很多，其一般原因在于：原始数据本身未包含太多线粒体基因，线粒体基因变异较大。针对这种情况，您可以适当增加原始数据量大小以及选择更近源的参考序列.

```shell
#example
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rmito mito.gb  -n 15000000
```

### 5.3.3 提取核基因：

​        在测序深度较低的测序数据中，GeneMiner也能挖掘出中高拷贝数基因，如rDNAs

```shell
#example
#ITS
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rn ITS_ref.fasta -o ITS_out
#18S
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rn 18S_ref.fasta -o 18S_out
#ETS
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rn ETS_ref.fasta -o ETS_out
#26S
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rn 26S_ref.fasta -o 26S_out
```

# 6.方法

GeneMiner 核心流程主要分为三个步骤：

![](https://gitee.com/xiepulin/picgo_xpl/raw/master/GeneMiner_picture/流程图改改改.svg)



## **6.1数据过滤：**

​		区别传统的先拼接序列再映射到参考基因组的方法，GeneMiner选择了先将序列映射到参考基因组再组装序列。首先将格式为fastq的二代数据作为原始的输入数据，根据实际情况，选择大小合适的数据量。接着采用Ukkonen’s算法将用户提供的参考基因组构建为一棵后缀树。随后，将原始数据中长度为L的reads拆分为（L-k+1）个长度为k的子序列，其中L的大小取决于测序方式，k的大小取决于用户的设定的kmer。最后，如果检测到某一个reads的kmer在后缀树中出现，则把该条reads保留。

​		由于采用的先过滤后拼接的策略，需要过滤的reads数据量相当大，需要拼接的reads则大幅度减少，这意味着拼接不再是计算瓶颈，过滤这一步骤才是本方法需要突破的核心。幸运的是，我们将过滤这一步骤的算法复杂读降低到O（n），大幅度压缩了时间成本。



## 6.2 拼接与校正：

### 6.2.1数据拼接

​       将过滤后的reads拼接为contigs，我们选择minia3作序列拼接软件。主要原因在于，区别于其他的一些序列拼接软件如velvet、spades、soapdenovo,minia3准确度高，速度快，消耗的内存更小，更加适用于短序列的拼接。

(Ⅰ)make kmers

把过滤出的reads拆分为更小的片段k-mers

(Ⅱ)remove low quality kmers

绘制kmer频数分布曲线，横坐标为kmer的频数，纵坐标为该频数的kmer的总数，对曲线做平滑处理，第一次拐点的位置的横坐标即为作为limit的阈值。凡是kmer频次低于阈值的kmer,将作为低质量的kmer剔除。由于每一个目标loci的kmer频数分布曲线不尽相同，所以，每个目标基因都有其自适应的limit

(Ⅲ)choose seed

(Ⅳ)build graph

(Ⅴ)walk graph and output contig

### 6.2.2数据校验

​      minia3可能产生多条contigs,而这些序列不都是我们需要的，因此GeneMiner对这些contigs做了进一步的处理。

**（1）长度校正**

​		根据参考的序列的平均长度，过滤掉长度太长或者太短的contigs。

**（2）方向校正**

​		GeneMiner利用BLAST+套件中的makeblastdb和balstn工具，将与参考序列的方向一致的contigs保留，与参考序列方向不一致的contigs做反向互补处理。



## 6.3自展检测

​		GeneMiner创新性的提出了基于自展检测的校验方法，可以在不依赖参考序列的情况下，对结果进行准确性评估，并对目标序列进行重复验证，输出更加可靠的一致性序列。

(Ⅰ)首次获取目的基因后，与参考序列进行比对，得到变异率v。为了防止结果中有多条contigs或有多条参考序列导致差异度不方便计算，GeneMiner根据序列之间的一致度与覆盖度，选择了契合度最高的contig与参考序列作为计算差异度的标准。

(Ⅱ)将目的基因每个位点以变异率v进行随机重采样，获取变异率同样为v的ref1, ref2.....refn

(Ⅲ)使用ref1, ref2.....refn作为参考序列重新运行整个过程，得到新的目的基因target1, target2.... targetn

(Ⅳ)对target1, target2.... targetn获取一致性序列，位点的一致性比例即为该位点的支持率		

# 7. 常见问题

## 7.1 如何验证结果的可靠性？

​		我们的软件设置了基于自展检测的方法对结果进行检验，如果您对您的结果产生了质疑，可以使用-bn选项，并设置自展次数，但自展检验需要反复的调用GeneMiner内部脚本，因此，建议您设置一个合适的值，以免耗费过多时间。同时也可以利用NCBI上的二代测序数据或者真实的一代数据作为补充验证

## 7.2 没有得到结果怎么办？

​		影响基因挖掘成功率的的主要原因大概有如下方面：第一，原始数据质量的好坏。相较于raw_data，clean_data的效果更好。第二，数据量的大小，数据量太小可能导致基因丰度不够无法拼接；数据量太大可能导致多种拼接情况的出现。我们推荐将-n设置为100w~1000w。第三，参考序列的选择。选择近缘属或者同属不同种的序列作为参考序列，往往能够得到较好的结果。同时，当存在多条序列作为参考序列的时候，参考序列之间的差异度不要太大。

## 7.3  -n 设置多少合适？

​		一般来说使用默认参数n=10000000，即2500000条reads即可以很好的满足需要。但是对于单拷贝或者低拷贝序列，可以适当的增加数据量。我们不建议将全部数据输入，因为不仅会大大降低速度，而且由于拼接出的情形增多，反而得不到好的结果。

# 8.引用

当你使用GeneMiner的时候请引用：

GeneMiner : a software for extracting phylogenetic markers from next generation sequencing data

