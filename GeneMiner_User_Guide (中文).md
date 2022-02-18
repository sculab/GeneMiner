

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
python setup.py   #根据脚本提示完成，自动安装依赖并将GeneMiner写入环境变量
```
手动安装（自动安装遇到问题时使用）
```shell
tar -zxvf geneminer.tar.gz # 解压缩
cd geneminer
# 手动安装所需库
pip3 install  biopython --user
pip3 install  pandas --user
pip3 install  tqdm --user
pip3 install  openpyxl --user
pip3 install  pysimplegui --user
#也可以根据软件提供的依赖文件批量安装
pip3 install -r requirements.txt --user

# 配置环境变量
echo "export PATH=\$PATH:$(pwd)" >> ~/.bashrc
source ~/.bashrc
```
关闭并重启终端，检测是否配置成功
```shell
geneminer -h # 命令行界面
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
pip3 install  pandas --user
pip3 install  tqdm --user
pip3 install  openpyxl --user
pip3 install  pysimplegui --user
#也可以根据软件提供的依赖文件批量安装
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
geneminer -h # 命令行界面
```
## 3.3 For Windows users
### 3.3.1 为Windows安装wsl
由于GeneMiner需要在Windows中使用wsl，因此只支持Windows 10及其之后的操作系统。wsl并不是系统默认的功能，您首先需要在系统中安装wsl。对于Windows 10 版本 2004 及更高版本（内部版本 19041 及更高版本）或 Windows 11，您可以使用如下步骤安装wsl，安装过程中需要联入互联网：

- 使用管理员身份运行PowerShell：在开始菜单中找到 Windows PowerShell，右键单击选择“以管理员身份运行” 。
- 在打开的命令行窗口中，输入：
```shell
wsl --install
```
- 安装完成后在命令行窗口运行wsl以确认安装成功。
```shell
wsl
```
首次启动新安装的 Linux 发行版时，将打开一个控制台窗口，要求你等待将文件解压缩并存储到计算机上。 未来的所有启动时间应不到一秒。

对于旧版的Windows 10，建议您升级到最新版本，或者使用旧版的手动安装方式，具体参见微软公司的技术文档：
- https://docs.microsoft.com/zh-cn/windows/wsl/install-manual
### 3.3.2 安装打包好的图形界面应用程序（推荐）
下载GeneMiner对应版本的图形界面应用程序，双击即可运行。GeneMiner windows版maual详见：xxxxxxxxxxxxxxxxxx
### 3.3.3 配置命令行界面应用程序
配置命令行界面的GeneMiner需要您在系统中安装python 3.6以上的版本，您也可以安装Anaconda或者miniconda，具体安装方式请参阅python及其不同发行版的官方技术文档。您可以参考如下步骤配置GeneMinier：
自动安装

- 下载GeneMiner的Windows安装包，双击进行自动安装。
  如果自动安装遇到问题，可以执行手动安装：

手动安装

- 解压缩：将下载的windows安装包解压到指定的文件夹，例如geneminer。
- 打开命令行，手动安装所需库：
```shell
pip3 install  biopython --user
pip3 install  pandas --user
pip3 install  tqdm --user
pip3 install  openpyxl --user
pip3 install  pysimplegui --user
#也可以根据软件提供的依赖文件批量安装
pip3 install -r requirements.txt --user
```
- 将geneminer文件夹加入用户环境变量path中。
- 关闭并重启终端，检测是否配置成功
```shell
geneminer -h # 命令行界面
```

# 4.快速入门

​		在对您的测序数据进行挖掘之前，我们建议您先了解自己的测序数据，包括测序方式、深度、质量、数据量大小等等。我们的软件主要适用于Illumina平台所返回的二代测序数据。经大量的测试数据验证，即便对于较低的测序深度(10x以下)，GeneMiner也能从中挖掘到单拷贝核基因(nuclear gene)、叶绿体基因(cp gene)以及线粒体基因(mito gene)等，也可以用于从转录组数据中挖掘单拷贝基因。

​		用于练手的数据存放在 GeneMiner/example/

（1）挖掘单个目标基因

```shell
 geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta 
```
其中data1.fq和data2.fq是二代测序返回的测序数据，ref.fasta是近缘物种（同属或者同科）的同源基因。在本例中，
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

可以将多个fasta格式的文件放在一个文件夹中，便能够批量挖掘基因

```shell
 geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta
```

（3）批量挖掘不同类型的基因

```shell
 geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rcp cp.gb -rmito mito.gb
```

其中cp.gb是GenBank格式的叶绿体参考基因组，mito.gb是GenBank格式的线粒体参考基因组，-rcp和-rmito用于指定参考基因组类型


（4）评估可靠性

```shell
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta  -bn 10
```

-bn 自展检测的次数。采用基于Bootstrap思想的方法对结果准确性进行评估。设置的自展检测数值越大，计算耗时就越大



# 5.详细使用指南

## 5.1参数解读

使用geneminer -h ，计算机将显示所有选项，后面我将详细介绍它们的使用方法和一些小窍门

```shell
geneminer -h
usage: GeneMiner <-1 -2|-s|-12>  <-rn|rcp|-rmito|rt>  [options]

GeneMiner: a software for extracting phylogenetic markers from skimming genome
Version: 1.0
Copyright (C) 2021 Pu-lin Xie
Please contact <xiepulin@stu.edu.scu.cn>, if you have any bugs or questions

optional arguments:
  -h, --help            show this help message and exit

Basic option:
  -1                    one end of paired-end reads,support fastq/fastq.gz/fastq.bz2
  -2                    another end of paired-end reads,support fastq/fastq.gz/fastq.bz2
  -12                   interlaced forward and reverse paired-end reads,support fastq/fastq.gz/fastq.bz2
  -s , --single         single-read,support fastq/fastq.gz/fastq.bz2
  -o , --out            Specify the result folder [default='auto']
  -rcp <file|dir>       reference of chloroplast genome,only support GenBank-format
  -rmito <file|dir>     reference of mitochondrial genome,only support GenBank-format
  -rtfa <file|dir>      References of target genes, only support fasta format
  -rtgb <file|dir>      References of target genes, only support GenBank format

Advanced option:
  -n , --number         The number of rows of raw data from skimming genomes,default=1000000
  -k , --kmer           size of a kmer  [default =31]
  -max                  gene maximum length  [default =5000]
  -min                  gene minimum length  [default =300]
  -t , --thread         Specify the number of threads you want to run [default='auto']
  -b , --boundary       extend the length to both sides of the gene while extracting                             					 genes from  Genebank file [default=75]
  -sf                   Select the reference sequences to reduce the computation.
                        s1: do nothing;
                        s2,3,4: only use the reference sequence with the shortest/median/longest length;
                        s5: remove sequences with abnormal length.[default = 's1']

Gradient approximation option:
  -in , --iterative_number
                        Specify the number of iterative loop to gradually approximate the best results

Bootstrap option:
  -bn , --bootstrap_number 
                        Specify the bootstrap number.Evaluate the results based on the bootstrap method
```



### 5.1.1 基础参数

```shell
-1				
双末端测序数据其中一端的数据，支持fastq/fastq.gz/fastq.bz2格式。务必保留正确的文件扩展名。
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa  ref.fasta
    
-2   			
双末端测序数据中另一端的数据，支持fastq/fastq.gz/fastq.bz2格式。务必保留正确的文件扩展名。
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta

-12  			
交错合并的双末端测序数据，支持fastq/fastq.gz/fastq.bz2格式。务必保留正确的文件扩展名。
example: geneminer.py -12 data.fq.gz -rtfa  ref.fasta

-s   			
单端测序数据，支持fastq/fastq.gz/fastq.bz2格式。务必保留正确的文件扩展名。
example: geneminer.py -s data1.fq.gz -rtfa  ref.fasta

-rcp 	<file|dir>	
叶绿体参考基因组，只支持GenBank的格式。
可以输入一个GenBank文件，该文件里既可以包含一个种，也可以包含多个种；还可以将多个Genebank文件放入一个文件夹下。 
值得注意的是，用户可以修改GeneBank文件，仅保留自己感兴趣的部分基因。同时在选择参考基因组的时候尽量选择近缘参考基因组
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rcp cp.gb

-rmito    <file|dir>	
线粒体参考基因组，只支持GenBank的格式。具体用法同rcp
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rmito mito.gb

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
t4 species_d ITS
AGCTAGAA

-rtgb     <file|dir>  		
目标基因参考序列，仅支持GenBank格式。可用于寻找自己感兴趣的基因。
-rtgb的功能囊括了-rcp -rmito。之所以将rcp,rmito独立出来，是为了方便用户使用以及软件后期功能的扩展
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtgb gb_folder

-o , --out            
指定输出文件夹。如果不指定输出文件夹，geneminer.py将自动使用 ''GM+时间戳'' 作为输出文件夹名字.但是windows和mac的图形界面的版本必须指定输出文件夹。
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -o geneminer_out

-sf   
对参考序列进行筛选。当参考序列数目太多或者相互之间差异度太大的时候，不仅增加计算用时同时影响结果的准确性。GeneMiner目前提供了五种策略用于筛选参考序列。默认选择s1
strategy 1(s1):不做任何处理  
strategy 2(s2):仅保留长度最短的一条参考序列
strategy 3(s3):仅保留中位数长度的一条参考序列
strategy 4(s4):仅保留长度最长的一条参考序列
strategy 5(s5):剔除异常长度的序列，保留正常长度的参考序列
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -sf s5
```



### **5.1.2 高级参数**

```shell
-n , --number     
输入原始数据量的行数，默认值为1000万行。如果您需要将全部的原始数据输入，可以将-n 设置为all
经过测试，仅选择原始数据的一部分（100w~1000w）就已经能得到很好的结果，同时又大大减少了软件运行时间
对于读长为150bp的二代测序数据，1000w行的数据量大概在800~1000MB之间。
如果您对原始数据量的行数感兴趣，您可以使用以下命令查看您的数据
zcat your_data.fq.gz|wc -l 
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -n 2000000
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -n all

-k , --kmer   
指定k-mer长度，k-mer是de Bruijn图中节点的长度。kmer的取值严重依赖于数据集。
默认大小为31，kmer的取值范围在15~127之间
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -n 2000000 -k 43

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
当沿着挖掘出的基因的两侧延申的时候，随着延申长度的增加，碱基的准确率越来越低。但这种下降不是断崖式的，而是在某一定长度的缓冲区内逐渐下滑。我们将保留一段缓冲区，并将这段缓冲区称为软边界。推荐大小为0.5*reads的长度，软边界取值范围在0~200之间。该参数可以和-rcp,-rmito,-rtgb配合使用
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rcp cp.gb -n 2000000 -k 43 -max 3000 -min 500 -b 75
```



### 5.1.3 梯度逼近参数（可选）

```shell
-in , ----iterative_number
指定迭代次数，随着迭代次数的增加，结果逐渐趋近最优答案.
其基本原理为：将GeneMiner输出的结果作为下一次输入的参考序列，如果下一次输出的结果和上一次输出的结果高度一致，则停止。否则重复该过程。其中，重复的最大次数为用户指定的次数。[默认值=2]
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta  -cn 2
```



### 5.1.4自展检测参数（可选）

```shell
-bn,--bootstrap_number  
基于自展检测的校验方法，可以在不依赖参考序列的情况下，对结果进行准确性评估，并对目标序列进行重复验证
[默认值=10]
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta  -bn 20
```



## 5.2 结果解读

```shell
GeneMiner将产生大量的文件，主要分为三大部分:
#第一部分
过滤后的测序数据："data1.fq" 和 "data2.fq"

#第二部分
结果信息，包含各种统计信息的excel文件("results_information.xlsx")和日志("log.txt")

"results_information.xlsx":

| Column name			    | Example    |  Explanation                                                    
| ------------------------  | ---------- | ----------------------------------------------------------- 	
| nuclear_gene_name         | ITS        | 基因名                                                        
| filtered_reads_number     | 300        | 过滤后reads的数量                                                
| richness                  | 75.88      | 丰度，richness=(n*L1)/L2                                                                     				  n:reads数量,L1:reads平均长度,L2: 基因长度                   
| graph_construction        | 30.766     | graph construction step(minia)                               
| assembled_percentage      | 0.186      | 用于组装的reads的百分比                                         
| assembled_max_length      | 1181       | 组装后序列最大长度                                       
| result_max_length         | 1181       | 剪切对齐后序列最大长度                             
| identity_trimmed_sequence | 85.83%     | 参考序列与剪切对齐后的序列与的一致度                             
| coverage_trimmed_sequence | 98.52%     | 参考序列对于经剪切对齐后的序列的覆盖度                               
| Gene_extraction           | successful | 基因是否提取成功                                          
| Gene_aligned_cut          | successful | 基因是否对齐剪切成功                                   


#第三部分
结果汇总文件，不同类型的参考序列将生成不同名称的结果汇总文件。-rtfa,-rtgb,-rcp,-rmito分别对应
tfa_genes,tgb_genes,cp_genes, mito_gene  

结果汇总文件还可以细分为:
(1) Reference sequence folder  ("reference_database")
	GeneMiner对用户提供的参考序列做预处理后，将每个基因单独写为一个fasta格式的文件。
(2过滤文件夹 ("filtered_out")
	储存经filter脚本过滤后的的reads，格式为fastq
(3)组装文件夹 ("assembled_out")
	储存经软件minia组装后生成的contigs和untigs
(4)结果文件夹 ("GM_results")
	将所有被挖掘的基因汇总在一起。
	首先，如果未能挖掘到目标基因，将不产生任何结果文件；如果挖掘到目标基因，则将其保留为原始结果("xxx_raw.fasta")，此时未经过任何处理。xxx代表某个基因，后文表述同理。
	然后，如果挖掘到目标基因但该基因不满足后续各种筛选条件，GeneMiner会根据参考序列，从原始结果中挑选一条最佳结果，存储为最佳原始结果("xxx_raw_best.fasta");如果挖掘到目标基因且该基因满足后续各种筛选条件，a)当结果只有一条序列，则生成唯一结果("xxx.fasta")。  b)当结果有多条序列，则生成候选结果("xxx_options.fasta")。
	最后，GeneMiner将唯一结果(xxx.fasta)或候选结果("xxx_options.fasta")中的序列与参考序列比对切齐，并从中挑选一条最佳的序列作为剪切对齐结果("xxx_trimmed.fasta")。值得注意的是，不是所有结果序列都能比对切齐，用户也可以根据自己的具体需要自行比对切齐。

"xxx_raw.fasta"		:原始结果
"xxx_raw_best.fasta":最佳原始结果
"xxx.fasta"		    :唯一结果(推荐)
"xxx_options.fasta"	:候选结果
"xxx_trimmed.fasta"	:剪切对齐结果(推荐)

(5)自展检测文件夹 ("bootstrap_out",可选)
	存储基于自展的方法生成的结果，具体包括参考序列库（“reference_database”），过滤结果(“filtered_out”),组装结果(“assembled_out”),最终挖掘结果(“GM_results”)
	值得注意的是，只有比对切齐的序列才能够采用基于自展的方法进行结果校验。
(6)梯度逼近文件夹("interative_out",可选)
   存贮每一次迭代生成的结果

example:
geneminer.py -1 data1.fq  -2 data2.fq  -rtfa ref.fasta -bn 5 -cn 2 -o GeneMiner_out
tree -L 4 GeneMiner_result/   #查看结果文件目录结构
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

![](picture/流程图改改改.svg)



## **6.1数据过滤：**

​		区别传统的先拼接序列再映射到参考基因组的方法，GeneMiner选择了先将序列映射到参考基因组再组装序列。首先将格式为fastq的二代数据作为原始的输入数据，根据实际情况，选择大小合适的数据量。接着采用Ukkonen’s算法将用户提供的参考基因组构建为一棵后缀树。随后，将原始数据中长度为L的reads拆分为（L-k+1）个长度为k的子序列，其中L的大小取决于测序方式，k的大小取决于用户的设定的kmer。最后，如果检测到某一个reads的kmer在后缀树中出现，则把该条reads保留。

​		由于采用的先过滤后拼接的策略，需要过滤的reads数据量相当大，需要拼接的reads则大幅度减少，这意味着拼接不再是计算瓶颈，过滤这一步骤才是本方法需要突破的核心。幸运的是，我们将过滤这一步骤的算法复杂读降低到O（n），大幅度压缩了时间成本。

## 6.2 拼接与校正：

### 6.2.1数据拼接

​       将过滤后的reads拼接为contigs，我们选择minia3作序列拼接软件。主要原因在于，区别于其他的一些序列拼接软件如velvet、spades、soapdenovo,minia3准确度高，速度快，消耗的内存更小，更加适用于短序列的拼接。

### 6.2.2数据校验

​      minia3可能产生多条contigs,而这些序列不都是我们需要的，因此GeneMiner对这些contigs做了进一步的处理。

**（1）长度校正**

​		根据参考的序列的平均长度，过滤掉长度太长或者太短的contigs。

**（2）方向校正**

​		GeneMiner利用BLAST+套件中的makeblastdb和balstn工具，将与参考序列的方向一致的contigs保留，与参考序列方向不一致的contigs做反向互补处理。

**（3）对齐剪切**

​		为了方便用户直接使用GeneMiner的结果，我们将contigs剪切，GeneMiner做了如下处理：（1）取出一条contig，记为序列a（2）取出一条参考序列，记为序列b（3）将序列a与序列b局部比对，记录高分值片段(HSP，High Scoring Pair)的起始位点与终止位点。（3）如果高分值片段的长度接近参考基因序列的长度，则写入gene_trimmed.fasta，否则不生成gene_trimmed.fasta文件（5）重复 步骤（1)~(4)，直到将所有contig对齐后剪切    

**（4）获得最佳结果**

​		如果xxx_trimmed.fasta中有多条序列，则将这些序列与参考序列进行两两比对，保留一致度最高且长度最长的一条contig，认为该contig为最佳结果

### **6.2.3**梯度逼近

​        GeneMiner通过不断迭代的方法，使结果不断优化，我们将这一行为称之为梯度逼近。

​		在这里，我们将用户输入参考序列到GeneMiner生成最佳结果这一完整的过程定义为一次迭代。在下一次的迭代过程中，参考序列将被上一次迭代生成的结果所取代，并重复之前的数据拼接和校验的过程。这个过程将被重复预定义的迭代次数，或者直到没有新的序列被发现。伴随着迭代次数的增加，能够匹配并组装reads增多，所以每经过一轮迭代，生成的序列通常会增长。迭代的方法不仅可以延申序列长度，保留更多的系统发育信息，同时还可以弥补参考序列的部分缺陷，特别是当参考序列相对远缘或参考序列中间存在空洞。通过迭代的方法，是我们逐渐获得最优解成为可能

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

GeneMiner : GeneMiner: a software for extracting phylogenetic markers from next generation sequencing data

