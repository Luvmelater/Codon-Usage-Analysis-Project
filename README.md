# 跨物种密码子使用偏好性分析 (Codon Usage Bias Analysis)

这是一个用于分析跨物种密码子使用偏好性（Codon Usage Bias, CUB）的生物信息学项目。该项目通过计算和分析六个代表性模式生物的全基因组编码区（CDS）序列，旨在量化其CUB模式，鉴定最优密码子，并揭示这些模式与物种进化关系及翻译选择压力之间的联系。

## 主要分析功能

- **数据处理**: 从FASTA格式的CDS序列文件中读取数据，并进行有效的预处理。
- **RSCU计算**: 计算每个物种全基因组的相对同义密码子使用度（RSCU），作为标准化的偏好性度量。
- **ENC计算**: 从头实现了Wright (1990)的经典公式，为每个基因独立计算其有效密码子数（ENC）。
- **最优密码子鉴定**: 基于ENC值作为基因表达水平的代理，通过比较高、低表达基因集的RSCU差异，鉴定每个物种的最优密码子。
- **可视化**:
  - 绘制RSCU热图，宏观展示CUB模式。
  - 基于RSCU值进行PCA和层次聚类分析，揭示物种进化关系。
  - 绘制带注释的热图，可视化最优密码子及其选择强度（ΔRSCU）。

## 数据集 (Dataset)

本研究分析的CDS序列数据均来源于Ensembl及Ensembl Genomes数据库。由于原始数据文件体积过大，未包含在本代码仓库中。使用者可从以下来源自行下载所需物种的CDS FASTA文件。

**操作说明**:
1.  下载对应的 `.fa.gz` 压缩文件。
2.  使用解压软件（如7-Zip, WinRAR等）将其解压，得到 `.fa` 文件。
3.  将解压后的文件重命名为 `物种名.fasta` (例如, `Homo_sapiens.fasta`)。
4.  将所有重命名后的 `.fasta` 文件放入项目根目录下的 `data` 文件夹内，即可运行脚本复现本研究。

---

###  人类 (Homo sapiens)
- **来源**: Ensembl
- **下载页面**: [https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cds/](https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cds/)
- **直接下载命令**:
  ```bash
  wget [https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz](https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz)
  ```

###  小鼠 (Mus musculus)
- **来源**: Ensembl
- **下载页面**: [https://ftp.ensembl.org/pub/current_fasta/mus_musculus/cds/](https://ftp.ensembl.org/pub/current_fasta/mus_musculus/cds/)
- **直接下载命令**:
  ```bash
  wget [https://ftp.ensembl.org/pub/current_fasta/mus_musculus/cds/Mus_musculus.GRCm39.cds.all.fa.gz](https://ftp.ensembl.org/pub/current_fasta/mus_musculus/cds/Mus_musculus.GRCm39.cds.all.fa.gz)
  ```

###  果蝇 (Drosophila melanogaster)
- **来源**: Ensembl Genomes (Metazoa)
- **下载页面**: [https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/metazoa/fasta/drosophila_melanogaster/cds/](https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/metazoa/fasta/drosophila_melanogaster/cds/)
- **直接下载命令**:
  ```bash
  wget [https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/metazoa/fasta/drosophila_melanogaster/cds/Drosophila_melanogaster.BDGP6.32.cds.all.fa.gz](https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/metazoa/fasta/drosophila_melanogaster/cds/Drosophila_melanogaster.BDGP6.32.cds.all.fa.gz)
  ```

###  秀丽隐杆线虫 (Caenorhabditis elegans)
- **来源**: Ensembl Genomes (Metazoa)
- **下载页面**: [https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/metazoa/fasta/caenorhabditis_elegans/cds/](https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/metazoa/fasta/caenorhabditis_elegans/cds/)
- **直接下载命令**:
  ```bash
  wget [https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/metazoa/fasta/caenorhabditis_elegans/cds/Caenorhabditis_elegans.WBcel235.cds.all.fa.gz](https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/metazoa/fasta/caenorhabditis_elegans/cds/Caenorhabditis_elegans.WBcel235.cds.all.fa.gz)
  ```

###  拟南芥 (Arabidopsis thaliana)
- **来源**: Ensembl Genomes (Plants)
- **下载页面**: [https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/plants/fasta/arabidopsis_thaliana/cds/](https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/plants/fasta/arabidopsis_thaliana/cds/)
- **直接下载命令**:
  ```bash
  wget [https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/plants/fasta/arabidopsis_thaliana/cds/Arabidopsis_thaliana.TAIR10.cds.all.fa.gz](https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/plants/fasta/arabidopsis_thaliana/cds/Arabidopsis_thaliana.TAIR10.cds.all.fa.gz)
  ```

###  酵母 (Saccharomyces cerevisiae)
- **来源**: Ensembl Genomes (Fungi)
- **下载页面**: [https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/fungi/fasta/saccharomyces_cerevisiae/cds/](https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/fungi/fasta/saccharomyces_cerevisiae/cds/)
- **直接下载命令**:
  ```bash
  wget [https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/fungi/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz](https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/fungi/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz)
  ```

## 项目文件结构

```
.
├── data/                     # 存放所有物种的CDS FASTA文件
│   ├── Homo_sapiens.fasta
│   └── ...
├── plots/                    # 脚本运行后，所有生成的图表会保存在这里
│   ├── figure1_RSCU_heatmap.png
│   └── ...
├── analysis.py               # 主分析脚本
└── README.md                 # 项目说明文件
```

## 环境要求与安装

本项目使用Python 3.x，并依赖于以下第三方库。您可以通过pip来安装它们：

```bash
pip install pandas seaborn matplotlib scikit-learn scipy biopython numpy
```

## 如何使用

1.  **准备数据**: 按照上方“数据集”部分的说明，下载、解压并重命名所需的FASTA文件，放入`data`文件夹。
2.  **运行脚本**:
    - 打开终端或命令行，进入项目根目录。
    - 执行以下命令：
      ```bash
      python analysis.py
      ```
3.  **查看结果**:
    - 程序运行过程中，会在终端打印分析进度和最优密码子总结。
    - 所有生成的图表（PNG格式）会自动保存在一个新建的 `plots` 文件夹中。

## 主要发现摘要

本分析证实，密码子使用模式与物种的系统发育关系高度吻合。哺乳动物（智人、小鼠）不仅紧密聚类，还表现出强烈的GC偏好和更高的翻译选择压力。这些发现支持了CUB是物种在进化过程中，为适应翻译效率等选择压力而形成的稳定分子特征的观点。
