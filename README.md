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

1.  **准备数据**:
    - 在项目根目录下创建一个名为 `data` 的文件夹。
    - 将所有物种的CDS序列FASTA文件（以`.fasta`为后缀）放入`data`文件夹中。
2.  **运行脚本**:
    - 打开终端或命令行，进入项目根目录。
    - 执行以下命令：
      ```bash
      python analysis.py
      ```
3.  **查看结果**:
    - 程序运行过程中，会在终端打印分析进度和最优密码子总结。
    - 所有生成的图表（PNG格式）会自动保存在一个新建的 `plots` 文件夹中。
# Codon-Usage-Analysis-Project
# Codon-Usage-Analysis-Project
