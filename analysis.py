# -*- coding: utf-8 -*-
"""
================================================================================
跨物种密码子使用偏好性分析
================================================================================
本脚本旨在通过生物信息学方法，分析并比较多个物种的密码子使用偏好性（CUB）。

主要功能包括：
1. 从FASTA格式的CDS序列文件中读取数据。
2. 计算每个物种全基因组的相对同义密码子使用度（RSCU）。
3. 基于有效密码子数（ENC）作为基因表达水平的代理指标，通过比较高、低表达
   基因集的RSCU差异，鉴定每个物种的最优密码子。
4. 生成一系列可视化图表，包括RSCU热图、物种聚类图（PCA与层次聚类树）以及
   展示选择压力的最优密码子热图。

作者: YifanLiu
日期: 2025年6月20日
"""

import os
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio.SeqIO import parse
from collections import defaultdict
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage, dendrogram

# --- 1. 核心数据结构与计算模块 ---

# 标准遗传密码表，将密码子映射到氨基酸或终止符
CODON_MAP = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S',
    'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C',
    'TGA': '*', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P',
    'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
    'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A',
    'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G',
    'GGA': 'G', 'GGG': 'G',
}


def get_codon_dicts():
    """
    根据标准遗传密码表，生成两个辅助字典。

    Returns:
        tuple:
            - dict: 氨基酸到其同义密码子列表的映射。
            - dict: 密码子到其对应氨基酸的映射。
    """
    syn_codons = defaultdict(list)
    codon_to_aa = {}
    for codon, aa in CODON_MAP.items():
        codon_to_aa[codon] = aa
        if aa != '*':
            syn_codons[aa].append(codon)
    return dict(syn_codons), codon_to_aa


# 全局变量，存储密码子映射关系
SYNONYMOUS_CODONS, CODON_TO_AA = get_codon_dicts()

# 按简并度对氨基酸进行分类，用于ENC计算
DEGENERACY_MAP = {
    '2-fold': ['F', 'Y', 'C', 'H', 'Q', 'N', 'K', 'D', 'E'],
    '3-fold': ['I'],
    '4-fold': ['V', 'P', 'T', 'A', 'G'],
    '6-fold': ['L', 'S', 'R']
}


def calculate_enc_from_scratch(codon_counts):
    """
    根据 Wright (1990) 的经典公式，从头计算有效密码子数 (ENC)。
    ENC = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6

    Args:
        codon_counts (dict): 单个基因的密码子计数字典。

    Returns:
        float: 该基因的ENC值。
    """
    # 按氨基酸家族对密码子计数进行分组
    aa_codon_counts = defaultdict(dict)
    for codon, count in codon_counts.items():
        if count > 0 and CODON_TO_AA.get(codon, '*') != '*':
            aa = CODON_TO_AA[codon]
            aa_codon_counts[aa][codon] = count

    # 计算每个简并度家族的平均纯合度 (F_k)
    family_homozygosities = defaultdict(list)
    for aa, counts_dict in aa_codon_counts.items():
        total_codons_for_aa = sum(counts_dict.values())
        if total_codons_for_aa > 1:
            # 计算纯合度 F_a
            sum_of_sq_proportions = sum([(n / total_codons_for_aa) ** 2 for n in counts_dict.values()])
            homozygosity_F_a = (total_codons_for_aa * sum_of_sq_proportions - 1) / (total_codons_for_aa - 1)

            # 将F_a添加到对应的简并度列表中
            for deg, aa_list in DEGENERACY_MAP.items():
                if aa in aa_list:
                    family_homozygosities[deg].append(homozygosity_F_a)
                    break

    # 计算每个简并度家族的平均F值 (F_k)
    avg_F = {deg: np.mean(f_list) for deg, f_list in family_homozygosities.items() if f_list}

    # 根据公式计算ENC值，对缺失的家族进行容错处理
    enc = 2.0
    if '2-fold' in avg_F and avg_F['2-fold'] > 0: enc += 9.0 / avg_F['2-fold']
    if '3-fold' in avg_F and avg_F['3-fold'] > 0: enc += 1.0 / avg_F['3-fold']
    if '4-fold' in avg_F and avg_F['4-fold'] > 0: enc += 5.0 / avg_F['4-fold']
    if '6-fold' in avg_F and avg_F['6-fold'] > 0: enc += 3.0 / avg_F['6-fold']

    return enc


# --- 2. 数据处理与分析模块 ---

def calculate_rscu(codon_counts):
    """从密码子计数字典计算RSCU值。"""
    rscu_values = {}
    aa_counts = defaultdict(int)
    for codon, count in codon_counts.items():
        aa = CODON_TO_AA.get(codon, '*')
        if aa != '*':
            aa_counts[aa] += count

    for codon, count in codon_counts.items():
        aa = CODON_TO_AA.get(codon, '*')
        if aa != '*' and aa_counts.get(aa, 0) > 0:
            degeneracy = len(SYNONYMOUS_CODONS[aa])
            if degeneracy > 0:
                expected_count = aa_counts[aa] / degeneracy
                if expected_count > 0:
                    rscu_values[codon] = count / expected_count
    return rscu_values


def process_fasta_to_rscu(filepath):
    """
    读取FASTA文件，进行预处理，并计算全基因组的RSCU值。

    Args:
        filepath (str): FASTA文件的路径。

    Returns:
        tuple:
            - dict: 全基因组RSCU值字典。
            - list: 从文件中读取的所有合格CDS序列列表。
    """
    codon_counts = defaultdict(int)
    sequences = []
    for record in parse(filepath, "fasta"):
        seq = str(record.seq).upper()
        if len(seq) > 0 and len(seq) % 3 == 0 and all(c in 'ATCG' for c in seq):
            sequences.append(seq)
            for i in range(0, len(seq), 3):
                if len(seq[i:i + 3]) == 3:
                    codon_counts[seq[i:i + 3]] += 1

    genome_rscu = calculate_rscu(codon_counts)
    return genome_rscu, sequences


def determine_optimal_codons(sequences, num_genes=5):
    """
    通过比较高、低表达基因集的RSCU确定最优密码子。

    Args:
        sequences (list): 某物种的所有CDS序列列表。
        num_genes (int): 用于比较的高/低表达基因集的大小。

    Returns:
        dict: 包含最优密码子及其对应ΔRSCU值的字典。
    """
    if len(sequences) < 2 * num_genes:
        return {}

    # 为每个基因计算ENC值
    gene_enc_list = []
    for seq in sequences:
        codon_counts = defaultdict(int)
        for i in range(0, len(seq), 3):
            codon_counts[seq[i:i + 3]] += 1
        enc_val = calculate_enc_from_scratch(codon_counts)
        if enc_val > 0:
            gene_enc_list.append({'seq': seq, 'enc': enc_val})

    if len(gene_enc_list) < 2 * num_genes:
        return {}

    # 根据ENC值排序并筛选高/低表达基因集
    gene_enc_list.sort(key=lambda x: x['enc'])
    high_exp_seqs = [item['seq'] for item in gene_enc_list[:num_genes]]
    low_exp_seqs = [item['seq'] for item in gene_enc_list[-num_genes:]]

    # 分别计算高/低表达基因集的RSCU
    high_counts = defaultdict(int)
    for seq in high_exp_seqs:
        for i in range(0, len(seq), 3): high_counts[seq[i:i + 3]] += 1
    rscu_high = calculate_rscu(high_counts)

    low_counts = defaultdict(int)
    for seq in low_exp_seqs:
        for i in range(0, len(seq), 3): low_counts[seq[i:i + 3]] += 1
    rscu_low = calculate_rscu(low_counts)

    # 根据判定标准筛选最优密码子
    optimal_codons_info = {}
    for codon in rscu_high:
        if codon in rscu_low:
            delta_rscu = rscu_high[codon] - rscu_low[codon]
            if delta_rscu >= 0.08 and rscu_high[codon] > 1.0:
                optimal_codons_info[codon] = delta_rscu

    return optimal_codons_info


# --- 3. 可视化模块 ---

def plot_rscu_heatmap(rscu_df, output_dir):
    """绘制RSCU热图。"""
    codons_to_exclude = ['ATG', 'TGG']
    plot_df = rscu_df.drop(columns=codons_to_exclude, errors='ignore')
    plt.figure(figsize=(20, 8))
    sns.heatmap(plot_df, cmap='coolwarm', center=1.0)
    plt.title('Relative Synonymous Codon Usage (RSCU) Across Species', fontsize=16)
    plt.xlabel('Codon Triplets (excluding Met, Trp, Stop)', fontsize=12)
    plt.ylabel('Species', fontsize=12)
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    output_path = os.path.join(output_dir, 'figure1_RSCU_heatmap.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"RSCU heatmap saved to {output_path}")
    plt.show()


def plot_clustering_from_rscu(rscu_df, output_dir):
    """基于RSCU值绘制PCA图和层次聚类树。"""
    # PCA
    scaled_data = StandardScaler().fit_transform(rscu_df)
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(scaled_data)
    pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'], index=rscu_df.index)
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x='PC1', y='PC2', data=pca_df, s=100, hue=pca_df.index, legend=False)
    for i in range(pca_df.shape[0]):
        plt.text(x=pca_df['PC1'].iloc[i] + 0.1, y=pca_df['PC2'].iloc[i] + 0.1, s=pca_df.index[i])
    variance_ratio = pca.explained_variance_ratio_
    plt.xlabel(f'Principal Component 1 ({variance_ratio[0]:.2%})')
    plt.ylabel(f'Principal Component 2 ({variance_ratio[1]:.2%})')
    plt.title('Species Clustering based on RSCU (PCA)', fontsize=16)
    plt.grid()
    output_path_pca = os.path.join(output_dir, 'figure2_RSCU_pca_clustering.png')
    plt.savefig(output_path_pca, dpi=300)
    print(f"PCA plot based on RSCU saved to {output_path_pca}")
    plt.show()

    # Dendrogram
    linked = linkage(rscu_df, method='ward', metric='euclidean')
    plt.figure(figsize=(12, 7))
    dendrogram(
        linked,
        orientation='top',
        labels=rscu_df.index.to_list(),
        distance_sort='descending',
        leaf_rotation=30,
        leaf_font_size=10
    )
    plt.title('Hierarchical Clustering Dendrogram based on RSCU', fontsize=16)
    plt.ylabel('Distance (Ward)', fontsize=12)
    plt.xlabel('Species', fontsize=12)
    plt.tight_layout()
    output_path_dendro = os.path.join(output_dir, 'figure3_RSCU_dendrogram.png')
    plt.savefig(output_path_dendro, dpi=300)
    print(f"Dendrogram based on RSCU saved to {output_path_dendro}")
    plt.show()


def plot_optimal_codon_heatmap(optimal_codon_data, output_dir):
    """绘制最优密码子热图，背景颜色代表ΔRSCU。"""
    records = []
    for species, codons_dict in optimal_codon_data.items():
        for codon, delta_rscu in codons_dict.items():
            records.append({
                'species': species, 'codon': codon, 'amino_acid': CODON_TO_AA.get(codon, '?'), 'delta_rscu': delta_rscu
            })
    if not records:
        print("未找到可供绘图的最优密码子。")
        return

    long_df = pd.DataFrame(records)

    # 创建用于颜色和文字注释的两个透视表
    delta_rscu_pivot = long_df.pivot_table(index='amino_acid', columns='species', values='delta_rscu', aggfunc='mean')
    codon_text_pivot = long_df.pivot_table(index='amino_acid', columns='species', values='codon',
                                           aggfunc=lambda x: ', '.join(sorted(x)))

    aa_order = sorted(long_df['amino_acid'].unique())
    delta_rscu_pivot = delta_rscu_pivot.reindex(aa_order)
    codon_text_pivot = codon_text_pivot.reindex(aa_order)

    plt.figure(figsize=(14, 10))
    sns.heatmap(
        delta_rscu_pivot,
        annot=codon_text_pivot.fillna(''),
        fmt='s',
        cmap='viridis',
        linewidths=.5,
        linecolor='black',
        cbar_kws={'label': 'ΔRSCU (Average Selection Strength)'}
    )
    plt.title('Optimal Codons and Their Selection Strength', fontsize=16)
    plt.xlabel('Species', fontsize=12)
    plt.ylabel('Amino Acid', fontsize=12)
    plt.xticks(rotation=30, ha='right')
    plt.tight_layout()
    output_path = os.path.join(output_dir, 'figure4_optimal_codon_heatmap.png')
    plt.savefig(output_path, dpi=300)
    print(f"Optimal codons heatmap saved to {output_path}")
    plt.show()


# --- 4. 主程序入口 ---
def main():
    """
    主函数，执行整个分析流程。
    """
    # --- 配置参数 ---
    # 将包含FASTA文件的文件夹命名为 "data"，并与此脚本放在同一目录下
    DATA_FOLDER = 'data'
    # 所有生成的图表将保存在 "plots" 文件夹中
    OUTPUT_FOLDER = 'plots'

    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)

    fasta_files = glob.glob(os.path.join(DATA_FOLDER, '*.fa'))
    if not fasta_files:
        print(f"错误：在'{DATA_FOLDER}'中未找到.fasta文件。请检查文件路径和后缀。")
        return

    all_species_rscu = {}
    all_optimal_codons = {}

    print("--- 开始分析 ---")
    for f in fasta_files:
        species_name = os.path.basename(f).split('.')[0].replace('_', ' ').title()
        print(f"\n正在处理物种: {species_name}...")

        # 计算RSCU和最优密码子
        genome_rscu, sequences = process_fasta_to_rscu(f)
        all_species_rscu[species_name] = genome_rscu

        optimal_codons_info = determine_optimal_codons(sequences)
        all_optimal_codons[species_name] = optimal_codons_info
        print(f"找到 {len(optimal_codons_info)} 个最优密码子。")

    # 准备用于绘图的DataFrame
    rscu_df = pd.DataFrame(all_species_rscu).T.fillna(0)
    sense_codons = [c for c, a in CODON_MAP.items() if a != '*']
    # 确保所有密码子列都存在
    for codon in sense_codons:
        if codon not in rscu_df.columns:
            rscu_df[codon] = 0
    rscu_df = rscu_df[sorted(sense_codons)]

    # --- 生成可视化图表 ---
    print("\n--- 生成可视化图表 ---")
    plot_rscu_heatmap(rscu_df, OUTPUT_FOLDER)
    plot_clustering_from_rscu(rscu_df, OUTPUT_FOLDER)
    plot_optimal_codon_heatmap(all_optimal_codons, OUTPUT_FOLDER)

    # --- 打印最终总结 ---
    print("\n--- 分析完成 ---")
    print("最优密码子总结:")
    for species, codons_dict in all_optimal_codons.items():
        codons_str = ', '.join(sorted(codons_dict.keys()))
        print(f"- {species} ({len(codons_dict)}): {codons_str}")


if __name__ == '__main__':
    main()