#!/usr/bin/env python3
"""
处理featureCounts输出文件的脚本
将其拆分为：基因长度文件、基因count文件和FPKM文件
"""

import sys
import pandas as pd
import numpy as np
import argparse
import os

def calculate_fpkm(counts, gene_lengths, total_reads):
    """
    计算FPKM值
    FPKM = (10^9 * C) / (N * L)
    其中：
    C = 比对到该基因的reads数
    N = 总的比对reads数（以百万为单位）
    L = 基因长度（以kb为单位）
    """
    # 避免除零错误
    gene_lengths_kb = gene_lengths / 1000.0
    total_reads_millions = total_reads / 1e6
    
    fpkm = (counts * 1e9) / (total_reads_millions * gene_lengths_kb * 1e9)
    return fpkm

def process_featurecounts(input_file, output_prefix):
    """
    处理featureCounts输出文件
    """
    # 读取文件
    print(f"正在读取文件: {input_file}")
    
    # featureCounts输出文件的前两行是注释，跳过
    df = pd.read_csv(input_file, sep='\t', comment='#')
    
    # 提取基因信息列
    gene_info_cols = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']
    gene_info = df[gene_info_cols].copy()
    
    # 提取样本的count列（从第7列开始）
    sample_cols = [col for col in df.columns if col not in gene_info_cols]
    count_data = df[['Geneid'] + sample_cols].copy()
    
    # 创建基因长度文件
    gene_length_file = f"{output_prefix}_gene_lengths.txt"
    print(f"正在创建基因长度文件: {gene_length_file}")
    gene_lengths = df[['Geneid', 'Length']].copy()
    gene_lengths.to_csv(gene_length_file, sep='\t', index=False)
    
    # 创建count文件
    count_file = f"{output_prefix}_counts.txt"
    print(f"正在创建count文件: {count_file}")
    count_data.to_csv(count_file, sep='\t', index=False)
    
    # 计算FPKM
    fpkm_file = f"{output_prefix}_fpkm.txt"
    print(f"正在计算FPKM并创建文件: {fpkm_file}")
    
    fpkm_data = pd.DataFrame()
    fpkm_data['Geneid'] = df['Geneid']
    
    # 对每个样本计算FPKM
    for sample_col in sample_cols:
        # 计算该样本的总reads数
        total_reads = df[sample_col].sum()
        print(f"样本 {sample_col} 的总reads数: {total_reads:,}")
        
        # 计算FPKM
        fpkm_values = calculate_fpkm(df[sample_col], df['Length'], total_reads)
        fpkm_data[sample_col] = fpkm_values
    
    # 保存FPKM文件
    fpkm_data.to_csv(fpkm_file, sep='\t', index=False, float_format='%.3f')
    
    # 生成汇总统计
    summary_file = f"{output_prefix}_summary_stats.txt"
    print(f"\n正在生成汇总统计: {summary_file}")
    
    with open(summary_file, 'w') as f:
        f.write("样本统计信息\n")
        f.write("=" * 50 + "\n")
        f.write(f"总基因数: {len(df)}\n")
        f.write(f"样本数: {len(sample_cols)}\n")
        f.write(f"样本名称: {', '.join(sample_cols)}\n\n")
        
        f.write("各样本reads统计:\n")
        for sample_col in sample_cols:
            total_reads = df[sample_col].sum()
            mapped_genes = (df[sample_col] > 0).sum()
            f.write(f"{sample_col}:\n")
            f.write(f"  总reads数: {total_reads:,}\n")
            f.write(f"  有reads的基因数: {mapped_genes:,}\n")
            f.write(f"  占总基因数比例: {mapped_genes/len(df)*100:.2f}%\n\n")
    
    print("\n处理完成！")
    print(f"生成的文件：")
    print(f"  - {gene_length_file}")
    print(f"  - {count_file}")
    print(f"  - {fpkm_file}")
    print(f"  - {summary_file}")

def main():
    parser = argparse.ArgumentParser(
        description='处理featureCounts输出文件，生成基因长度、count和FPKM文件',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  python process_featurecounts.py counts.txt -o results
  python process_featurecounts.py counts.txt -o results/project1
        """
    )
    
    parser.add_argument('input_file', 
                        help='featureCounts输出文件')
    parser.add_argument('-o', '--output_prefix', 
                        default='output',
                        help='输出文件前缀 (默认: output)')
    
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.input_file):
        print(f"错误: 输入文件 '{args.input_file}' 不存在！")
        sys.exit(1)
    
    # 处理文件
    try:
        process_featurecounts(args.input_file, args.output_prefix)
    except Exception as e:
        print(f"错误: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()

