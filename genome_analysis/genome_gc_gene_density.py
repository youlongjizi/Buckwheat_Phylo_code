import argparse
from Bio import SeqIO
from collections import defaultdict

def calculate_gc_content(sequence):
    """计算给定序列的GC含量"""
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    return (g_count + c_count) / len(sequence) if len(sequence) > 0 else 0

def parse_gff(gff_file):
    """从GFF文件中解析基因的位置及其所在的染色体"""
    genes = defaultdict(list)
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):  # 跳过注释行
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            chrom = parts[0]
            feature_type = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            if feature_type == "gene":
                genes[chrom].append((start, end))
    return genes

def calculate_gene_density(sequence, genes, window_size):
    """计算GC含量和基因密度"""
    length = len(sequence)
    gc_contents = []
    gene_densities = []
    
    # 统计窗口内的基因数量
    window_gene_count = [0] * ((length // window_size) + 1)
    
    for start, end in genes:
        for i in range(start // window_size, (end - 1) // window_size + 1):
            if i < len(window_gene_count):
                window_gene_count[i] += 1
    
    # 计算每个窗口的GC含量和基因密度
    for i in range(len(window_gene_count)):
        window_start = i * window_size
        window_end = min(window_start + window_size, length)
        window_sequence = sequence[window_start:window_end]
        
        gc_content = calculate_gc_content(window_sequence)
        gene_count = window_gene_count[i]
        
        gc_contents.append(gc_content)
        gene_densities.append(gene_count)

    return gc_contents, gene_densities

def main():
    parser = argparse.ArgumentParser(description='计算基因组移动窗口的GC含量和基因密度')
    parser.add_argument('-g', '--genome', required=True, help='基因组序列文件（FASTA格式）')
    parser.add_argument('-f', '--annotation', required=True, help='基因组注释文件（GFF格式）')
    parser.add_argument('-w', '--window_size', type=int, required=True, help='窗口大小')
    parser.add_argument('-o', '--output_prefix', required=True, help='输出文件前缀')

    args = parser.parse_args()

    # 读取基因组序列
    genome_seq = {}
    for record in SeqIO.parse(args.genome, "fasta"):
        genome_seq[record.id] = str(record.seq)

    # 解析GFF文件以获取基因位置及其所在的染色体
    genes = parse_gff(args.annotation)

    # 处理每个染色体
    for chrom, chrom_genes in genes.items():
        if chrom not in genome_seq:
            print(f"警告：染色体 {chrom} 在基因组序列文件中未找到，跳过该染色体。")
            continue

        sequence = genome_seq[chrom]

        # 计算GC含量和基因密度
        gc_contents, densities = calculate_gene_density(sequence, chrom_genes, args.window_size)

        # 输出GC含量结果
        with open(f"{args.output_prefix}_gc.txt", 'a') as f:
            for i in range(len(gc_contents)):
                window_start = i * args.window_size
                window_end = min(window_start + args.window_size, len(sequence))
                f.write(f"{chrom}\t{window_start}\t{window_end}\t{gc_contents[i]:.4f}\n")

        # 输出基因密度结果
        with open(f"{args.output_prefix}_gden.txt", 'a') as f:
            for i in range(len(densities)):
                window_start = i * args.window_size
                window_end = min(window_start + args.window_size, len(sequence))
                f.write(f"{chrom}\t{window_start}\t{window_end}\t{densities[i]}\n")

if __name__ == '__main__':
    main()

