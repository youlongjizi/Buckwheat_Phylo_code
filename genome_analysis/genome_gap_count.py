import argparse
from Bio import SeqIO

def count_gaps(fasta_file, gap_char):
    gap_info = []  # 存储gap信息，每个元素为一个元组 (染色体, 左端位置, 右端位置, 长度)

    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        start_position = 0  # 当前序列的起始位置

        while start_position < len(sequence):
            # 查找gap字符的起始位置
            gap_start = sequence.find(gap_char, start_position)
            if gap_start == -1:
                break  # 没有找到更多的gap字符，退出循环
            
            # 查找gap字符的结束位置
            gap_end = gap_start
            while gap_end < len(sequence) and sequence[gap_end] == gap_char:
                gap_end += 1
            
            # 计算gap的长度
            gap_length = gap_end - gap_start
            
            # 记录gap信息
            gap_info.append((record.id, gap_start, gap_end - 1, gap_length))
            
            # 更新起始位置
            start_position = gap_end  # 移动到下一个位置

    return gap_info

def main():
    parser = argparse.ArgumentParser(description="Count gaps in a genome assembly.")
    parser.add_argument("input", help="Input FASTA file containing the genome assembly.")
    parser.add_argument("output", help="Output text file to save the gap statistics.")
    parser.add_argument("--gap_char", default='N', help="Character representing gaps in the sequence (default: 'N').")

    args = parser.parse_args()

    gap_info = count_gaps(args.input, args.gap_char)

    with open(args.output, 'w') as f:
        f.write("Chromosome\tLeft_Position\tRight_Position\tGap_Length\n")  # 输出表头
        for info in gap_info:
            f.write(f"{info[0]}\t{info[1]}\t{info[2]}\t{info[3]}\n")

if __name__ == "__main__":
    main()

