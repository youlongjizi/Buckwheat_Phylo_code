#!/bin/bash

# 设置默认路径
MAFFT_DEFAULT="mafft"
TRIMAL_DEFAULT="trimAl"
FASTTREE_DEFAULT="FastTree"

# 检查输入参数
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 -i <input_fasta_file> -O <output_folder> [-mafft <path_to_mafft>] [-trimal <path_to_trimAl>]"
    exit 1
fi

# 解析输入参数
while getopts "i:O:mafft:trimal:" opt; do
  case $opt in
    i) input_file="$OPTARG"
    ;;
    O) output_folder="$OPTARG"
    ;;
    mafft) mafft_path="$OPTARG"
    ;;
    trimal) trimal_path="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
       exit 1
    ;;
  esac
done

# 如果未指定路径，则使用默认路径
mafft_path=${mafft_path:-$MAFFT_DEFAULT}
trimal_path=${trimal_path:-$TRIMAL_DEFAULT}
fasttree_path=${fasttree_path:-$FASTTREE_DEFAULT}

# 检查输出文件夹是否存在，如果不存在则创建
if [ ! -d "$output_folder" ]; then
  mkdir -p "$output_folder"
fi

# 输出的比对后的序列文件
aligned_file="${output_folder}/$(basename ${input_file%.*})_aligned.fasta"

# 输出的修剪后的序列文件
trimmed_file="${output_folder}/$(basename ${input_file%.*})_trimmed.fasta"

# 输出的树文件
tree_file="${output_folder}/$(basename ${input_file%.*}).treefile"

# 使用mafft进行序列比对
echo "Aligning sequences with mafft..."
$mafft_path --auto "$input_file" > "$aligned_file"

# 检查mafft是否成功执行
if [ $? -ne 0 ]; then
    echo "mafft failed to execute. Exiting."
    exit 1
fi

# 使用trimAl修剪序列
echo "Trimming conserved regions with trimAl..."
$trimal_path -in "$aligned_file" -out "$trimmed_file" -gt 0.8 -st 0.001 -cons 80

# 检查trimAl是否成功执行
if [ $? -ne 0 ]; then
    echo "trimAl failed to execute. Exiting."
    exit 1
fi

# 使用FastTree构建进化树
echo "Building phylogenetic tree with FastTree..."
$fasttree_path -nt "$trimmed_file" > "$tree_file"

# 检查FastTree是否成功执行
if [ $? -ne 0 ]; then
    echo "FastTree failed to execute. Exiting."
    exit 1
fi

# 打印完成消息
echo "The tree has been built and saved to $tree_file"

