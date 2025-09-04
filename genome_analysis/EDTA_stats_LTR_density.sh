#!/bin/bash

# 默认参数
genome=""
TE_anno=""
window_size_kb=1000
output_prefix=""

# 处理命令行参数
while getopts "g:f:w:o:" opt; do
  case ${opt} in
    g ) genome=${OPTARG};;
    f ) TE_anno=${OPTARG};;
    w ) window_size_kb=${OPTARG};;
    o ) output_prefix=${OPTARG};;
    \? ) echo "Invalid option: -${OPTARG}" >&2; exit 1;;
    : ) echo "Option -${OPTARG} requires an argument." >&2; exit 1;;
  esac
done

# 检查必需参数
if [[ -z "$genome" || -z "$TE_anno" || -z "$output_prefix" ]]; then
  echo "Usage: $0 -g <genome_file> -f <TE_annotation_gff3> -w <window_size_kb> -o <output_prefix>"
  exit 1
fi

# 将窗口大小从 Kb 转换为 bp
window_size=$(($window_size_kb * 1000))

# 提取前缀
prefix=$(basename "${genome}" .fa)
if [[ "$prefix" == "${genome}" ]]; then
  prefix=$(basename "${genome}" .fasta)
fi

# 生成染色体大小文件
samtools faidx "${genome}"
awk '{print $1"\t"$2}' "${genome}.fai" > "${output_prefix}.txt"

# 生成滑动窗口文件 (bed3 和 bed6)
bedtools makewindows -g "${output_prefix}.txt" -w "${window_size}" -s "${window_size}" > "${output_prefix}_sw.bed3"
awk '{print $1"\t"$2"\t"$3"\t"NR"\t.\t."}' "${output_prefix}_sw.bed3" > "${output_prefix}_sw.bed6"

# 统计 Copia LTR 转座子的密度
awk '$3=="Copia_LTR_retrotransposon"' "${TE_anno}" | \
  bedtools coverage -a "${output_prefix}_sw.bed6" -b - -counts -F 0.5|awk '{print $1,$2,$3,$NF}' > "${output_prefix}_${window_size_kb}kb_wd.Copia.density.txt"

# 统计 Gypsy LTR 转座子的密度
awk '$3=="Gypsy_LTR_retrotransposon"' "${TE_anno}" | \
  bedtools coverage -a "${output_prefix}_sw.bed6" -b - -counts -F 0.5|awk '{print $1,$2,$3,$NF}' > "${output_prefix}_${window_size_kb}kb_wd.Gypsy.density.txt"

# 输出完成信息
echo "Copia LTR density calculated and saved to ${output_prefix}_sw.Copia.density.txt"
echo "Gypsy LTR density calculated and saved to ${output_prefix}_sw.Gypsy.density.txt"
