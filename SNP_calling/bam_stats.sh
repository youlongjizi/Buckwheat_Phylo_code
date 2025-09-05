#!/bin/bash

# 检查输入参数
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <BAM_DIRECTORY> <REFERENCE_GENOME_SIZE> <OUTPUT_FILE>"
    exit 1
fi

BAM_DIRECTORY=$1
REFERENCE_GENOME_SIZE=$2
OUTPUT_FILE=$3

# 初始化输出文件
echo -e "Filename\tTotal Reads\tMapped Reads\tMapping Rate\tAverage Depth\tCoverage Rate (>=1X)\tCoverage Rate (>=4X)\tCoverage Rate (>=5X)\tCoverage Rate (>=10X)\tCoverage Rate (>=20X)" > "$OUTPUT_FILE"

# 遍历文件夹中的每个 BAM 文件
for BAM_FILE in "$BAM_DIRECTORY"/*.bam; do
    # 检查文件是否存在
    if [ ! -f "$BAM_FILE" ]; then
        echo "No BAM files found in the directory."
        exit 1
    fi

    # 获取文件名
    FILENAME=$(basename "$BAM_FILE")

    # 计算总的 reads 数
    TOTAL_READS=$(samtools view -c "$BAM_FILE")

    # 计算 Mapped Reads 数
    MAPPED_READS=$(samtools view -c -F 4 "$BAM_FILE")

    # 计算比对率
    MAPPING_RATE=$(echo "scale=2; ($MAPPED_READS / $TOTAL_READS) * 100" | bc)

    # 使用 samtools depth 计算覆盖度
    samtools depth "$BAM_FILE" > depth.txt

    # 计算 Coverage Rates
    COVERAGE_RATE_1X=$(awk -v threshold=1 '$3 >= threshold {count++} END {printf "%.2f", (count / '"$REFERENCE_GENOME_SIZE"') }' depth.txt)
    COVERAGE_RATE_4X=$(awk -v threshold=4 '$3 >= threshold {count++} END {printf "%.2f", (count / '"$REFERENCE_GENOME_SIZE"') }' depth.txt)
    COVERAGE_RATE_5X=$(awk -v threshold=5 '$3 >= threshold {count++} END {printf "%.2f", (count / '"$REFERENCE_GENOME_SIZE"') }' depth.txt)
    COVERAGE_RATE_10X=$(awk -v threshold=10 '$3 >= threshold {count++} END {printf "%.2f", (count / '"$REFERENCE_GENOME_SIZE"') }' depth.txt)
    COVERAGE_RATE_20X=$(awk -v threshold=20 '$3 >= threshold {count++} END {printf "%.2f", (count / '"$REFERENCE_GENOME_SIZE"') }' depth.txt)

    # 计算平均深度
    AVERAGE_DEPTH=$(awk '{sum+=$3} END {printf "%.2f", sum/NR}' depth.txt)

    # 将结果写入输出文件
    echo -e "$FILENAME\t$TOTAL_READS\t$MAPPED_READS\t$MAPPING_RATE\t$AVERAGE_DEPTH\t$COVERAGE_RATE_1X\t$COVERAGE_RATE_4X\t$COVERAGE_RATE_5X\t$COVERAGE_RATE_10X\t$COVERAGE_RATE_20X\t$AVERAGE_DEPTH" >> "$OUTPUT_FILE"

    # 删除临时文件
    rm depth.txt
done

echo "Results have been written to $OUTPUT_FILE"
