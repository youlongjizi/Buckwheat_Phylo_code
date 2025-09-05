#!/bin/bash

# 解析命令行参数
while getopts "r:i:o:l:b:t:" opt; do
    case $opt in
        r) REFERENCE="$OPTARG" ;;
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        l) SAMPLE_LIST="$OPTARG" ;;
        b) BATCH_SIZE="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        *) echo "Unknown parameter passed: -$OPTARG"; exit 1 ;;
    esac
done

# 检查必需参数
if [[ -z "$REFERENCE" || -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$SAMPLE_LIST" || -z "$BATCH_SIZE" || -z "$THREADS" ]]; then
    echo "Usage: $0 -r <reference_genome> -i <input_directory> -o <output_directory> -l <sample_list> -b <batch_size> -t <threads>"
    exit 1
fi

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 读取样本列表
SAMPLES=()
while IFS=$'\t' read -r sample; do
    SAMPLES+=("$sample")
done < "$SAMPLE_LIST"

# 定义函数
find_fastq_files() {
    local sample=$1
    local fastq1=""
    local fastq2=""

    # 尝试匹配不同的文件命名模式
    for pattern in "_1.fq.gz" ".R1.fq.gz" "_1.fastq.gz" ".R1.fastq.gz"; do
        if [[ -f "${INPUT_DIR}/${sample}${pattern}" ]]; then
            fastq1="${INPUT_DIR}/${sample}${pattern}"
            fastq2="${fastq1/${pattern}/${pattern/1/2}}"
            break
        fi
    done

    if [[ -z "$fastq1" || -z "$fastq2" ]]; then
        echo "Error: Could not find FASTQ files for sample $sample"
        exit 1
    fi

    echo "$fastq1 $fastq2"
}

align_with_bwa() {
    local sample=$1
    read fastq1 fastq2 < <(find_fastq_files "$sample")
    local output_bam="${OUTPUT_DIR}/${sample}.sorted.bam"

    if [[ ! -f "$output_bam" ]]; then
        echo "Aligning and sorting BAM for sample $sample"
        bwa mem -t "$THREADS" -M -R "@RG\tID:$sample\tSM:$sample\tPL:illumina" "$REFERENCE" "$fastq1" "$fastq2" | \
        samtools sort -@ "$THREADS" -o "$output_bam"
    else
        echo "BAM file already exists for sample $sample, skipping alignment and sorting."
    fi
}

index_bam() {
    local sample=$1
    local bam_file="${OUTPUT_DIR}/${sample}.sorted.bam"
    local bai_file="${bam_file}.bai"
    local csi_file="${bam_file}.csi"

    if [[ ! -f "$bam_file" ]]; then
        echo "BAM file does not exist for sample $sample, skipping indexing."
        return
    fi

    if [[ ! -f "$bai_file" && ! -f "$csi_file" ]]; then
        echo "Indexing BAM file for sample $sample"
        samtools index "$bam_file"
    else
        echo "BAM index file already exists for sample $sample, skipping."
    fi
}

mark_duplicates() {
    local sample=$1
    local input_bam="${OUTPUT_DIR}/${sample}.sorted.bam"
    local output_bam="${OUTPUT_DIR}/${sample}.marked_duplicates.bam"
    local metrics_file="${OUTPUT_DIR}/${sample}.metrics.txt"

    if [[ ! -f "$output_bam" ]]; then
        echo "Marking duplicates for sample $sample"
        gatk MarkDuplicates -I "$input_bam" -O "$output_bam" -M "$metrics_file"
    else
        echo "Marked duplicates BAM file already exists for sample $sample, skipping."
    fi
}

index_marked_duplicates_bam() {
    local sample=$1
    local bam_file="${OUTPUT_DIR}/${sample}.marked_duplicates.bam"
    local bai_file="${bam_file}.bai"
    local csi_file="${bam_file}.csi"

    if [[ ! -f "$bam_file" ]]; then
        echo "Marked duplicates BAM file does not exist for sample $sample, skipping indexing."
        return
    fi

    if [[ ! -f "$bai_file" && ! -f "$csi_file" ]]; then
        echo "Indexing marked duplicates BAM file for sample $sample"
        samtools index "$bam_file"
    else
        echo "Marked duplicates BAM index file already exists for sample $sample, skipping."
    fi
}

call_variants() {
    local sample=$1
    local input_bam="${OUTPUT_DIR}/${sample}.marked_duplicates.bam"
    local output_gvcf="${OUTPUT_DIR}/${sample}.g.vcf.gz"

    # 确保 BAM 文件已索引
    index_marked_duplicates_bam "$sample"

    if [[ ! -f "$output_gvcf" ]]; then
        echo "Calling variants for sample $sample"
        gatk HaplotypeCaller -R "$REFERENCE" -I "$input_bam" -O "$output_gvcf" -ERC GVCF
    else
        echo "GVCF file already exists for sample $sample, skipping variant calling."
    fi
}

split_gvcfs_by_chromosome() {
    local sample=$1
    local input_gvcf="${OUTPUT_DIR}/${sample}.g.vcf.gz"
    local output_dir="${OUTPUT_DIR}/chromosomes/${sample}"

    mkdir -p "$output_dir"

    # 获取染色体列表
    chromosomes=$(bcftools index -s "$input_gvcf" | cut -f1)

    for chrom in $chromosomes; do
        local output_gvcf="${output_dir}/${sample}.${chrom}.g.vcf.gz"
        if [[ ! -f "$output_gvcf" ]]; then
            echo "Splitting chromosome $chrom for sample $sample"
            bcftools view -r "$chrom" "$input_gvcf" -Oz -o "$output_gvcf"
            tabix -p vcf "$output_gvcf"
        else
            echo "Chromosome $chrom gVCF file already exists for sample $sample, skipping."
        fi
    done
}

generate_chromosome_lists() {
    local list_dir="${OUTPUT_DIR}/chromosome_lists"
    mkdir -p "$list_dir"

    # 获取所有染色体列表
    chromosomes=$(bcftools index -s "${OUTPUT_DIR}/${SAMPLES[0]}.g.vcf.gz" | cut -f1)

    for chrom in $chromosomes; do
        local list_file="${list_dir}/${chrom}.list"
        > "$list_file"
        for sample in "${SAMPLES[@]}"; do
            local input_gvcf="${OUTPUT_DIR}/chromosomes/${sample}/${sample}.${chrom}.g.vcf.gz"
            if [[ -f "$input_gvcf" ]]; then
                echo "$input_gvcf" >> "$list_file"
                echo "Adding $input_gvcf to the list for chromosome $chrom"
            else
                echo "Warning: gVCF file for sample $sample and chromosome $chrom does not exist. Skipping."
            fi
        done
    done
}

combine_gvcfs_by_chromosome() {
    local list_dir="${OUTPUT_DIR}/chromosome_lists"
    local chromosomes=$(ls "$list_dir"/*.list | sed 's/\.list$//' | xargs -n1 basename)

    for chrom in $chromosomes; do
        local list_file="${list_dir}/${chrom}.list"
        local output_gvcf="${OUTPUT_DIR}/combined.${chrom}.g.vcf.gz"

        if [[ ! -f "$output_gvcf" ]]; then
            echo "Combining GVCF files for chromosome $chrom"
            gatk CombineGVCFs -R "$REFERENCE" --variant "$list_file" -O "$output_gvcf"
        else
            echo "Combined GVCF file for chromosome $chrom already exists, skipping."
        fi
    done
}

merge_combined_gvcfs() {
    local combined_gvcfs=()
    for chrom in $(bcftools index -s "${OUTPUT_DIR}/combined.chr1.g.vcf.gz" | cut -f1); do
        local input_gvcf="${OUTPUT_DIR}/combined.${chrom}.g.vcf.gz"
        if [[ -f "$input_gvcf" ]]; then
            combined_gvcfs+=("-V $input_gvcf")
            echo "Adding $input_gvcf to the list for final merging"
        else
            echo "Warning: Combined gVCF file for chromosome $chrom does not exist. Skipping."
        fi
    done

    if [[ ! -f "${OUTPUT_DIR}/final_combined.g.vcf.gz" ]]; then
        echo "Merging combined GVCF files"
        gatk CombineGVCFs -R "$REFERENCE" "${combined_gvcfs[@]}" -O "${OUTPUT_DIR}/final_combined.g.vcf.gz"
    else
        echo "Final combined GVCF file already exists, skipping."
    fi
}

joint_genotyping() {
    if [[ ! -f "${OUTPUT_DIR}/joint_calls.vcf.gz" ]]; then
        echo "Performing joint genotyping"
        gatk GenotypeGVCFs -R "$REFERENCE" -V "${OUTPUT_DIR}/final_combined.g.vcf.gz" -O "${OUTPUT_DIR}/joint_calls.vcf.gz"
    else
        echo "Joint calls VCF file already exists, skipping."
    fi
}

hard_filter_snps() {
    if [[ ! -f "${OUTPUT_DIR}/filtered_snps.vcf.gz" ]]; then
        echo "Selecting and filtering SNPs"
        gatk SelectVariants -R "$REFERENCE" -V "${OUTPUT_DIR}/joint_calls.vcf.gz" -select-type SNP -O "${OUTPUT_DIR}/joint_calls.snps.vcf.gz"
        gatk VariantFiltration -R "$REFERENCE" -V "${OUTPUT_DIR}/joint_calls.snps.vcf.gz" -O "${OUTPUT_DIR}/filtered_snps.vcf.gz" --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || ReadPosRankSum < -8.0" --filter-name "SNP_HARD_FILTER"
    else
        echo "Filtered SNPs VCF file already exists, skipping."
    fi
}

# 检查并创建参考基因组索引
index_reference() {
    local ref_prefix="${REFERENCE%.fa*}.fa"
    local index_files=("${ref_prefix}.amb" "${ref_prefix}.ann" "${ref_prefix}.bwt" "${ref_prefix}.pac" "${ref_prefix}.sa")
    local all_indices_exist=true

    for index_file in "${index_files[@]}"; do
        if [[ ! -f "$index_file" ]]; then
            all_indices_exist=false
            break
        fi
    done

    if ! $all_indices_exist; then
        echo "Creating BWA index for the reference genome..."
        bwa index "$REFERENCE"
    else
        echo "BWA index already exists for the reference genome, skipping."
    fi

    # 创建参考基因组的 fai 索引文件
    local fai_file="${REFERENCE}.fai"
    if [[ ! -f "$fai_file" ]]; then
        echo "Creating Fasta index file for the reference genome..."
        samtools faidx "$REFERENCE"
    else
        echo "Fasta index file already exists for the reference genome, skipping."
    fi

    # 创建参考基因组的 dict 文件
    local dict_file="${REFERENCE%.fa*}.dict"
    if [[ ! -f "$dict_file" ]]; then
        echo "Creating sequence dictionary for the reference genome..."
        gatk CreateSequenceDictionary -R "$REFERENCE" -O "$dict_file"
    else
        echo "Sequence dictionary file already exists for the reference genome, skipping."
    fi
}

# 创建参考基因组索引
index_reference

# 使用 GNU Parallel 并行处理每个样本
export -f align_with_bwa
export -f index_bam
export -f mark_duplicates
export -f call_variants
export -f split_gvcfs_by_chromosome

parallel -j "$THREADS" 'align_with_bwa {} && index_bam {} && mark_duplicates {} && call_variants {} && split_gvcfs_by_chromosome {}' ::: "${SAMPLES[@]}"

# 生成每个染色体的 gVCF 文件列表
generate_chromosome_lists

# 分染色体合并不同样本的 gVCF 文件
combine_gvcfs_by_chromosome

# 合并不同染色体的 gVCF 文件
merge_combined_gvcfs

# 联合基因型分析
joint_genotyping

# 对SNP进行硬过滤
hard_filter_snps

echo "SNP calling and filtering process completed. Results are in ${OUTPUT_DIR}"

