##parameter
SAMPLE=$1
Reference=$2
Gff=$3
raw_data="./raw_fastq"
filter_data="./clean_fastq_dir"

# 确保输出目录存在
mkdir -p ${raw_data} ${filter_data}

##build index
if [ -f "${Reference}.fai" ];then
    echo -e "${Reference}.fai 文件存在，skip index genome"
else
    samtools faidx $Reference
fi


if [ -f "${Reference%.*}.1.ht2" ];then
    echo -e "skip hisat index ${Reference}"
else
    /data/software/RNA-seq/hisat2-2.2.1/hisat2-build $Reference ${Reference%.*}
fi

##gtf
if [ -f "${Gff%.*}.gtf" ];then
    echo  "skip build gtf"
else
    gffread ${Gff} -T -o ${Gff%.*}.gtf
fi

##filter
fastp -i ${raw_data}/${SAMPLE}_1.fastq.gz -I ${raw_data}/${SAMPLE}_2.fastq.gz -w 20 -o ${filter_data}/${SAMPLE}_R1_filtered.fq.gz -O ${filter_data}/${SAMPLE}_R2_filtered.fq.gz


##align
echo "将RNA-seq的测序reads使用hisat2比对到参考基因租组-- ${SAMPLE}"
/data/software/RNA-seq/hisat2-2.2.1/hisat2 \
    -q -x ${Reference%.*} \
    -1 ${filter_data}/${SAMPLE}_R1_filtered.fq.gz \
    -2 ${filter_data}/${SAMPLE}_R2_filtered.fq.gz \
    -S ${SAMPLE}.sam  -p 20 -t --bowtie2-dp 2 ;
samtools sort -@ 20 -o ${SAMPLE}.bam ${SAMPLE}.sam ;

rm ${SAMPLE}.sam
