for i in {1..19}
do 
    pop1=$(shuf -n 10 FD.sample.txt|sed ':a;N;s/\n/,/g;ta')   ##修改文件名
    pop2=$(shuf -n 10 FTW1.sample.txt|sed ':a;N;s/\n/,/g;ta')  ##修改文件名
    mkdir ./FD.data.${i}          ##修改文件名
    mkdir ./FTW1.data.${i}        ##修改文件名
    mkdir ./merge.data.${i}      ##修改文件名
    for k in {1..8}
    do
        smc++ vcf2smc --cores 30 --ignore-missing /data2/wild.genome_compare/DDX_analysis/population_analysis/SNP_population_dir/part_xizang_jinqiao_kuqiao_dir/filter_data/merge.jinqiao.kuqiao.snp.id.vcf.gz ./FD.data.${i}/FD.${k}.smc.gz chr0${k} FD:${pop1}     ##修改文件名
        smc++ vcf2smc --cores 30 --ignore-missing /data2/wild.genome_compare/DDX_analysis/population_analysis/SNP_population_dir/part_xizang_jinqiao_kuqiao_dir/filter_data/merge.jinqiao.kuqiao.snp.id.vcf.gz ./FTW1.data.${i}/FTW1.${k}.smc.gz chr0${k} FTW1:${pop2}     ##修改文件名
        smc++ vcf2smc --cores 30 --ignore-missing /data2/wild.genome_compare/DDX_analysis/population_analysis/SNP_population_dir/part_xizang_jinqiao_kuqiao_dir/filter_data/merge.jinqiao.kuqiao.snp.id.vcf.gz ./merge.data.${i}/FD_FTW1.${k}.smc.gz chr0${k} FD:${pop1} FTW1:${pop2}    ##修改文件名
        smc++ vcf2smc --cores 30 --ignore-missing /data2/wild.genome_compare/DDX_analysis/population_analysis/SNP_population_dir/part_xizang_jinqiao_kuqiao_dir/filter_data/merge.jinqiao.kuqiao.snp.id.vcf.gz ./merge.data.${i}/FTW1_FD.${k}.smc.gz chr0${k} FTW1:${pop2} FD:${pop1}   ##修改文件名
    done
    smc++ estimate --cores 30 -o ./FD.data.${i} 7e-9 ./FD.data.${i}/FD.*.smc.gz       ##修改文件名
    smc++ estimate --cores 30 -o ./FTW1.data.${i} 7e-9 ./FTW1.data.${i}/FTW1.*.smc.gz      ##修改文件名
    smc++ split -o FD_FTW1.split.${i}/ ./FD.data.${i}/model.final.json ./FTW1.data.${i}/model.final.json ./merge.data.${i}/*.smc.gz   ##修改文件名
    smc++ plot --csv --cores 20 FD_FTW1.split.${i}.pdf ./FD_FTW1.split.${i}/model.final.json       ##修改文件名
done

