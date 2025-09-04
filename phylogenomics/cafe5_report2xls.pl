#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

# 定义变量
my ($species_name, $orthogroups_file, $cafe_file, $output_file, $help);

# 获取命令行参数
GetOptions(
    'species|s=s'     => \$species_name,
    'orthogroups|g=s' => \$orthogroups_file,
    'cafe|c=s'        => \$cafe_file,
    'output|o=s'      => \$output_file,
    'help|h'          => \$help,
) or die "Error in command line arguments\n";

# 显示帮助信息
if ($help || !$species_name || !$orthogroups_file || !$cafe_file || !$output_file) {
    print <<EOF;
Usage: $0 [options]

Options:
    -s, --species      Species name to analyze
    -g, --orthogroups  Orthogroups file path
    -c, --cafe         CAFE5 output file path
    -o, --output       Output file path
    -h, --help         Show this help message

Example:
    perl $0 -s species_name -g orthogroups.txt -c cafe5_output.txt -o expansion_decrease.xls

EOF
    exit;
}

my (%hash, $orthologs_join);

# 读取orthogroups文件
open (IN1, $orthogroups_file) || die "Cannot open $orthogroups_file: $!";
while(<IN1>) {
    chomp;
    next if /^\s*$/; # 跳过空行
    $orthologs_join = "";
    
    my ($orthologs, $id) = split(/:\s*/, $_, 2);
    next unless defined $id; # 确保有ID部分
    
    my @d = split(/\s+/, $id);
    foreach my $i (@d) {
        if ($i =~ /$species_name/) {
            $orthologs_join = join(' ', $orthologs_join, $i);
        }
    }
    $orthologs_join =~ s/^\s+//g; # 去除开头空格
    $hash{$orthologs} = $orthologs_join if $orthologs_join;
}
close IN1;

# 打开输出文件
open (OUT, ">$output_file") || die "Cannot create output file $output_file: $!";

my ($species, $Progenitor);
my $header_printed = 0;

# 读取CAFE5输出文件
open (IN2, $cafe_file) || die "Cannot open $cafe_file: $!";
while(<IN2>) {
    chomp;
    next if /^\s*$/; # 跳过空行
    next if /^Tree/;
    next if /^Lambda/;
    next if /^# IDs of nodes/;
    next if /^# Output format for/;
    next if /^Average Expansion/;
    next if /^Expansion/;
    next if /^nRemain/;
    next if /^nDecrease/;
    
    my @a = split(/\t/, $_);
    next unless @a >= 3; # 确保至少有3列
    
    # 打印表头
    if ($a[2] =~ /Family-wide P-value/ || $a[2] =~ /P-value/ || !$header_printed) {
        if (!$header_printed) {
            print OUT "ID\tExpansion_Decrease\tn$species_name\tnProgenitor\tFamily-wide P-value\tgeneid\n";
            $header_printed = 1;
        }
        next if $a[2] =~ /P-value/; # 跳过表头行
    }
    
    # 解析物种数据 - 适配CAFE5格式
    $species = undef;
    $Progenitor = undef;
    
    # 尝试多种可能的CAFE5输出格式
    if ($a[1] =~ /\b${species_name}_(\d+):(\d+)/) {
        $species = $1;
        # 寻找祖先节点信息
        if ($a[1] =~ /\)_(\d+)/) {
            $Progenitor = $1;
        }
    }
    
    # 如果上述方法失败，尝试其他格式
    if (!defined $species) {
        # 尝试匹配格式: species_number:count
        if ($a[1] =~ /\b${species_name}_?(\d+):?(\d+)/) {
            $species = $2; # 使用count作为species数量
            # 尝试提取祖先信息
            if ($a[1] =~ /<(\d+)>/ || $a[1] =~ /\((\d+)\)/ || $a[1] =~ /_(\d+)$/) {
                $Progenitor = $1;
            }
        }
    }
    
    # 如果仍然无法解析，尝试简单的数字提取
    if (!defined $species && $a[1] =~ /$species_name/) {
        my @numbers = ($a[1] =~ /(\d+)/g);
        if (@numbers >= 2) {
            $species = $numbers[-2]; # 倒数第二个数字
            $Progenitor = $numbers[-1]; # 最后一个数字
        }
    }
    
    # 检查是否在hash中存在该基因家族
    if (defined $species && defined $Progenitor && exists $hash{$a[0]}) {
        # 确保p-value是数值且显著
        my $pvalue = $a[2];
        if ($pvalue =~ /^[\d\.e\-\+]+$/ && $pvalue < 0.05) {
            if ($species > $Progenitor) {
                print OUT "$a[0]\tExpansion\t$species\t$Progenitor\t$a[2]\t$hash{$a[0]}\n";
            } elsif ($species < $Progenitor) {
                print OUT "$a[0]\tDecrease\t$species\t$Progenitor\t$a[2]\t$hash{$a[0]}\n";
            }
        }
    }
}
close IN2;
close OUT;

# 输出统计信息
if (!$header_printed) {
    warn "Warning: No valid data found. Please check:\n";
    warn "1. Species name '$species_name' matches the format in CAFE5 output\n";
    warn "2. CAFE5 output file format is correct\n";
    warn "3. Orthogroups file format is correct\n";
} else {
    print "Analysis completed successfully!\n";
    print "Output file: $output_file\n";
}
