#!/bin/bash
# process_mcmctree_perl_fixed.sh

# 默认倍数
MULTIPLIER=100

show_usage() {
    echo "用法: $0 [选项] <输入文件> <输出文件>"
    echo ""
    echo "选项:"
    echo "  -m, --multiplier NUM   时间节点乘数 (默认: 100)"
    echo "  -h, --help            显示此帮助信息"
    echo ""
    echo "示例:"
    echo "  $0 FigTree.tre tree.txt"
    echo "  $0 -m 1000 FigTree.tre tree.txt"
}

# 解析参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -m|--multiplier)
            MULTIPLIER="$2"
            shift 2
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        -*)
            echo "错误: 未知选项 $1"
            show_usage
            exit 1
            ;;
        *)
            break
            ;;
    esac
done

if [ $# -ne 2 ]; then
    echo "错误: 需要输入文件和输出文件参数"
    show_usage
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="$2"

if [ ! -f "$INPUT_FILE" ]; then
    echo "错误: 输入文件 '$INPUT_FILE' 不存在"
    exit 1
fi

# 修正的perl代码
grep "UTREE 1 =" "$INPUT_FILE" | \
sed -E -e "s/\[[^]]*\]//g" -e "s/[ \t]//g" -e "/^$/d" -e "s/UTREE1=//" | \
perl -pe '
    s/:([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)/
        my $num = $1;
        my $result = $num * '"$MULTIPLIER"';
        if($result == int($result)) {
            ":" . sprintf("%.0f", $result);
        } else {
            my $formatted = sprintf("%.6f", $result);
            $formatted =~ s\/0+$\/\/;
            $formatted =~ s\/\.$\/\/;
            ":" . $formatted;
        }
    /ge;
' > "$OUTPUT_FILE"

echo "处理完成: $INPUT_FILE -> $OUTPUT_FILE (倍数: $MULTIPLIER)"

