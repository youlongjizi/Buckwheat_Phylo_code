# 加载必要的库
library(optparse)
library(tidyverse)
library(bedtoolsr)
library(openxlsx)

# 定义命令行选项
option_list <- list(
  make_option(c("-t", "--type"), type = "character", default = "fst",
              help = "Type of annotation to perform: 'fst' or 'pi'", metavar = "character"),
  make_option(c("-g", "--gene_file"), type = "character", default = NULL,
              help = "Path to the gene position annotation file", metavar = "character"),
  make_option(c("-i", "--input_file"), type = "character", default = NULL,
              help = "Path to the input FST or Pi data file", metavar = "character"),
  make_option(c("-o", "--output_prefix"), type = "character", default = "output",
              help = "Prefix for output files", metavar = "character")
)

# 解析命令行参数
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 读取基因位置注释文件
gene <- read.delim(opt$gene_file, header = TRUE, sep = "\t")

# 定义函数用于FST注释
perform_fst_annotation <- function(fst_data, gene, prefix) {
  compare <- fst_data %>% select('pop1', 'pop2') %>% distinct()
  
  for (k in 1:nrow(compare)) {
    P1 <- compare[k, 1]
    P2 <- compare[k, 2]
    fst <- fst_data %>% filter(pop1 == P1, pop2 == P2)
    fst_top5 <- fst %>% slice_max(avg_wc_fst, prop = 0.05) %>%
      select(chromosome, window_pos_1, window_pos_2, avg_wc_fst) %>% 
      set_names(c("Chr", "Start", "End", "Fst")) %>%
      arrange(Chr, Start, End)
    fst_top5_merge <- bt.merge(i = fst_top5, c = 4, o = "max")
    fst_top5_merge_anno <- bt.intersect(fst_top5_merge, gene, loj = TRUE) %>%
      set_names(c("Chr", "Start", "End", "Fst", names(gene))) %>%
      filter(!is.na(get(names(gene)[1])))
    
    write.csv(fst, file = paste0(prefix, "_", P1, "_", P2, ".all.csv"), row.names = FALSE)
    write.table(fst_top5_merge_anno, file = paste0(prefix, "_", P1, "_", P2, ".top5.anno.csv"), 
                row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    write.xlsx(fst_top5_merge_anno, file = paste0(prefix, "_", P1, "_", P2, ".top5.anno.xlsx"))
  }
}

# 定义函数用于Pi注释
perform_pi_annotation <- function(pi_data, gene, prefix) {
  compare <- pi_data %>% select(pop1 = 1, pop2 = 2) %>% distinct()
  
  for (k in 1:nrow(compare)) {
    p1 <- compare[k, 1]
    p2 <- compare[k, 2]
    
    prefix_out <- paste0(prefix, "_", p1, "_", p2, ".pi")
    
    Pi_diff <- pi_data %>% select(chromosome, window_pos_1, window_pos_2, !!p1, !!p2) %>%
      set_names(c("Chr", "Start", "End", "pop1_pi", "pop2_pi")) %>%
      mutate(pi_diff = abs(pop1_pi - pop2_pi))
    
    Pi_diff_top5 <- Pi_diff %>% slice_max(pi_diff, prop = 0.05) %>% arrange(Chr, Start)
    
    write.table(Pi_diff, file = paste0(prefix_out, ".all.csv"), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    write.csv(Pi_diff_top5, file = paste0(prefix_out, ".top5.csv"), row.names = FALSE)
    
    Pi_top5_merge <- bt.merge(i = Pi_diff_top5, c = 6, o = "max") %>% arrange(Chr, Start)
    Pi_top5_merge_anno <- bt.intersect(Pi_top5_merge, gene, loj = TRUE) %>%
      set_names(c("Chr", "Start", "End", "pi_diff", names(gene))) %>%
      filter(!is.na(get(names(gene)[1])))
    
    write.xlsx(Pi_top5_merge_anno, file = paste0(prefix_out, ".top5.anno.xlsx"))
  }
}

# 执行根据选项的不同逻辑
if (opt$type == "fst") {
  pixy <- read.delim(opt$input_file, header = TRUE, sep = "\t")
  perform_fst_annotation(pixy, gene, opt$output_prefix)
} else if (opt$type == "pi") {
  pi_data <- read.delim(opt$input_file, header = TRUE, sep = "\t") %>%
    select(pop, chromosome, window_pos_1, window_pos_2, avg_pi) %>%
    pivot_wider(names_from = pop, values_from = avg_pi)
  perform_pi_annotation(pi_data, gene, opt$output_prefix)
} else {
  stop("Unknown type: ", opt$type)
}

