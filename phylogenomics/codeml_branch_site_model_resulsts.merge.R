#!/usr/bin/env Rscript

# 设置CRAN镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# 加载必要的包
if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}

library(readr)

# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)

# 检查参数数量
if (length(args) != 4) {
  stop("Usage: Rscript sim.R -d <directory> -o <output_file>")
}

# 解析参数
dir <- args[2]
output <- args[4]

# 打印参数，用于调试
print(dir)
print(output)

# 获取目录列表
dirs <- list.files(dir, full.names = TRUE)

# 初始化一个空的数据框来收集数据
df <- data.frame()

# 遍历目录
for (dir_path in dirs) {
  # 获取目录中的.tsv文件
  files <- list.files(dir_path, pattern = "\\.tsv$", full.names = TRUE)
  
  # 遍历文件
  for (file in files) {
    # 读取文件
    temp_df <- read.delim(file,header=T,sep="\t")
    
    # 确保所有列的类型一致
    if (nrow(df) == 0) {
      df <- temp_df
    } else {
      # 使用rbind合并数据框
      df <- rbind(df, temp_df)
    }
  }
}

# 写入结果到文件
write_tsv(df, output)

