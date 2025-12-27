# Load thư viện cần thiết
# vcfR: Để đọc file VCF chuyên dụng
library(vcfR)
# ggplot2: Để vẽ biểu đồ
library(ggplot2)

# Nhận tham số từ dòng lệnh (Do Snakemake truyền vào)
# args[1] sẽ là file VCF
# args[2] sẽ là file PDF đầu ra
args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[1]
output_file <- args[2]

# In ra màn hình để kiểm tra
print(paste("Processing:", vcf_file))
print(paste("Results at", output_file))
