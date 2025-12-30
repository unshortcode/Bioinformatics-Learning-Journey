# Load thư viện cần thiết
# vcfR: Để đọc file VCF chuyên dụng
library(vcfR)

# Nhận tham số từ dòng lệnh (Do Snakemake truyền vào)
# args[1] sẽ là file VCF
# args[2] sẽ là file PDF đầu ra
args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[1]
output_file <- args[2]

# In ra màn hình để kiểm tra
print(paste("Processing:", vcf_file))

# 1. Đọc file VCF
vcf <- read.vcfR(vcf_file, verbose = FALSE)

# 2. Đọc file fasta của reference genome
seq_ref <- ape::read.dna("refs/ecoli_ref.fasta", format = "fasta")

# 3. Tạo chrom object
chrom <- create.chromR(
  name = "E. Coli Analysis (filtered)",
  vcf = vcf,
  seq = seq_ref,
  verbose = FALSE
)

# 4. Tạo mask để chọn lọc dữ liệu hiển thị
chrom_mask <- masker(
  chrom,
  min_QUAL = 0,
  min_DP = 0,
  max_DP = 150,
  min_MQ = 30,
  max_MQ = 60.5
)

# Mở file PDF
pdf(output_file)

# In biểu đồ vào file pdfpdf
plot(chrom_mask)
chromoqc(chrom_mask)

# Đóng file
dev.off()

print(paste("Results at", output_file)) # hiển thị kết quả
