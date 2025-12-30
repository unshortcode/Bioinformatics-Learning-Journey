# Load thư viện cần thiết
# vcfR: Để đọc file VCF chuyên dụng
library(vcfR)
library(ggplot2)

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

# --- XỬ LÝ BIỂU ĐỒ KIỂU THAY THẾ (Substitution Type) ---
# Lấy thông tin cố định (CHROM, POS, ID, REF, ALT, QUAL, FILTER)
vcf_fix <- as.data.frame(getFIX(vcf))
vcf_fix$QUAL <- as.numeric(as.character(vcf_fix$QUAL))

# Chỉ lấy các biến thể là SNP đơn (độ dài REF và ALT đều bằng 1)
# Loại bỏ Indels (thêm/bớt)
snps <- vcf_fix[
  nchar(as.character(vcf_fix$REF)) == 1 & nchar(as.character(vcf_fix$ALT)) == 1,
]

# Tạo cột kiểu thay thế (ví dụ: "A>G")
snps$Type <- paste(snps$REF, snps$ALT, sep = ">")

# Đếm số lượng từng kiểu
type_counts <- as.data.frame(table(snps$Type))
colnames(type_counts) <- c("Type", "Count")

p <- ggplot(type_counts, aes(x = Type, y = Count)) +
  geom_bar(stat = "identity", fill = "#E69F00", color = "white") +
  theme_minimal() +
  labs(
    title = "SNP Substitution Spectrum",
    subtitle = "Distribution of observed SNP mutation classes",
    x = "Substitution (REF > ALT)",
    y = "Frequency"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- TẠO FILE PDF LƯU CÁC BIỂU ĐỒ ---
# Mở file PDF
pdf(output_file)

# In biểu đồ vào file pdf
plot(chrom_mask)
chromoqc(chrom_mask)
print(p)

# Đóng file
dev.off()

print(paste("Results at", output_file)) # Hiển thị kết quả
