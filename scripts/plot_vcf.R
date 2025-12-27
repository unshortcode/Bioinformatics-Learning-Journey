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

# 1. Đọc file VCF
vcf <- read.vcfR(vcf_file, verbose = FALSE)

# 2. Lấy dữ liệu cơ bản (CHROM, POS, ID, REF, ALT, QUAL, FILTER)
# getFIX là hàm của vcfR để lấy các cột thông tin cố định
vcf_data <- getFIX(vcf)

# Chuyển nó thành bảng dữ liệu (Data Frame) quen thuộc
df <- as.data.frame(vcf_data)

# 3. Làm sạch dữ liệu
# Cột QUAL (Chất lượng) đang ở dạng chữ (Character), phải chuyển sang số (Numeric) để vẽ
df$QUAL <- as.numeric(as.character(df$QUAL))

# Kiểm tra xem có bao nhiêu dòng
print(paste("Number of Variants:", nrow(df)))

# Tạo canvas vẽ từ dữ liệu df, trục x là cột QUAL
p <- ggplot(data = df, aes(x = QUAL)) +
  # Vẽ biểu đồ cột (Histogram)
  geom_histogram(fill = "steelblue", color = "white", bins = 30) +
  # Gắn nhãn cho dễ hiểu
  labs(
    title = "Phan bo diem chat luong (Quality Score)",
    x = "Diem Chat Luong (QUAL)",
    y = "So luong bien the"
  ) +
  # Chọn giao diện nền trắng đơn giản
  theme_minimal()

print(paste("Results at", output_file))
