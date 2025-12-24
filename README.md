# Bioinformatics-Learning-Journey

## Raw Data

Dữ liệu Raw FASTQ quá lớn để lưu trên Github

Để chạy dự án này, hãy chạy lệnh sau để tải read 1 và 2 của **SRR2584863** về máy

### Read 1 (Đầu trái)

`wget -O data/reads_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_1.fastq.gz`

### Read 2 (Đầu phải)

`wget -O data/reads_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_2.fastq.gz`

### Kiểm tra dung lượng file (Phải tầm 130MB - 140MB)
ls -lh data/

### Kiểm tra nội dung bên trong
zcat data/reads_1.fastq.gz | head -n 4

