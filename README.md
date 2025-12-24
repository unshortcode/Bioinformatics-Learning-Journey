# Bioinformatics-Learning-Journey

## Reference Genome

### Download and rename to 'ecoli_ref.fasta.gz'

`wget -O refs/ecoli_ref.fasta.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz`

### Unzip file .gz to .fasta

`gunzip refs/ecoli_ref.fasta.gz`

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

