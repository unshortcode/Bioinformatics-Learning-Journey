# Bioinformatics-Learning-Journey

## Pretiques

### Tools

### Directory

`mkdir data refs results scripts`

-   data để lưu trữ dữ liệu thô của các read
-   refs để lưu trữ Reference Genome của E. coli K12
-   results để lưu trữ các kết quả và hình ảnh
-   scripts để lưu trữ code Python, R, Bash

## Reference Genome

### Download and rename to 'ecoli_ref.fasta.gz'

`wget -O refs/ecoli_ref.fasta.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz`

### Unzip file .gz to .fasta

`gunzip refs/ecoli_ref.fasta.gz`

## Raw Data

-   Dữ liệu Raw FASTQ quá lớn để lưu trên Github
-   Để chạy dự án này, hãy chạy lệnh sau để tải read 1 và 2 của **SRR2584863** về máy

### Read 1 (Đầu trái)

`wget -O data/reads_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_1.fastq.gz`

### Read 2 (Đầu phải)

`wget -O data/reads_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_2.fastq.gz`

### Kiểm tra dung lượng file (Phải tầm 130MB - 140MB)

`ls -lh data/`

### Kiểm tra nội dung bên trong

`zcat data/reads_1.fastq.gz \| head -n 4`

## Quality Control (QC)

`fastqc data/reads_1.fastq.gz data/reads_2.fastq.gz -o results/qc/`

-   Công cụ: FastQC v0.12.1
-   Dữ liệu đầu vào (Input): SRR2584863 (Illumina MiSeq, Paired-end)
-   Kết quả read 1:
    -   Chiều dài 150 bp, 50 % GC content, phù hợp với E. Coli (50.8 %).
    -   Chất lượng cao (\>Q30) cho hầu hết các base, có phần hướng xuống nhẹ ở đầu 3'.
    -   Adapter Content: Độ sạch cao, đường biểu diễn nằm sát đáy biểu đồ, tuy nhiên có phần cong nhẹ đi lên ở cuối (\>140 bp).
-   Kết quả read 2:
    -   Chiều dài 150 bp, 50 % GC content giống với read 1.
    -   Chất lương trung bình có xu hướng giảm nhanh từ 120 bp, các đoạn cuối có chất lượng thấp ("\<"Q20).
    -   Tương tự read 1, bị nhiễu ở adapter đoạn cuối

--\> Cắt bỏ (Trimming) đoạn cuối bị nhiễu và Adapter

## Trimming

-   Tool: Trimmomatic

### Tải file TruSeq3-PE.fa (Adapter chuẩn cho Illumina Paired-End)

`wget -O refs/adapters.fa https://github.com/timflutre/trimmomatic/raw/master/adapters/TruSeq3-PE.fa`

### Cập nhật Rule Trimming cho Snakefile.

-   Kết quả nằm ở file logs/trimmomatic.log, cho thấy:
    -   Tỉ lệ Both Surviving là 87.56 % (cao), forward only surviving là 11.01 %, reverse only Surviving là 0.72 %, Drop là 0.72 % --\> cho thấy dữ liệu gốc khá sạch.

## Alignment

-   Tools: bwa

### Cập nhật Rule index và Mapping cho Snakefile

-   Chạy lệnh: snakemake -c4 (Dùng 4 core).
-   Dùng lệnh sau để kiểm tra xem file BAM có thực sự chứa dữ liệu không (Lệnh này sẽ báo cáo bao nhiêu % reads đã được map thành công).

`samtools flagstat results/mapped/aligned.bam`

Kết quả:

```         
2585326 + 0 mapped (94.40% : N/A) # bwa tìm được vị trí chính xác cho 94.4 % số lượng reads trên reference genome của E. coli K12
...
2479400 + 0 properly paired (91.16% : N/A) # điều này chứng minh read 1 và 2 đều nằm trên cùng nhiễm sắc thể và hướng vào nhau (ngược chiều)
...
19219 + 0 singletons (0.71% : N/A) # 0.71 % read 1 map được, trimming đã loại bỏ các cặp chất lượng kém
```

## Variant Calling
- Tools: samtools, bcftools

### Cập nhật Rule Variant Calling cho Snakefile

### Đếm số đột biến

`grep -v "^#" results/variants/raw_variants.vcf | wc -l`

- Giải thích đoạn code trên:
    - `^` bắt đầu dòng với ký tự `#` (trong vcf thì "\#" đại diện cho các dòng KHÔNG có dữ liệu đột biến)
    - `-v` là invert match
    --> tìm những dòng CÓ dữ liệu đột biến
- Kết quả là: 34834 khác biệt
    --> Quá nhiều đột biến, có khả năng là nhiễu hoặc chủng E. coli khác
    --> Tiến hành lọc

## Filter Variants Calling

### Cập nhật Rule FILTERING cho Snakefile
- Chạy lệnh `grep -v "^#" results/variants/filtered_variants.vcf | wc -l` để đếm số lượng đột biến đã qua filter
- Kết quả: 34044
--> Khả năng cao là khác chủng.

## VCF Stats

### Cập nhật Rule Thống kê biến thể (stats) cho Snakefile
- 
