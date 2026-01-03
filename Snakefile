rule all:
    input:
        "results/qc/reads_1_fastqc.html",
        "results/qc/reads_2_fastqc.html",
        "results/mapped/aligned.bam",
        #"results/variants/vcf_stats.txt",
        "results/variants/filtered_targeted.vcf"
        #"results/plots/filtered_variants_visualization.pdf"
        "results/plots/target_gene_visualization.pdf"

# RULE: TRIMMOMATIC
rule trim_reads:
    input:
        r1 = "data/reads_1.fastq.gz",
        r2 = "data/reads_2.fastq.gz",
        adapters = "refs/adapters.fa" # Input phụ thuộc vào file adapter
    output:
        # File quan trọng dùng để đi tiếp
        r1_paired = "results/trimmed/reads_1_paired.fastq.gz",
        r2_paired = "results/trimmed/reads_2_paired.fastq.gz",
        # File bị loại, giữ lại để đối chứng
        r1_unpaired = "results/trimmed/reads_1_unpaired.fastq.gz",
        r2_unpaired = "results/trimmed/reads_2_unpaired.fastq.gz"
    log:
        "logs/trimmomatic.log" # BẮT BUỘC PHẢI CÓ LOG
    params:
        # Cấu hình trim
        trimmer = "ILLUMINACLIP:refs/adapters.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36"
    shell:
        # Lệnh chạy (2> {log} để ghi thông báo vào file log thay vì hiện lên màn hình)
        "trimmomatic PE {input.r1} {input.r2} "
        "{output.r1_paired} {output.r1_unpaired} "
        "{output.r2_paired} {output.r2_unpaired} "
        "{params.trimmer} 2> {log}"

# RULE: Index cho reference genome
rule bwa_index:
    input:
        "refs/ecoli_ref.fasta"
    output:
        "refs/ecoli_ref.fasta.bwt"
    shell:
        "bwa index {input}"

# RULE: Alignment
rule bwa_map:
    input:
        ref = "refs/ecoli_ref.fasta",
        index = "refs/ecoli_ref.fasta.bwt", # Bắt buộc phải có index mới chạy được
        r1 = "results/trimmed/reads_1_paired.fastq.gz", # Lấy từ output của Trimmomatic
        r2 = "results/trimmed/reads_2_paired.fastq.gz"
    output:
        "results/mapped/aligned.bam"
    threads: 4  # Dùng 4 nhân CPU 
    shell:
        # Giải thích lệnh:
        # bwa mem: Thuật toán mapping
        # -t {threads}: Số luồng CPU
        # | samtools view -Sb -: Nén ngay lập tức từ SAM sang BAM (Tiết kiệm dung lượng)
        "bwa mem -t {threads} {input.ref} {input.r1} {input.r2} | samtools view -b - > {output}"

# RULE: Target Filter
rule filter_targets:
    input:
        bam = "results/mapped/aligned.bam",
        bed = "refs/targets.bed"
    output:
        "results/mapped/filtered_targets.bam"
    shell:
        # bedtools intersect chỉ giữ lại những read chồng lên với vùng liệt kê trong file .bed
        # -a là file .bam đầu vào
        # -b là file .bed chứa vị trí của target
        # -header giữ lại header của file .bam
        "bedtools intersect -a {input.bam} -b {input.bed} -header > {output}"

# RULE: VARIANT CALLING
rule call_variants:
    input:
        ref = "refs/ecoli_ref.fasta",
        # bam = "results/mapped/aligned.bam"
        bam = "results/mapped/filtered_targets.bam"
    output:
        # "results/variants/raw_variants.vcf"
        "results/variants/raw_targeted.vcf"
    log:
        "logs/bcftools.log"
    shell:
        "(samtools sort {input.bam} -o - | "  # sắp xếp các đoạn reads theo đúng thứ tự NST đi từ 1
                                              # -o output cho vào ống dẫn `-`
        "bcftools mpileup -f {input.ref} -Ou - | " # mpileup lấy các reads chồng lên ref xem vị trí khớp và lệch
                                                   # -f file ref
                                                   # -O output format là `u` uncompressed BCF để máy đọc nhanh hơn
        "bcftools call -mv -Ov -o {output}) 2>{log}"  # variant calling với -m là Multicallectic caller cho phép 1 vị trí có nhiều kiểu đột biến và `v` loại bỏ những dòng không có đột biến
                                                      # -O output format là `v` VCF 
                                                      # -o ghi kết quả vào file output
                                                      # 2>{log} chỉ ghi lại lỗi vào file log

# RULE: FILTERING
rule filter_variants:
    input:
        "results/variants/raw_targeted.vcf"
    output:
        "results/variants/filtered_targeted.vcf"
    shell:
        "bcftools filter -O v -o {output} -e 'QUAL<20 || DP<10' {input}"

# RULE: Thống kê biến thể (stats)
rule vcf_stats:
    input:
        #"results/variants/filtered_variants.vcf"
        "results/variants/filtered_targeted.vcf"
    output:
        #"results/variants/vcf_stats.txt"
        "results/variants/vcf_stats_targeted.txt"
    shell:
        # bcftools stats: Tính toán thống kê
        "bcftools stats {input} > {output}"

# RULE: Vẽ biểu đồ bằng R scripts
rule plot_quality:
    input:
        #"results/variants/filtered_variants.vcf"        
        "results/variants/filtered_targeted.vcf"        
    output:
        #"results/plots/filtered_variants_visualization.pdf"
        "results/plots/target_gene_visualization.pdf"
    shell:
        "Rscript scripts/plot_vcf.R {input} {output}"