rule all:
    input:
        "results/qc/reads_1_fastqc.html",
        "results/qc/reads_2_fastqc.html",

        "results/trimmed/reads_1_paired.fastq.gz",

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
        # File bị loại, giữ lại để đối chứng)
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