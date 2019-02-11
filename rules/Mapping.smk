rule bismark_map:
    input:
        genome="data/genome.fa",
        sample1="data/samples/A.fq.gz"
        sample2="data/samples/B.fq.gz"
    output:
        "mapped_reads/A.bam"
    threads: 8
    shell:
        "bismark --multicore {threads} --genome {input.genome} -1 {input.sample1} -2 {input.sample2}[2] -o {output}"
