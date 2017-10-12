from snakemake.io import glob_wildcards, expand
import sys
import os

directories, files, pair = glob_wildcards("data/rna_seq/SID{dir}/{file}_{P}_001.fastq.gz")
print(directories, files, pair)

rule all:
    input:
        expand("/Users/rewatitappu/Documents/SID{dir}/{file}_bam.gtf", zip, dir=directories, file=files)

rule hisat_map:
    input:
        r1 = "data/rna_seq/SID{dir}/{file}_R1_001.fastq.gz",
        r2 = "data/rna_seq/SID{dir}/{file}_R2_001.fastq.gz"
    output:
        aln = "/Users/rewatitappu/Documents/SID{dir}/{file}.bam"
    log: "hisat_out.log"
    params:
        idx = "data/indexes/chrX_tran"
    threads: 8
    shell:
        "./hisat2 -p {threads} --dta -x {params.idx} -1 {input.r1} -2 {input.r2} | ./samtools view -Sb -o {output.aln}"

rule sambamba_sort:
    input:
        ubf = "/Users/rewatitappu/Documents/SID{dir}/{file}.bam"
    output:
        sbf = "/Users/rewatitappu/Documents/SID{dir}/{file}_sorted.bam"
    params:
        mem = "4G"
    shell:
        "./sambamba_v0.6.6 sort -m {params.mem} {input.ubf} -o {output.sbf}"

rule stringtie:
    input:
        bam_file = "/Users/rewatitappu/Documents/SID{dir}/{file}_sorted.bam"
    output:
        gtf_file = "/Users/rewatitappu/Documents/SID{dir}/{file}.gtf"
    params:
        g = "data/chrX.gtf"
    shell:
        "./stringtie {input.bam_file} -G {params.g} -o {output.gtf_file}"

rule stringtie_merge:
    input:
        gtf1 = expand("/Users/rewatitappu/Documents/SID{dir}/{file}.gtf", zip, dir=directories, file=files)
    output:
        "data/merged.txt"
    run:
        All = []
        for i in input.gtf1:
            All.append(i)
        Str = " ".join(All)
        os.system("./stringtie --merge -G data/chrX.gtf -o " + output[0] + " " + Str)

rule stringtie_ballgown:
    input:
        file = "data/merged.txt",
        bamf = "/Users/rewatitappu/Documents/SID{dir}/{file}_sorted.bam"
    output:
        outf = "/Users/rewatitappu/Documents/SID{dir}/{file}_bam.gtf"
    shell:
        "./stringtie -B {input.bamf} -G {input.file} -o {output.outf}"