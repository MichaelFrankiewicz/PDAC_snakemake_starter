###############################################################################
# PDAC RNA-seq pipeline: QC → STAR align → BAM index → BigWig + MultiQC
###############################################################################
import csv
from pathlib import Path

SAMPLESHEET = config["samplesheet"]
REF = config["reference"]
SEQ = config["sequencing"]

# Parse samples.csv
SAMPLES = []
ROWS = {}
with open(SAMPLESHEET) as f:
    reader = csv.DictReader(f)
    for row in reader:
        s = row["sample"].strip()
        if s.startswith("#") or s == "":
            continue
        SAMPLES.append(s)
        ROWS[s] = row

def is_paired(sample):
    layout = ROWS[sample].get("layout", "PE").strip().upper()
    return layout != "SE"

def fastq_paths(sample):
    """Return (r1, r2 or None) and whether they need to be downloaded from SRA."""
    r1 = ROWS[sample].get("fastq_1", "").strip()
    r2 = ROWS[sample].get("fastq_2", "").strip()
    sra = ROWS[sample].get("sra", "").strip()

    if r1:
        return (r1, r2 if is_paired(sample) else None, False, sra)
    # else derive from SRA accession
    r1_local = f"data/fastq/{sample}_R1.fastq.gz"
    r2_local = f"data/fastq/{sample}_R2.fastq.gz" if is_paired(sample) else None
    return (r1_local, r2_local, True, sra)

def bam_path(sample):
    return f"results/bam/{sample}.bam"

def bigwig_path(sample):
    return f"results/bigwig/{sample}.bw"

def star_index_dir():
    return REF["star_index"]

rule all:
    input:
        expand("results/fastqc/{sample}_R1_fastqc.html", sample=SAMPLES),
        expand(bam_path, sample=SAMPLES),
        expand(lambda s: f"{bam_path(s)}.bai", s=SAMPLES),
        expand(bigwig_path, sample=SAMPLES),
        "results/multiqc/multiqc_report.html"

###############################################################################
# DOWNLOAD (SRA → FASTQ) if needed
###############################################################################
rule sra_to_fastq_pe:
    input:
        lambda wildcards: ROWS[wildcards.sample]["sra"]
    output:
        r1 = lambda wc: fastq_paths(wc.sample)[0],
        r2 = lambda wc: fastq_paths(wc.sample)[1],
    params:
        outdir = "data/fastq"
    threads: 4
    resources:
        mem_mb=12000, time="02:00:00"
    conda:
        "envs/sra.yaml"
    shell:
        r"""
        set -euo pipefail
        accession="{input}"
        test -n "$accession"
        mkdir -p {params.outdir}
        fasterq-dump --threads {threads} --outdir {params.outdir} "$accession"
        pigz -p {threads} {params.outdir}/{accession}_1.fastq
        pigz -p {threads} {params.outdir}/{accession}_2.fastq
        ln -sf {params.outdir}/{accession}_1.fastq.gz {output.r1}
        ln -sf {params.outdir}/{accession}_2.fastq.gz {output.r2}
        """

rule sra_to_fastq_se:
    input:
        lambda wildcards: ROWS[wildcards.sample]["sra"]
    output:
        r1 = lambda wc: fastq_paths(wc.sample)[0]
    params:
        outdir = "data/fastq"
    threads: 4
    resources:
        mem_mb=8000, time="01:00:00"
    conda:
        "envs/sra.yaml"
    shell:
        r"""
        set -euo pipefail
        accession="{input}"
        test -n "$accession"
        mkdir -p {params.outdir}
        fasterq-dump --threads {threads} --outdir {params.outdir} "$accession"
        pigz -p {threads} {params.outdir}/{accession}.fastq
        ln -sf {params.outdir}/{accession}.fastq.gz {output.r1}
        """

###############################################################################
# FASTQC on raw reads
###############################################################################
rule fastqc:
    input:
        r1 = lambda wc: fastq_paths(wc.sample)[0],
        r2 = lambda wc: fastq_paths(wc.sample)[1] if is_paired(wc.sample) else None
    output:
        htmls = lambda wc: [
            f"results/fastqc/{wc.sample}_R1_fastqc.html"
        ] + ([f"results/fastqc/{wc.sample}_R2_fastqc.html"] if is_paired(wc.sample) else [])
    threads: 2
    resources:
        mem_mb=4000, time="00:30:00"
    conda:
        "envs/fastqc.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/fastqc
        fastqc -t {threads} -o results/fastqc {input.r1} {input.r2 if input.r2 else ""}
        # rename outputs to consistent sample names
        for f in results/fastqc/*_fastqc.html; do
          base=$(basename "$f")
          if [[ "$base" != "{wildcards.sample}"* ]]; then
            # handle cases where original file names were different; we just keep them
            true
          fi
        done
        """

###############################################################################
# STAR genome index (built if missing)
###############################################################################
rule star_index:
    input:
        fasta = REF["fasta"],
        gtf   = REF["gtf"]
    output:
        directory(star_index_dir())
    threads: 8
    resources:
        mem_mb=64000, time="06:00:00"
    conda:
        "envs/star.yaml"
    params:
        sjdbOverhang = lambda wc: int(SEQ["read_length"]) - 1
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output}
        STAR --runThreadN {threads} \
             --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang {params.sjdbOverhang}
        """

###############################################################################
# Alignment (STAR) → sorted BAM
###############################################################################
rule star_align:
    input:
        idx = rules.star_index.output,
        r1  = lambda wc: fastq_paths(wc.sample)[0],
        r2  = lambda wc: fastq_paths(wc.sample)[1] if is_paired(wc.sample) else None
    output:
        bam = temp(lambda wc: f"results/bam/{wc.sample}.sortedByCoord.out.bam")
    threads: 8
    resources:
        mem_mb=64000, time="12:00:00"
    conda:
        "envs/star.yaml"
    params:
        prefix = lambda wc: f"results/bam/{wc.sample}."
    shell:
        r"""
        set -euo pipefail
        readcmd="--readFilesIn {input.r1} {input.r2}"

        if [[ "{input.r2}" == "None" ]]; then
          readcmd="--readFilesIn {input.r1}"
        fi

        zcmd=""
        if [[ "{input.r1}" == *.gz ]]; then
          zcmd="--readFilesCommand zcat"
        fi

        STAR --runThreadN {threads} \
             --genomeDir {input.idx} \
             $zcmd $readcmd \
             --outFileNamePrefix {params.prefix} \
             --outSAMtype BAM SortedByCoordinate \
             --twopassMode Basic \
             --outSAMattributes NH HI AS nM XS
        """

rule bam_symlink:
    input:
        "results/bam/{sample}.sortedByCoord.out.bam"
    output:
        bam = "results/bam/{sample}.bam"
    threads: 1
    resources:
        mem_mb=1000, time="00:05:00"
    shell:
        r"""
        ln -sf {input} {output.bam}
        """

rule bam_index:
    input:
        bam = "results/bam/{sample}.bam"
    output:
        bai = "results/bam/{sample}.bam.bai"
    threads: 2
    resources:
        mem_mb=2000, time="00:20:00"
    conda:
        "envs/samtools.yaml"
    shell:
        r"""
        samtools index -@ {threads} {input.bam}
        """

###############################################################################
# BigWig coverage (deepTools)
###############################################################################
rule bam_to_bigwig:
    input:
        bam = "results/bam/{sample}.bam",
        bai = "results/bam/{sample}.bam.bai"
    output:
        bw = "results/bigwig/{sample}.bw"
    threads: 4
    resources:
        mem_mb=16000, time="02:00:00"
    conda:
        "envs/deeptools.yaml"
    params:
        binsize = 10,
        norm = "CPM"   # choices: CPM, RPKM, BPM; adjust as needed
    shell:
        r"""
        mkdir -p results/bigwig
        bamCoverage -b {input.bam} -o {output.bw} \
          --binSize {params.binsize} \
          --normalizeUsing {params.norm} \
          --numberOfProcessors {threads}
        """

###############################################################################
# MultiQC
###############################################################################
rule multiqc:
    input:
        expand("results/fastqc/{sample}_R1_fastqc.html", sample=SAMPLES)
    output:
        report = "results/multiqc/multiqc_report.html"
    threads: 2
    resources:
        mem_mb=4000, time="00:20:00"
    conda:
        "envs/multiqc.yaml"
    shell:
        r"""
        mkdir -p results/multiqc
        multiqc -o results/multiqc .
        """

# Choose which SRA rule applies dynamically when sra is provided
use rule sra_to_fastq_pe as _sra_to_fastq_pe with:
    input: lambda wc: ROWS[wc.sample]["sra"]
    output:
        r1 = lambda wc: fastq_paths(wc.sample)[0],
        r2 = lambda wc: fastq_paths(wc.sample)[1]

use rule sra_to_fastq_se as _sra_to_fastq_se with:
    input: lambda wc: ROWS[wc.sample]["sra"]
    output:
        r1 = lambda wc: fastq_paths(wc.sample)[0]
