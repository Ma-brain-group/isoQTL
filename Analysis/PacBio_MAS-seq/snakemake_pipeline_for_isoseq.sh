# experimental design
CCS_READS = {
    "Flp-d100-Meth-Fixed": "data/umass/Flp-d100-Meth-Fixed/2_B01/ccs_consensus-11/outputs/m54328U_230520_091738.hifi_reads.bam",
    "Cre-d100-Meth-Fixed": "data/umass/Cre-d100-Meth-Fixed/1_A01/ccs_consensus-11/outputs/m54328U_230518_222347.hifi_reads.bam",
    "1C-d100-Meth-Fixed":  None,
    "2C-d100-Meth-Fixed":  "data/umass/C9632_d101/2_B01/ccs_consensus-11/cromwell-executions/pb_ccs/7697386f-38eb-4d81-b82d-87f574ef1077/call-export_bam/execution/m54328U_231220_070308.hifi_reads.bam",
    "3C-d100-Meth-Fixed":  "data/umass/C8905_d101/3_C01/ccs_consensus-11/cromwell-executions/pb_ccs/f581b0cb-df05-4793-85b1-c855e8baf9a3/call-export_bam/execution/m54328U_231210_155246.hifi_reads.bam",
    "4C-d100-Meth-Fixed":  "data/umass/4C-d100-Meth-Fixed/1_A01/ccs_consensus-11/cromwell-executions/pb_ccs/5b2644e2-8e5b-449c-834f-a8aaec49f3d2/call-export_bam/execution/m54328U_231207_180149.hifi_reads.bam",
    "1N-d100-Meth-Fixed":  None,
    "2N-d100-Meth-Fixed":  "data/umass/N9540a_d101/4_D01/ccs_consensus-11/cromwell-executions/pb_ccs/7b636d28-3548-433c-9d41-47fb529af4b9/call-export_bam/execution/m54328U_231212_024858.hifi_reads.bam",
    "3N-d100-Meth-Fixed":  "data/umass/N1884a_d101/1_A01/ccs_consensus-11/cromwell-executions/pb_ccs/86bff6f2-8e2e-4167-9521-796bd2a82ba7/call-export_bam/execution/m54328U_231218_200815.hifi_reads.bam",
    "4N-d100-Meth-Fixed":  "data/umass/4N-d100-Meth-Fixed/2_B01/ccs_consensus-11/cromwell-executions/pb_ccs/80f75478-6f68-4407-af47-7e54e3216a29/call-export_bam/execution/m54328U_231209_045618.hifi_reads.bam",
    "A3277-d103-live-seq": "data/umass/1_A3277-d103-live-seq/1_A01/ccs_consensus-11/cromwell-executions/pb_ccs/9e5ea921-1556-4132-b2da-966c8ed5fc38/call-export_bam/execution/m54328U_240229_164729.hifi_reads.bam",
    "A4927-d103-live-seq": "data/umass/2_A4927-d103-live-seq/2_B01/ccs_consensus-11/cromwell-executions/pb_ccs/6e693c28-1925-4c27-92a3-8e700e8545f3/call-export_bam/execution/m54328U_240302_034228.hifi_reads.bam",
    "A7880-d106-live-seq": "data/umass/3_A7880-d106-live-seq/3_C01/ccs_consensus-11/cromwell-executions/pb_ccs/db004e71-1a6a-4d59-9e8a-ef86ea1dfe11/call-export_bam/execution/m54328U_240303_143919.hifi_reads.bam"
}
SAMPLE_GROUPS = {
    "merge1": ["Flp-d100-Meth-Fixed", "Cre-d100-Meth-Fixed"],
    "merge2": ["Flp-d100-Meth-Fixed", "Cre-d100-Meth-Fixed", "2C-d100-Meth-Fixed", "3C-d100-Meth-Fixed", "4C-d100-Meth-Fixed", "2N-d100-Meth-Fixed", "3N-d100-Meth-Fixed", "4N-d100-Meth-Fixed"],
    "merge3": ["Flp-d100-Meth-Fixed", "Cre-d100-Meth-Fixed", "1C-d100-Meth-Fixed", "2C-d100-Meth-Fixed", "3C-d100-Meth-Fixed", "4C-d100-Meth-Fixed", "1N-d100-Meth-Fixed", "2N-d100-Meth-Fixed", "3N-d100-Meth-Fixed", "4N-d100-Meth-Fixed"],
    "merge4": ["Flp-d100-Meth-Fixed", "Cre-d100-Meth-Fixed", "1C-d100-Meth-Fixed", "2C-d100-Meth-Fixed", "3C-d100-Meth-Fixed", "4C-d100-Meth-Fixed", "1N-d100-Meth-Fixed", "2N-d100-Meth-Fixed", "3N-d100-Meth-Fixed", "4N-d100-Meth-Fixed", "A3277-d103-live-seq", "A4927-d103-live-seq", "A7880-d106-live-seq"]
}
MAS_ADAPTERS = "ref/mas16_primers.fasta"
PCR_PRIMERS  = "ref/10x_primers.fasta"
TAG_DESIGN   = "T-12U-16B"
BARCODE_SET  = "ref/3M-february-2018-REVERSE-COMPLEMENTED.txt.gz"

# reference genome and annotation
REF_GENOME = "ref/GRCh38.primary_assembly.genome.fa" # 2023 and minimap2 still doesn't support gzipped fasta
REF_ANNOT  = "ref/gencode.v44.primary_assembly.annotation.gtf"

# additional references used by sqanti3/pigeon
CAGE_PEAK_PATH   = "ref/refTSS_v3.3_human_coordinate.hg38.sorted.bed"
POLYA_MOTIF_PATH = "ref/polyA.list.txt"
POLYA_ATLAS_PATH = "ref/atlas.clusters.2.0.GRCh38.96.bed" # sqanti3 only
JUNCTION_PATH    = "ref/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified2.sorted.tsv"

# paths to software
SQANTI3_DIR = "SQANTI3"
CUPCAKE_DIR = "cDNA_Cupcake"

# merged collapse -> sqanti3 qc -> sqanti3 filter -> make_seurat notebook -> sqanti3 qc again to get report

rule filtered_reads_all:
    input:
        expand(
            "results/{sample_group}/{sample}.mapped.sorted.min_cells_{min_cells}.bam",
            sample_group="merge4",
            sample=SAMPLE_GROUPS["merge4"],
            min_cells=[0, 10]
        )

rule sqanti3_a_few_min_cells:
    input:
        expand(
            "results/merge4/merge4.filtered.min_cells_{min_cells}.SQANTI3/merge4_classification.txt",
            min_cells=[3, 5, 7, 10]
        )

rule sqanti3_all:
    input:
        "results/merge4/merge4.filtered.SQANTI3/merge4_classification.txt"
#        "results/merge4/merge4.filtered.min_cells_3.SQANTI3/merge4_classification.txt"
rule map_and_qc_all:
    input:
        "results/multiqc_report.html"

rule bcstats_all:
    input:
        expand("results/{sample}/{sample}.knee.png", sample=CCS_READS.keys())

rule index_all:
    input:
#        expand("results/{sample}/{sample}.dedup.bam.bai", sample=CCS_READS.keys())
        expand("results/merge4/{sample}.mapped.sorted.newids.bam.bai", sample=SAMPLE_GROUPS["merge4"]),
        "results/merge4/merge4.mapped.sorted.bam.bai"

rule map_all:
    input:
        expand("results/{sample}/{sample}.knee.png", sample=CCS_READS.keys()),
        expand("results/{sample}/{sample}.mapped.sorted.bam", sample=CCS_READS.keys())

rule old_all:
    input:
        expand("results/{sample}/{sample}.knee.png", sample=CCS_READS.keys()),
        expand("results/{sample}/{sample}.mapped.sorted.bam", sample=CCS_READS.keys()),
        "results/multiqc_report.html",
        "results/merged.collapse.gff",
        "results/merged.collapse.sorted.filtered_lite.gff",
        "results/merged.collapse.SQANTI3/merged.collapse_RulesFilter_result_classification.txt",
        "results/merged.collapse.filtered.SQANTI3/merged.collapse.filtered_classification.txt",
#        "results/merged.collapse.sorted.filtered_lite.gff",
#        "results/pigeon_make_seurat",
        expand("results/{sample}/{sample}.non_real_cells.mapped.sorted.bam", sample=CCS_READS.keys()),
        "results/merged.non_real_cells.collapse.sorted.filtered_lite.gff",
#        "results/merged.non_real_cells.collapse.filtered.min_cells_3.SQANTI3/merged.non_real_cells.collapse.filtered.min_cells_3_classification.txt"
        "results/merged.non_real_cells.from_seurat.SQANTI3/merged.non_real_cells.from_seurat_classification.txt"

# We start by running skera to split parent reads into segmented reads (S-reads).
# https://skera.how/

rule skera:
    input:
        hifi_reads=lambda wildcards: CCS_READS[wildcards.sample],
        adapters=MAS_ADAPTERS
    output:
        multiext(
            "results/{sample}/{sample}.skera",
            ".bam",
            ".bam.pbi",
            ".found_adapters.csv.gz",
            ".ligations.csv",
            ".non_passing.bam",
            ".non_passing.bam.pbi",
            ".read_lengths.csv",
            ".summary.csv",
            ".summary.json"
        )
    threads: 8
    resources:
        mem_mb=16384,
        runtime=60
    conda:
        "envs/pb.yaml"
    shell:
        "skera split -j{threads} {input.hifi_reads} {input.adapters} {output[0]}"

# We then proceed with the standard isoseq CLI single-cell workflow, starting with lima.
# https://lima.how/
# https://isoseq.how/umi/cli-workflow.html

# Lima removes primers and outputs a variable number of bam files depending on the
# combinations of primer pairs it detects. Since we are using non-barcoded primers,
# we can predict that lima will only output a single primer pair, .5p--3p.bam.
rule lima:
    input:
        segmented_reads="results/{sample}/{sample}.skera.bam",
        primers=PCR_PRIMERS
    output:
        multiext(
            "results/{sample}/{sample}.fl",
            ".5p--3p.bam",
            ".5p--3p.bam.pbi",
            ".5p--3p.consensusreadset.xml",
            ".consensusreadset.xml",
            ".json",
            ".lima.clips",
            ".lima.counts",
            ".lima.report",
            ".lima.summary"
        )
    params:
        prefix="results/{sample}/{sample}.fl.bam"
    conda:
        "envs/pb.yaml"
    threads:
        4
    resources:
        mem_mb=8192,
        runtime=120
    shell:
        "lima -j{threads} --isoseq --dump-clips {input.segmented_reads} {input.primers} {params.prefix}"

# Clip UMIs and cell barcodes from the reads and store them in bam tags.
rule tag:
    input:
        fl_reads="results/{sample}/{sample}.fl.5p--3p.bam"
    output:
        "results/{sample}/{sample}.flt.bam",
        "results/{sample}/{sample}.flt.bam.pbi"
    params:
        design=TAG_DESIGN
    conda:
        "envs/pb.yaml"
    threads:
        4
    resources:
        mem_mb=8192,
        runtime=120
    shell:
        "isoseq tag -j{threads} --design {params.design} {input.fl_reads} {output[0]}"

# Trim polyA tails.
rule refine:
    input:
        flt_reads="results/{sample}/{sample}.flt.bam",
        primers=PCR_PRIMERS
    output:
        "results/{sample}/{sample}.fltnc.bam"
    conda:
        "envs/pb.yaml"
    threads:
        4
    resources:
        mem_mb=8192,
        runtime=120
    shell:
        "isoseq refine -j{threads} --require-polya {input.flt_reads} {input.primers} {output[0]}"

# Correct errors in cell barcodes using a whitelist of known candidate barcodes.
rule correct:
    input:
        fltnc_reads="results/{sample}/{sample}.fltnc.bam",
        barcodes=BARCODE_SET
    output:
        "results/{sample}/{sample}.corrected.bam"
    conda:
        "envs/pb.yaml"
    threads:
        4
    resources:
        mem_mb=16384,
        runtime=120
    shell:
        "isoseq correct -j{threads} --barcodes {input.barcodes} {input.fltnc_reads} {output[0]}"

# Generate a table of barcode statistics.
rule barcode_stats:
    input:
        "results/{sample}/{sample}.corrected.bam"
    output:
        "results/{sample}/{sample}.bcstats_report.tsv",
        "results/{sample}/{sample}.bcstats_report.json"
    conda:
        "envs/pb.yaml"
    threads:
        4
    resources:
        mem_mb=16384,
        runtime=60
    shell:
        "isoseq bcstats -j{threads} -o {output[0]} -R {output[1]} {input[0]}"

# Re-run bcstats with a custom percentile cutoff. Note that the BAM tags from
# isoseq correct will be wrong unless you also re-run that with --percentile,
# but it's ok because we ran groupdedup downstream with --keep-non-real-cells.
rule barcode_stats_custom:
    input:
        "results/{sample}/{sample}.corrected.bam"
    output:
        "results/{sample}/{sample}.bcstats_report.cutoff{percentile}.tsv",
        "results/{sample}/{sample}.bcstats_report.cutoff{percentile}.json"
    conda:
        "envs/pb.yaml"
    threads:
        4
    resources:
        mem_mb=16384,
        runtime=60
    shell:
        "isoseq bcstats -j{threads} --method percentile --percentile {wildcards.percentile} -o {output[0]} -R {output[1]} {input[0]}"


# Use a script from PacBio to generate the barcode rank plot (knee plot) used in cell calling.
rule plot_knees:
    input:
        "results/{sample}/{sample}.bcstats_report.tsv"
    output:
        "results/{sample}/{sample}.knee.png"
    conda:
        "datasci"
    localrule: True
    shell:
        "python3 scripts/plot_knees.py -t {input[0]} -o results/{wildcards.sample}/{wildcards.sample}"

# Sort the bam file by cell barcode for use with groupdedup.
rule sort_by_barcode:
    input:
        "results/{sample}/{sample}.corrected.bam"
    output:
        "results/{sample}/{sample}.corrected.sorted.bam"
    params:
        mem_sort=lambda w, threads, resources: int(resources.mem_mb // threads * .75)
    threads:
        4
    resources:
        mem_mb=16384,
        runtime=60
    shell:
        "samtools sort -@{threads} -m{params.mem_sort}M -t CB -o {output[0]} {input[0]}"

# Cluster reads originating from the same founder molecule based on UMI (and cell
# barcode). This eliminates PCR duplicates (but number of reads associated with a
# molecule is still available in the ic BAM tag).

# We use samtools to copy over the BAM header lines because currently groupdedup
# crudely overwrites them, losing read group information. However we must retain
# one crucial modification that groupdedup makes to the @RG line, which is to
# remove READTYPE=CCS from the DS tag and replace it with READTYPE=TRANSCRIPT.
# Otherwise pbmm2 will refuse to map any of our reads.

# Run dedup, but retain reads from non-real cells or where barcode correction
# failed. https://isoseq.how/umi/isoseq-correct.html
rule dedup:
    input:
        corrected_reads="results/{sample}/{sample}.corrected.sorted.bam"
    output:
        "results/{sample}/{sample}.dedup.bam",
        "results/{sample}/{sample}.dedup.fasta"
    conda:
        "envs/pb.yaml"
    threads:
        4
    resources:
        mem_mb=16384,
        runtime=120
    shell:
        """
        isoseq groupdedup \\
            -j{threads} \\
            --keep-non-real-cells \\
            --keep-fail-barcodes \\
            {input.corrected_reads} \\
            {output[0]}
        samtools reheader --no-PG <( \\
            samtools view -H --no-PG {input.corrected_reads} | sed '/^@RG/s/DS:[^\\t]*/DS:READTYPE=TRANSCRIPT/'; \\
            samtools view -H --no-PG {output[0]} | grep '^@PG' \\
        ) {output[0]} > {output[0]}.tmp.bam
        mv {output[0]}.tmp.bam {output[0]}
        samtools addreplacerg -@{threads} --no-PG -R $( \\
            samtools view -H --no-PG {input.corrected_reads} | \\
                grep '^@RG' | \\
                tr '\\t' '\\n' | \\
                sed -n '/^ID:/s/^ID://p' | \\
                head -n1 \\
        ) -o {output[0]}.tmp.bam {output[0]}
        mv {output[0]}.tmp.bam {output[0]}
        pbindex -j{threads} {output[0]}
        """

rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    threads: 4
    resources:
        mem_mb=8192,
        runtime=60
    shell:
        "samtools index -@{threads} {input[0]}"

# Map reads with pbmm2, a PacBio wrapper for minimap2 that enforces recommended parameters and
# also applies some post-processing (repeated matches trimming). Note: pbmm2 doesn't output XS
# or ts BAM tags that are required by StringTie for spliced alignments (and has no option to
# do so).
rule pbmm2:
    input:
        REF_GENOME,
        "results/{sample}/{sample}.dedup.bam"
    output:
        "results/{sample}/{sample}.mapped.sorted.bam"
    log:
        "results/{sample}/{sample}.mapped.sorted.bam.log"
    conda:
        "envs/pb.yaml"
    threads:
        32 # was 16
    resources:
        mem_mb=65536, # need to narrow this down; seems to need > 32G though
        runtime=480 # was 240
    shell:
        "pbmm2 align -j{threads} --preset ISOSEQ --sort --unmapped --log-level INFO {input} > {output[0]} 2> {log}"

# We now perform some QC on the mapped reads.

# These three rules prepare input for picard metrics.
localrules: gtf_to_genepred, genepred_to_refflat, gtf_to_ribo_intervals

rule gtf_to_genepred:
    input:
        REF_ANNOT
    output:
        REF_ANNOT+".genepred"
    shell:
        "ref/gtfToGenePred {input[0]} {output[0]}"

rule genepred_to_refflat:
    input:
        REF_ANNOT+".genepred"
    output:
        REF_ANNOT+".refflat"
    shell:
        "awk -vOFS='\\t' '{{print $1,$0}}' {input[0]} > {output[0]}"

rule gtf_to_ribo_intervals:
    input:
        REF_ANNOT,
        "results/"+next(iter(CCS_READS.keys()))+".mapped.sorted.bam" # for SQ lines
    output:
        REF_ANNOT+".rRNA.interval_list"
    conda:
        "datasci_py310"
    shell:
        "python3 scripts/gtf_to_ribo_intervals.py {input} {output[0]}"

# Picard CollectRnaSeqMetrics - produces nice coverage plots
# https://gatk.broadinstitute.org/hc/en-us/articles/9570471757083-CollectRnaSeqMetrics-Picard-
rule picard_metrics:
    input:
        "results/{sample}/{sample}.mapped.sorted.bam",
        REF_ANNOT+".refflat",
        REF_ANNOT+".rRNA.interval_list"
    output:
        "results/{sample}/{sample}.RNA_Metrics",
        "results/{sample}/{sample}.coverage.pdf"
    log:
        "results/{sample}/{sample}.RNA_Metrics.log"
    threads: 1
    resources:
        mem_mb=8192,
        runtime=240 # was 120
    envmodules:
        "GATK/4.3.0.0-GCCcore-11.3.0-Java-11",
        "R-/4.3.1"
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_mb}M" CollectRnaSeqMetrics \\
            -I {input[0]} \\
            -O {output[0]} \\
            -REF_FLAT {input[1]} \\
            -STRAND FIRST_READ_TRANSCRIPTION_STRAND \\
            -RIBOSOMAL_INTERVALS {input[2]} \\
            -CHART_OUTPUT {output[1]} > {log} 2>&1
        """

rule samtools_stats:
    input:
        "results/{sample}/{sample}.mapped.sorted.bam",
        REF_GENOME
    output:
        "results/{sample}/{sample}.mapped.sorted.bam.bamstat"
    threads: 2
    resources:
        mem_mb=8192,
        runtime=120 # was 60
    shell:
        "samtools stats -@{threads} --ref-seq {input[1]} {input[0]} > {output[0]}"

# for N50 - multiQC does not actually understand its report
rule bbmap_stats:
    input:
        "results/{sample}/{sample}.mapped.sorted.bam"
    output:
        "results/{sample}/{sample}.mapped.sorted.bam.bbstats"
    threads: 2
    resources:
        mem_mb=8192,
        runtime=120 # was 60
    shell:
        "samtools fasta -@{threads} {input[0]} | ref/bbmap/stats.sh in=stdin > {output[0]}"

# Automatically detects QC outputs that it recognizes and generates a nice HTML
# report. We have to exclude lima outputs which it currently doesn't interpret
# correctly.
rule multiqc:
    input:
        dirs_to_scan=expand("results/{sample}", sample=CCS_READS.keys()),
        stats_files=expand("results/{sample}/{sample}.{ext}", sample=CCS_READS.keys(), ext=[
            "RNA_Metrics",
            "mapped.sorted.bam.bamstat",
            "mapped.sorted.bam.bbstats"
        ])
    output:
        "results/multiqc_report.html"
    conda:
        "envs/multiqc.yaml"
    threads: 1
    resources:
        mem_mb=8192,
        runtime=60
    shell:
        "multiqc -o results -f -x '*.lima*' {input.dirs_to_scan}"

# Moving on from QC of mapped reads, we now call unique isoforms from the reads.

# isoseq collapse - collapses redundant transcripts (based on exonic structure)
# into unique isoforms.
# https://isoseq.how/classification/isoseq-collapse.html

# isoseq collapse doesn't really seem to support multiple samples. however, at
# this stage we still have BAM files per sample and we want to create a single
# shared transcriptome (where isoforms shared between samples have the same ID).
# we will merge the samples, but i want "molecule/" IDs and cell barcodes to
# remain unique between samples so that the abundance.txt and group.txt produced
# by isoseq collapse remain somewhat interpretable. we can do this by suffixing
# a _0, _1, etc. to the IDs in each BAM file based on the order of the samples
# in the dict at the top of this Snakefile.
rule prepare_for_collapse:
    input:
        "results/{sample}/{sample}.mapped.sorted.bam"
    params:
        num=lambda w: SAMPLE_GROUPS[w.sample_group].index(w.sample)
    output:
        "results/{sample_group}/{sample}.mapped.sorted.newids.bam"
    threads: 4
    resources:
        mem_mb=8192,
        runtime=60
    shell:
        """
        samtools view -@1 --no-PG -h {input[0]} | \\
            sed -r \\
                -e '/^[^@]/s/(^[^\\t]*)/\\1_{params.num}/' \\
                -e '/^[^@]/s/(CB:Z:[^\\t]*)/\\1_{params.num}/' | \\
            samtools view -@1 --no-PG -bo {output[0]}
        """

# merge the samples
rule merge_for_collapse:
    input:
        lambda w: expand("results/{{sample_group}}/{sample}.mapped.sorted.newids.bam", sample=SAMPLE_GROUPS[w.sample_group])
    output:
        "results/{sample_group}/{sample_group}.mapped.sorted.bam"
    threads: 4
    resources:
        mem_mb=8192,
        runtime=120 # was 60
    shell:
        "samtools merge -@{threads} -o {output[0]} {input}"

# run collapse

# set yet another flag, this time for the collapse command, to retain non-real
# cells. note also the --min-aln-coverage and --min-aln-identity thresholds that
# may present all reads from ending up in the .read_stat.txt.
rule collapse:
    input:
        "results/{sample}/{sample}.mapped.sorted.bam" # can be sample or sample_group, aren't i clever
    output:
        "results/{sample}/{sample}.collapse.gff",
        "results/{sample}/{sample}.collapse.abundance.txt",
        "results/{sample}/{sample}.collapse.fasta",
        "results/{sample}/{sample}.collapse.flnc_count.txt",
        "results/{sample}/{sample}.collapse.group.txt",
        "results/{sample}/{sample}.collapse.read_stat.txt",
        "results/{sample}/{sample}.collapse.report.json"
    log:
        "results/{sample}/{sample}.collapse.gff.log"
    conda:
        "envs/pb.yaml"
    threads:
        32 # was 16
    resources:
        mem_mb=65536, # was 32768
        runtime=480 # was 240
    shell:
        "isoseq collapse -j{threads} --log-level INFO --keep-non-real-cells {input[0]} {output[0]} 2> {log}"

# run sqanti3 qc (instead of pigeon)
rule sqanti3:
    input:
        "results/{sample}/{sample}.collapse.gff",
        REF_ANNOT,
        REF_GENOME,
        CAGE_PEAK_PATH,
        POLYA_MOTIF_PATH,
        POLYA_ATLAS_PATH,
        JUNCTION_PATH
    output:
        "results/{sample}/{sample}.SQANTI3/{sample}_classification.txt",
        "results/{sample}/{sample}.SQANTI3/{sample}_corrected.fasta",
        "results/{sample}/{sample}.SQANTI3/{sample}_corrected.gtf"
#        "results/{sample}/{sample}.SQANTI3/{sample}_corrected.faa" # no longer generated bc --skipORF
    log:
        "results/{sample}/{sample}.SQANTI3/{sample}_SQANTI3_log.txt"
    params:
        prefix="results/{sample}/{sample}.SQANTI3"
    conda:
        "envs/sqanti3.yaml"
    threads:
        8 # only uses 1 core with skipORF
    resources:
        mem_mb=262144, # used 175G for merge3 and 201GB for merge4
        runtime=1440 # merge3 and merge4 took 12hr; could be much longer with ORF prediction?
    shell:
        """
        export PYTHONPATH="${{PYTHONPATH:+$PYTHONPATH:}}{CUPCAKE_DIR}/sequence:{CUPCAKE_DIR}"
        python -u {SQANTI3_DIR}/sqanti3_qc.py \\
            {input[0]} \\
            {input[1]} \\
            {input[2]} \\
            --CAGE_peak {input[3]} \\
            --polyA_motif_list {input[4]} \\
            --polyA_peak {input[5]} \\
            -c {input[6]} \\
            --skipORF \\
            -o {wildcards.sample} \\
            -d {params.prefix} \\
            --cpus {threads} \\
            --report both > {log} 2>&1
        """

# run sqanti3 filter - in latest versions it gives a filter report PDF specific
# to this step, but doesn't re-generate the full QC report.
# gives us about 900,000 isoforms.
rule sqanti3_filter:
    input:
        "results/{sample}/{sample}.SQANTI3/{sample}_classification.txt",
        "results/{sample}/{sample}.SQANTI3/{sample}_corrected.fasta",
        "results/{sample}/{sample}.SQANTI3/{sample}_corrected.gtf"
#        "results/{sample}/{sample}.SQANTI3/{sample}_corrected.faa"
    output:
        "results/{sample}/{sample}.SQANTI3/{sample}_RulesFilter_result_classification.txt",
        "results/{sample}/{sample}.SQANTI3/{sample}.filtered.gtf"
    params:
        prefix="results/{sample}/{sample}.SQANTI3"
    conda:
        "envs/sqanti3.yaml"
    threads:
        1
    resources:
        mem_mb=65536,
        runtime=960 # was 480
    shell:
        """
        python -u {SQANTI3_DIR}/sqanti3_filter.py rules \\
            {input[0]} \\
            --isoforms {input[1]} \\
            --gtf {input[2]} \\
            -o {wildcards.sample} \\
            -d {params.prefix}
        """
#            --faa {input[3]} \\

# since latest versions of sqanti3 are unable to re-generate the full
# report as part of the filter step
rule sqanti3_on_filtered:
    input:
        "results/{sample}/{sample}.SQANTI3/{sample}.filtered.gtf",
        REF_ANNOT,
        REF_GENOME,
        CAGE_PEAK_PATH,
        POLYA_MOTIF_PATH,
        POLYA_ATLAS_PATH,
        JUNCTION_PATH
    output:
        "results/{sample}/{sample}.filtered.SQANTI3/{sample}_classification.txt",
        "results/{sample}/{sample}.filtered.SQANTI3/{sample}_corrected.fasta",
        "results/{sample}/{sample}.filtered.SQANTI3/{sample}_corrected.gtf"
#        "results/{sample}/{sample}.filtered.SQANTI3/{sample}_corrected.faa" # no longer generated bc --skipORF
    log:
        "results/{sample}/{sample}.filtered.SQANTI3/{sample}_SQANTI3_log.txt"
    params:
        prefix="results/{sample}/{sample}.filtered.SQANTI3"
    conda:
        "envs/sqanti3.yaml"
    threads:
        8
    resources:
        mem_mb=204800, # was 163840
        runtime=2160 # was 1440
    shell:
        """
        export PYTHONPATH="${{PYTHONPATH:+$PYTHONPATH:}}{CUPCAKE_DIR}/sequence:{CUPCAKE_DIR}"
        python -u {SQANTI3_DIR}/sqanti3_qc.py \\
            {input[0]} \\
            {input[1]} \\
            {input[2]} \\
            --CAGE_peak {input[3]} \\
            --polyA_motif_list {input[4]} \\
            --polyA_peak {input[5]} \\
            -c {input[6]} \\
            --skipORF \\
            -o {wildcards.sample} \\
            -d {params.prefix} \\
            --cpus {threads} \\
            --report both > {log} 2>&1
        """

# run sqanti3 qc on isoforms filtered by presence in a minimum number of cells
# just as seurat does (and with an fl_counts csv for some additional plots).
# see make_seurat.ipynb
# e.g. {min_cells} = min_cells_3
rule sqanti3_on_min_cells:
    input:
        "results/{sample}/{sample}.filtered.{min_cells}.gff",
        REF_ANNOT,
        REF_GENOME,
        CAGE_PEAK_PATH,
        POLYA_MOTIF_PATH,
        POLYA_ATLAS_PATH,
        JUNCTION_PATH
#        "results/{sample}/{sample}.filtered.{min_cells}.fl_counts.csv"
    output:
        "results/{sample}/{sample}.filtered.{min_cells}.SQANTI3/{sample}_classification.txt",
        "results/{sample}/{sample}.filtered.{min_cells}.SQANTI3/{sample}_corrected.fasta",
        "results/{sample}/{sample}.filtered.{min_cells}.SQANTI3/{sample}_corrected.gtf",
        "results/{sample}/{sample}.filtered.{min_cells}.SQANTI3/{sample}_corrected.faa" # no longer generated if --skipORF
    log:
        "results/{sample}/{sample}.filtered.{min_cells}.SQANTI3/{sample}_SQANTI3_log.txt"
    params:
        prefix="results/{sample}/{sample}.filtered.{min_cells}.SQANTI3"
    conda:
        "envs/sqanti3.yaml"
    threads:
        8
    resources:
        mem_mb=204800, # was 163840
        runtime=1440 # was 1440 - only took 65 without ORF prediction
    shell:
        """
        export PYTHONPATH="${{PYTHONPATH:+$PYTHONPATH:}}{CUPCAKE_DIR}/sequence:{CUPCAKE_DIR}"
        python -u {SQANTI3_DIR}/sqanti3_qc.py \\
            {input[0]} \\
            {input[1]} \\
            {input[2]} \\
            --CAGE_peak {input[3]} \\
            --polyA_motif_list {input[4]} \\
            --polyA_peak {input[5]} \\
            -c {input[6]} \\
            -o {wildcards.sample} \\
            -d {params.prefix} \\
            --cpus {threads} \\
            --report both > {log} 2>&1
        """
#             -fl {input[7]} \\
#             --skipORF \\


rule filter_reads:
    input:
        bam="results/{sample}/{sample}.mapped.sorted.bam",
        collapse="results/{sample_group}/{sample_group}.collapse.group.txt",
        pbids="results/{sample_group}/{sample_group}.filtered.min_cells_{min_cells}.pbids.txt"
    output:
        "results/{sample_group}/{sample}.mapped.sorted.min_cells_{min_cells}.bam",
        "results/{sample_group}/{sample}.mapped.sorted.min_cells_{min_cells}.bam.bai"
    params:
        num=lambda w: SAMPLE_GROUPS[w.sample_group].index(w.sample)
    threads: 2
    resources:
        mem_mb=8192,
        runtime=60
    shell:
        """
        samtools view \\
            --write-index \\
            -@{threads} \\
            -N <(\\
                grep -Ff <(sed -e 's/"//g' -e 's/$/\\t/' {input.pbids}) {input.collapse} |
                cut -f2 |
                tr ',' '\\n' |
                grep '_{params.num}$' |
                sed -e 's/_{params.num}$//'\\
            ) \\
            -o {output[0]}##idx##{output[1]} \\
            {input.bam}
        """


rule ensg_of_interest:
    input:
        REF_ANNOT,
        "ref/genes_of_interest.txt"
    output:
        "ref/ensg_of_interest.txt"
    localrule: True
    shell:
        """
        awk '$3 == "gene"' {input[0]} | grep -Ff {input[1]} | sed -r 's/.*gene_id "([^."]*)\\.?.*";.*gene_name "([^"]*)";.*/\\2\\t\\1/' > {output[0]}
        """

rule pbids_of_interest:
    input:
        "results/{sample_group}/{sample_group}.SQANTI3/{sample_group}_RulesFilter_result_classification.txt",
        "ref/ensg_of_interest.txt"
    output:
        "results/{sample_group}/pbids_of_interest.txt"
    localrule: True
    shell:
        """
        awk '$49 == "Isoform"' {input[0]} | grep -Ff <(cut -f2 {input[1]}) | cut -f1 > {output[0]}
        """

SAMPLE_STATUS = {
    "merge4": {
        "Flp-d100-Meth-Fixed": "CTL",
        "Cre-d100-Meth-Fixed": "SCZ",
        "1C-d100-Meth-Fixed":  "CTL",
        "2C-d100-Meth-Fixed":  "CTL",
        "3C-d100-Meth-Fixed":  "CTL",
        "4C-d100-Meth-Fixed":  "CTL",
        "1N-d100-Meth-Fixed":  "SCZ",
        "2N-d100-Meth-Fixed":  "SCZ",
        "3N-d100-Meth-Fixed":  "SCZ",
        "4N-d100-Meth-Fixed":  "SCZ",
        "A3277-d103-live-seq": "ASD",
        "A4927-d103-live-seq": "ASD",
        "A7880-d106-live-seq": "ASD"
    }
}

rule small_bams:
    input:
        [f"results/merge4/small_bams/{SAMPLE_STATUS['merge4'][sample]}/{sample}.min_cells_{min_cells}.genes_of_interest.bam" for sample in SAMPLE_GROUPS["merge4"] for min_cells in [0, 10]]

rule filter_reads_by_gene:
    input:
        bam="results/{sample_group}/{sample}.mapped.sorted.min_cells_{min_cells}.bam",
        collapse="results/{sample_group}/{sample_group}.collapse.group.txt",
        pbids="results/{sample_group}/pbids_of_interest.txt"
    output:
        "results/{sample_group}/small_bams/{sample_status}/{sample}.min_cells_{min_cells}.genes_of_interest.bam",
        "results/{sample_group}/small_bams/{sample_status}/{sample}.min_cells_{min_cells}.genes_of_interest.bam.bai"
    params:
        num=lambda w: SAMPLE_GROUPS[w.sample_group].index(w.sample)
    threads: 2
    resources:
        mem_mb=8192,
        runtime=60
    shell:
        """
        samtools view \\
            --write-index \\
            -@{threads} \\
            -N <(\\
                grep -Ff <(sed -e 's/$/\\t/' {input.pbids}) {input.collapse} |
                cut -f2 |
                tr ',' '\\n' |
                grep '_{params.num}$' |
                sed -e 's/_{params.num}$//'\\
            ) \\
            -o {output[0]}##idx##{output[1]} \\
            {input.bam}
        """
