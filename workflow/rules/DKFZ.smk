rule count_matrix:
    input:
        expand(
            join(
                BASE_DATA_DIR,
                "{unit.sample_unit}",
                "{unit.tissue}",
                SEQ_TYPE,
                "merged-alignment/featureCounts/{unit.tissue}_{unit.sample_unit}.fpkm_tpm.featureCounts.tsv",
            ),
            unit=units.itertuples(),
        ),
    output:
        join(BASE_ANALYSIS_DIR, "counts/all.tsv"),
        join(BASE_ANALYSIS_DIR, "fpkm/true_fpkm.tsv"),
    params:
        samples=units["sample"].tolist(),
        strand=get_strandness(units),
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/count_matrix.log",
    script:
        "../scripts/count-matrix.py"


## RSEQC
rule rseqc_gtf2bed:
    input:
        config["ref"]["annotation"],
    output:
        bed=BASE_ANALYSIS_DIR + "qc/rseqc/annotation.bed",
        db=temp(BASE_ANALYSIS_DIR + "qc/rseqc/annotation.db"),
    log:
        BASE_ANALYSIS_DIR + "logs/rseqc_gtf2bed.log",
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/gtf2bed.py"


rule rseqc_junction_annotation:
    input:
        bam=join(
            BASE_DATA_DIR,
            "{sample}",
            "{tissue}",
            SEQ_TYPE,
            "merged-alignment",
            "{tissue}_{sample}_merged.mdup.bam",
        ),
        bed=BASE_ANALYSIS_DIR + "qc/rseqc/annotation.bed",
    output:
        BASE_ANALYSIS_DIR + "qc/rseqc/{sample}@{tissue}.junctionanno.junction.bed",
    priority: 1
    log:
        BASE_ANALYSIS_DIR + "logs/rseqc/rseqc_junction_annotation/{sample}@{tissue}.log",
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix=BASE_ANALYSIS_DIR + "qc/rseqc/{sample}@{tissue}.junctionanno",
    resources:
        mem_mb=4096,
        time_min=59,
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log[0]} 2>&1"


rule rseqc_junction_saturation:
    input:
        bam=join(
            BASE_DATA_DIR,
            "{sample}",
            "{tissue}",
            SEQ_TYPE,
            "merged-alignment",
            "{tissue}_{sample}_merged.mdup.bam",
        ),
        bed=BASE_ANALYSIS_DIR + "qc/rseqc/annotation.bed",
    output:
        BASE_ANALYSIS_DIR
        + "qc/rseqc/{sample}@{tissue}.junctionsat.junctionSaturation_plot.pdf",
    priority: 1
    log:
        BASE_ANALYSIS_DIR + "logs/rseqc/rseqc_junction_saturation/{sample}@{tissue}.log",
    params:
        extra="-q 255",
        prefix=BASE_ANALYSIS_DIR + "qc/rseqc/{sample}@{tissue}.junctionsat",
    resources:
        mem_mb=4096,
        time_min=59,
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log} 2>&1"


rule rseqc_stat:
    input:
        join(
            BASE_DATA_DIR,
            "{sample}",
            "{tissue}",
            SEQ_TYPE,
            "merged-alignment",
            "{tissue}_{sample}_merged.mdup.bam",
        ),
    output:
        BASE_ANALYSIS_DIR + "qc/rseqc/{sample}@{tissue}.stats.txt",
    priority: 1
    log:
        BASE_ANALYSIS_DIR + "logs/rseqc/rseqc_stat/{sample}@{tissue}.log",
    conda:
        "../envs/rseqc.yaml"
    resources:
        mem_mb=8192,
        time_min=59,
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"


rule rseqc_infer:
    input:
        bam=join(
            BASE_DATA_DIR,
            "{sample}",
            "{tissue}",
            SEQ_TYPE,
            "merged-alignment",
            "{tissue}_{sample}_merged.mdup.bam",
        ),
        bed=BASE_ANALYSIS_DIR + "qc/rseqc/annotation.bed",
    output:
        BASE_ANALYSIS_DIR + "qc/rseqc/{sample}@{tissue}.infer_experiment.txt",
    priority: 1
    log:
        BASE_ANALYSIS_DIR + "logs/rseqc/rseqc_infer/{sample}@{tissue}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_innerdis:
    input:
        bam=join(
            BASE_DATA_DIR,
            "{sample}",
            "{tissue}",
            SEQ_TYPE,
            "merged-alignment",
            "{tissue}_{sample}_merged.mdup.bam",
        ),
        bed=BASE_ANALYSIS_DIR + "qc/rseqc/annotation.bed",
    output:
        BASE_ANALYSIS_DIR
        + "qc/rseqc/{sample}@{tissue}.inner_distance_freq.inner_distance.txt",
    priority: 1
    log:
        BASE_ANALYSIS_DIR + "logs/rseqc/rseqc_innerdis/{sample}@{tissue}.log",
    params:
        prefix=BASE_ANALYSIS_DIR + "qc/rseqc/{sample}@{tissue}.inner_distance_freq",
    resources:
        mem_mb=4096,
        time_min=59,
    conda:
        "../envs/rseqc.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readdis:
    input:
        bam=join(
            BASE_DATA_DIR,
            "{sample}",
            "{tissue}",
            SEQ_TYPE,
            "merged-alignment",
            "{tissue}_{sample}_merged.mdup.bam",
        ),
        bed=BASE_ANALYSIS_DIR + "qc/rseqc/annotation.bed",
    output:
        BASE_ANALYSIS_DIR + "qc/rseqc/{sample}@{tissue}.readdistribution.txt",
    priority: 1
    log:
        BASE_ANALYSIS_DIR + "logs/rseqc/rseqc_readdis/{sample}@{tissue}.log",
    resources:
        mem_mb=4096,
        time_min=59,
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_readdup:
    input:
        join(
            BASE_DATA_DIR,
            "{sample}",
            "{tissue}",
            SEQ_TYPE,
            "merged-alignment",
            "{tissue}_{sample}_merged.mdup.bam",
        ),
    output:
        BASE_ANALYSIS_DIR + "qc/rseqc/{sample}@{tissue}.readdup.DupRate_plot.pdf",
    priority: 1
    log:
        BASE_ANALYSIS_DIR + "logs/rseqc/rseqc_readdup/{sample}@{tissue}.log",
    params:
        prefix=BASE_ANALYSIS_DIR + "qc/rseqc/{sample}@{tissue}.readdup",
    conda:
        "../envs/rseqc.yaml"
    resources:
        mem_mb=24000,
        time_min=59,
    threads: 1
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"


rule rseqc_readgc:
    input:
        join(
            BASE_DATA_DIR,
            "{sample}",
            "{tissue}",
            SEQ_TYPE,
            "merged-alignment",
            "{tissue}_{sample}_merged.mdup.bam",
        ),
    output:
        BASE_ANALYSIS_DIR + "qc/rseqc/{sample}@{tissue}.readgc.GC_plot.pdf",
    priority: 1
    log:
        BASE_ANALYSIS_DIR + "logs/rseqc/rseqc_readgc/{sample}@{tissue}.log",
    params:
        prefix=BASE_ANALYSIS_DIR + "qc/rseqc/{sample}@{tissue}.readgc",
    conda:
        "../envs/rseqc.yaml"
    resources:
        mem_mb=18000,
        time_min=59,
    threads: 1
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"


# rule fastqc:
#    input:
#        def_fq_yeet
#    output:
#        html=join(BASE_ANALYSIS_DIR, "qc/fastqc/{sample}@{tissue}.html"),
#        zip=join(BASE_ANALYSIS_DIR, "qc/fastqc/{sample}@{tissue}_fastqc.zip") # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
#   params: "--quiet"
#   log:
#       "logs/fastqc/{sample}_{tissue}.log"
#   threads: 1
#   wrapper:
#       "v1.0.0/bio/fastqc"
rule multiqc:
    input:
        expand(
            join(
                BASE_DATA_DIR,
                "{unit.sample_unit}",
                "{unit.tissue}",
                SEQ_TYPE,
                "merged-alignment",
                "{unit.tissue}_{unit.sample_unit}_merged.mdup.bam",
            ),
            unit=units.itertuples(),
        ),
        expand(
            BASE_ANALYSIS_DIR
            + "qc/rseqc/{unit.sample_unit}@{unit.tissue}.junctionanno.junction.bed",
            unit=units.itertuples(),
        ),
        expand(
            BASE_ANALYSIS_DIR
            + "qc/rseqc/{unit.sample_unit}@{unit.tissue}.junctionsat.junctionSaturation_plot.pdf",
            unit=units.itertuples(),
        ),
        expand(
            BASE_ANALYSIS_DIR
            + "qc/rseqc/{unit.sample_unit}@{unit.tissue}.infer_experiment.txt",
            unit=units.itertuples(),
        ),
        expand(
            BASE_ANALYSIS_DIR + "qc/rseqc/{unit.sample_unit}@{unit.tissue}.stats.txt",
            unit=units.itertuples(),
        ),
        expand(
            BASE_ANALYSIS_DIR
            + "qc/rseqc/{unit.sample_unit}@{unit.tissue}.inner_distance_freq.inner_distance.txt",
            unit=units.itertuples(),
        ),
        expand(
            BASE_ANALYSIS_DIR
            + "qc/rseqc/{unit.sample_unit}@{unit.tissue}.readdistribution.txt",
            unit=units.itertuples(),
        ),
        expand(
            BASE_ANALYSIS_DIR
            + "qc/rseqc/{unit.sample_unit}@{unit.tissue}.readdup.DupRate_plot.pdf",
            unit=units.itertuples(),
        ),
        expand(
            BASE_ANALYSIS_DIR
            + "qc/rseqc/{unit.sample_unit}@{unit.tissue}.readgc.GC_plot.pdf",
            unit=units.itertuples(),
        ),
        expand(
            BASE_ANALYSIS_DIR
            + "logs/rseqc/rseqc_junction_annotation/{unit.sample_unit}@{unit.tissue}.log",
            unit=units.itertuples(),
        ),
    # expand(join(BASE_ANALYSIS_DIR, "qc/fastqc/{unit.sample}@{unit.tissue}_fastqc.zip"), unit = units.itertuples())
    output:
        BASE_ANALYSIS_DIR + "/qc/multiqc_report.html",
        directory(BASE_ANALYSIS_DIR + "/qc/multiqc_data"),
    log:
        BASE_ANALYSIS_DIR + "/logs/multiqc.log",
    params:
        extra="--data-dir",
    resources:
        mem_mb=8192 * 2,
        time_min=59 * 4,
    wrapper:
        "v7.2.0/bio/multiqc"


rule export_xlsx:
    input:
        fpkm=join(BASE_ANALYSIS_DIR, "fpkm/true_fpkm.tsv"),
    output:
        xlsx=join(BASE_ANALYSIS_DIR, "fpkm/all.xlsx"),
    conda:
        "../envs/pandas.yaml"
    threads: 1
    resources:
        mem_mb=8192,
        time_min=30,
    log:
        BASE_ANALYSIS_DIR + "logs/export_xlsx.log",
    script:
        "../scripts/export_fpkm.py"
