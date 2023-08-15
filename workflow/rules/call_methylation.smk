if TECH == "ont":
    rule modkit:
        input:
            unpack(get_cpg_bams),
            ref=get_reference
        output:
            methyl_bed_gz="results/ont/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_{suffix}.bed.gz"
        log:
            "results/ont/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_{suffix}.log",
        threads: config["methylation"]["modkit"]["threads"]
        resources:
            mem=lambda wildcards, attempt: attempt * config["methylation"]["modkit"]["mem"],
            hrs=config["methylation"]["modkit"]["hrs"],
        envmodules:
            "modules",
            "modules-init",
            "modules-gs/prod",
            "modules-eichler/prod",
            f"modkit/{MODKIT_VERSION}"
        shell:
            """
            outname=$(echo {output.methyl_bed_gz} | sed 's/\.gz//')
            modkit pileup \
                --ref {input.ref} \
                --preset traditional \
                --only-tabs \
                --log-filepath {log} \
                --threads {threads} \
                {input.bam} \
                $outname && gzip $outname
            """
elif TECH == "hifi":
    rule call_cpg_hifi:
        input:
            unpack(get_cpg_bams),
            ref=get_reference
        output:
            methyl_bed_gz = "results/hifi/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_{suffix}.combined.bed.gz",
            methyl_bigwig = "results/hifi/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_{suffix}.combined.bw",
        wildcard_constraints:
            suffix="cpg-pileup|hap1_cpg-pileup|hap2_cpg-pileup"
        threads: config["methylation"]["pb-CpG-tools"]["threads"]
        resources:
            mem=lambda wildcards, attempt: attempt * config["methylation"]["pb-CpG-tools"]["mem"],
            hrs=config["methylation"]["pb-CpG-tools"]["hrs"],
        params:
            output_prefix="results/hifi/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_{suffix}"
        envmodules:
            "modules",
            "modules-init",
            "modules-gs/prod",
            "modules-eichler/prod",
            f"pb-CpG-tools/{pbCpGtools_VERSION}",
        log:
            "results/hifi/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_{suffix}.log",
        shell:
            """
            outname=$(echo {output.methyl_bed_gz} | sed 's/\.gz//')
            aligned_bam_to_cpg_scores \
                --bam {input.bam} \
                --output-prefix {params.output_prefix} \
                --model $PB_MODEL/pileup_calling_model.v1.tflite \
                --threads {threads} \
                2> {log} && gzip $outname
            """