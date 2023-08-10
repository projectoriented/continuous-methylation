if TECH == "ont":
    rule modkit_trio:
        input:
            bam="results/ont/{ref}/align/phased/{phase_type}/{sample}/{sample}_{hap}_sorted-linked.bam",
            bai="results/ont/{ref}/align/phased/{phase_type}/{sample}/{sample}_{hap}_sorted-linked.bam.bai",
            ref=get_reference
        output:
            methyl_bed_gz="results/ont/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_{hap}_cpg-pileup.bed.gz"
        log:
            "results/ont/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_{hap}_cpg-pileup.log",
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

    rule modkit_unphased:
        input:
            unpack(get_modkit_unphased_inputs),
            ref=get_reference
        output:
            methyl_bed_gz="results/ont/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_cpg-pileup.bed.gz"
        log:
            "results/ont/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_cpg-pileup.log",
        threads: config["methylation"]["modkit"]["threads"]
        resources:
            mem=lambda wildcards, attempt: attempt * config["methylation"]["modkit"]["mem"],
            hrs=config["methylation"]["modkit"]["hrs"],
        envmodules:
            "modules",
            "modules-init",
            "modules-gs/prod",
            "modules-eichler/prod",
            f"modkit/{MODKIT_VERSION}",
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
            unpack(get_call_cpg_hifi_inputs),
            ref=get_reference
        output:
            methyl_bed_gz = "results/hifi/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_{suffix}.combined.bed.gz",
            methyl_bigwig = "results/hifi/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_{suffix}.combined.bw",
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
                --threads {threads} 2>&1 {log}
            && gzip $outname   
            """