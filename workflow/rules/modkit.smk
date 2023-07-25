rule modkit_trio:
    input:
        bam="results/{tech}/{ref}/align/phased/{phase_type}/{sample}/{sample}_{hap}_sorted-linked.bam",
        bai="results/{tech}/{ref}/align/phased/{phase_type}/{sample}/{sample}_{hap}_sorted-linked.bam.bai",
        ref=get_reference
    output:
        methyl_bed_gz="results/{tech}/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_{hap}_cpg-pileup.bed.gz"
    log:
        "results/{tech}/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_{hap}_cpg-pileup.log",
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
        methyl_bed_gz="results/{tech}/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_cpg-pileup.bed.gz"
    log:
        "results/{tech}/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_cpg-pileup.log",
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
