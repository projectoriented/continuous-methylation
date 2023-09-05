if TECH == "ont":
    rule modkit:
        input:
            unpack(get_cpg_bams),
            ref=get_reference
        output:
            methyl_bed_gz="results/ont/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_{suffix}.bed.gz"
        wildcard_constraints:
            suffix="cpg-pileup|hap1_cpg-pileup|hap2_cpg-pileup"
        log:
            "results/ont/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_{suffix}.log",
        threads: 16
        resources:
            mem=lambda wildcards, attempt: attempt * 4,
            hrs=72,
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

    rule bedgraph:
        input:
            methyl_bed_gz=rules.modkit.output.methyl_bed_gz
        output:
            bedgraph=temp("results/ont/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_{suffix}.bedgraph")
        wildcard_constraints:
            suffix="cpg-pileup|hap1_cpg-pileup|hap2_cpg-pileup"
        threads: 1
        resources:
            mem=lambda wildcards, attempt: attempt * 8,
            hrs=72,
        shell:
            """
            # Columns grabbed are based on this documentation: https://github.com/nanoporetech/modkit/#bedmethyl-column-descriptions
            
            # chrom, start, n_valid, n_mod
            zcat {input.methyl_bed_gz} | awk '{{print $1,$2,$3,$11}}' FS='\\t' OFS='\\t' > {output.bedgraph}
            """

    rule bedgraph_to_bigwig:
        input:
            bedgraph=rules.bedgraph.output.bedgraph,
            chrom_sizes="results/ont/{ref}/chrom.sizes"
        output:
            bigwig="results/ont/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_{suffix}.bw"
        wildcard_constraints:
            suffix="cpg-pileup|hap1_cpg-pileup|hap2_cpg-pileup"
        threads: 1
        resources:
            mem=lambda wildcards, attempt: attempt * 8,
            hrs=72,
        envmodules:
            "modules",
            "modules-init",
            "modules-gs/prod",
            "modules-eichler/prod",
            f"ucsc/{UCSC_VERSION}",
        shell:
            """
            bedGraphToBigWig {input.bedgraph} {input.chrom_sizes} {output.bigwig}
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
        threads: 16
        resources:
            mem=lambda wildcards, attempt: attempt * 4,
            hrs=72,
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