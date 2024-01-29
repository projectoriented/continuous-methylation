if TECH == "ont":
    rule modkit:
        input:
            unpack(get_5mC_bams),
            ref=get_reference
        output:
            methyl_ungrouped_bed_gz="results/ont/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_cpg-pileup_unknown.bed.gz",
            methyl_hap1_bed_gz="results/ont/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_cpg-pileup_hap1.bed.gz",
            methyl_hap2_bed_gz="results/ont/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_cpg-pileup_hap2.bed.gz",
        params:
            search_string = "results/ont/{ref}/methylation/phased/{phase_type}/{sample}/{sample}"
        log:
            "results/ont/{ref}/methylation/phased/{phase_type}/{sample}/{sample}.log",
        threads: 10
        resources:
            mem=lambda wildcards, attempt: attempt * 2,
            hrs=72,
        envmodules:
            "modules",
            "modules-init",
            "modules-gs/prod",
            "modules-eichler/prod",
            f"modkit/{MODKIT_VERSION}"
        shell:
            """
            standard_params="--preset traditional --only-tabs --log-filepath {log} --ref {input.ref} --threads {threads}"
            
            outdir=$(dirname {output.methyl_hap1_bed_gz})
            modkit pileup {input.bam} $outdir \
                --prefix {wildcards.sample} \
                --partition-tag HP \
                $standard_params
            
            prefix={params.search_string}_cpg-pileup
            mv {params.search_string}_ungrouped.bed ${{prefix}}_unknown.bed && gzip ${{prefix}}_unknown.bed
            mv {params.search_string}_1.bed ${{prefix}}_hap1.bed && gzip ${{prefix}}_hap1.bed
            mv {params.search_string}_2.bed ${{prefix}}_hap2.bed && gzip ${{prefix}}_hap2.bed
            """

    rule modkit_combined:
        input:
            unpack(get_5mC_bams),
            ref=get_reference
        output:
            methyl_combined_bed_gz="results/ont/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_cpg-pileup_combined.bed.gz"
        log:
            "results/ont/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_combined.log",
        threads: 10
        resources:
            mem=lambda wildcards, attempt: attempt * 2,
            hrs=72,
        envmodules:
            "modules",
            "modules-init",
            "modules-gs/prod",
            "modules-eichler/prod",
            f"modkit/{MODKIT_VERSION}"
        shell:
            """
            standard_params="--preset traditional --only-tabs --log-filepath {log} --ref {input.ref} --threads {threads}"
            
            # get non-haplotyped bed
            outname=$(echo {output.methyl_combined_bed_gz} | sed 's/\.gz//')
            modkit pileup {input.bam} $outname $standard_params && gzip $outname
            """

    rule bedgraph:
        input:
            methyl_bed_gz="results/ont/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_{suffix}.bed.gz"
        output:
            bedgraph=temp("results/ont/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_{suffix}.bedgraph")
        wildcard_constraints:
            suffix="cpg-pileup_combined|cpg-pileup_hap1|cpg-pileup_hap2|cpg-pileup_unknown"
        threads: 1
        resources:
            mem=calc_mem_gb,
            hrs=72,
        shell:
            """
            # Columns grabbed are based on this documentation: https://github.com/nanoporetech/modkit/#bedmethyl-column-descriptions
            
            # Fixating the locale is necessary for bedGraphToBigWig to not freak out
            export LC_COLLATE=C
            
            # column 4 is single letter code for modified base
            # chrom, start, end, fraction
            zcat {input.methyl_bed_gz} | awk '$4 ~ /m/ {{print $1,$2,$3,$11}}' FS=' ' OFS='\\t' | sort -k 1,1 -k2,2n > {output.bedgraph}
            """

    rule bedgraph_to_bigwig:
        input:
            bedgraph=rules.bedgraph.output.bedgraph,
            chrom_sizes="results/ont/{ref}/chrom.sizes"
        output:
            bigwig="results/ont/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_{suffix}.bw"
        wildcard_constraints:
            suffix="cpg-pileup|cpg-pileup_hap1|cpg-pileup_hap2|cpg-pileup_unknown|cpg-pileup_combined"
        threads: 1
        resources:
            mem=calc_mem_gb,
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
            unpack(get_5mC_bams),
            ref=get_reference
        output:
            methyl_out = expand("results/hifi/{{ref}}/methylation/phased/{{phase_type}}/{{sample}}/{{sample}}_cpg-pileup.{out_type}.{out_ext}", out_type=["combined", "hap1", "hap2"], out_ext=["bed.gz", "bw"])
        threads: 16
        resources:
            mem=lambda wildcards, attempt: attempt * 4,
            hrs=72,
        params:
            output_prefix="results/hifi/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_cpg-pileup"
        envmodules:
            "modules",
            "modules-init",
            "modules-gs/prod",
            "modules-eichler/prod",
            f"pb-CpG-tools/{pbCpGtools_VERSION}",
        log:
            "results/hifi/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_cpg-pileup.log",
        shell:
            """
            aligned_bam_to_cpg_scores \
                --bam {input.bam} \
                --output-prefix {params.output_prefix} \
                --model $PB_MODEL/pileup_calling_model.v1.tflite \
                --threads {threads} \
                2> {log}
            
            for suffix in combined hap1 hap2
            do
                gzip {params.output_prefix}.${{suffix}}.bed
            done 
            """