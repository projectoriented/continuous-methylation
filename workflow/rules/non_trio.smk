rule align:
    input:
        cell_fastq=get_tech_files(which_one="fofn"),
        ref=get_reference,
    output:
        cell_bam=temp(
            "results/{tech}/{ref}/align/phased/{phase_type}/minimap2/{sample}/{sample}_{cell}_sorted.bam"
        ),
        cell_bam_bai=temp(
            "results/{tech}/{ref}/align/phased/{phase_type}/minimap2/{sample}/{sample}_{cell}_sorted.bam.bai"
        ),
    params:
        tech_arg = lambda wildcards: f"map-{wildcards.tech}"
    threads: config["align"]["minimap2"]["threads"]
    resources:
        mem=lambda wildcards, attempt: attempt * config["align"]["minimap2"]["mem"],
        hrs=config["align"]["minimap2"]["hrs"],
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"minimap2/{MINIMAP2_VERSION}",
        f"sambamba/{SAMBAMBA_VERSION}",
    shell:
        """
        minimap2 -t {threads} --MD --eqx -ax {params.tech_arg} {input.ref} {input.cell_fastq} | sambamba view --sam-input --format bam /dev/stdin | sambamba sort --nthreads {threads} --out {output.cell_bam} /dev/stdin && sambamba index --nthreads {threads} {output.cell_bam}
        """

rule link_meth_tags_non_trio:
    input:
        cell_bam=rules.align.output.cell_bam,
        cell_bam_bai=rules.align.output.cell_bam_bai,
        unmapped_bam=get_tech_files(which_one="unmapped_bam_fofn"),
    output:
        cell_linked_bam=temp("results/{tech}/{ref}/align/phased/{phase_type}/minimap2/{sample}/{sample}_{cell}_sorted-linked.bam"),
        cell_linked_bam_bai=temp("results/{tech}/{ref}/align/phased/{phase_type}/minimap2/{sample}/{sample}_{cell}_sorted-linked.bam.bai")
    threads: config["methylation"]["methylink"]["threads"]
    resources:
        mem=lambda wildcards, attempt: attempt * config["methylation"]["methylink"]["mem"],
        hrs=config["methylation"]["methylink"]["hrs"],
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"methylink/{METHYLINK_VERSION}",
    shell:
        """
        methylink \
          --threads {threads} \
          --aln {input.cell_bam} \
          --sample {wildcards.sample}_{wildcards.cell} \
          --methyl_bams "$(ls {input.unmapped_bam})" \
          --output {output.cell_linked_bam}
        """

rule merge_align:
    input:
        cell_bams=gather_tech_bams(which_one="not_hap_specific"),
    output:
        merged_bam=temp("results/{tech}/{ref}/align/phased/{phase_type}/minimap2/{sample}/{sample}_sorted-linked.bam"),
        merged_bam_bai=temp("results/{tech}/{ref}/align/phased/{phase_type}/minimap2/{sample}/{sample}_sorted-linked.bam.bai"),
    threads: config["align"]["sambamba"]["threads"]
    resources:
        mem=lambda wildcards, attempt: attempt * config["align"]["sambamba"]["mem"],
        hrs=config["align"]["sambamba"]["hrs"],
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"sambamba/{SAMBAMBA_VERSION}",
    shell:
        """
        if [[ $( echo "{input.cell_bams}" | tr ' ' '\\n' | wc -l ) -eq 1 ]]; then
            mv {input.cell_bams} {output.merged_bam} && sambamba index --nthreads {threads} {output.merged_bam}
        else
            sambamba merge --nthreads {threads} {output.merged_bam} {input.cell_bams} && sambamba index --nthreads {threads} {output.merged_bam}
        fi
        """


rule split_bam_by_chrom:
    input:
        bam="results/{tech}/{ref}/align/phased/non-trio/minimap2/{sample}/{sample}_sorted-linked.bam",
        bai="results/{tech}/{ref}/align/phased/non-trio/minimap2/{sample}/{sample}_sorted-linked.bam.bai",
    output:
        chr_bam=temp(
            "results/{tech}/{ref}/align/phased/non-trio/minimap2/{sample}/{sample}_sorted-linked-{chr}.bam"
        ),
        chr_bam_bai=temp(
            "results/{tech}/{ref}/align/phased/non-trio/minimap2/{sample}/{sample}_sorted-linked-{chr}.bam.bai"
        ),
    threads: config["variant_call"]["samtools"]["threads"]
    resources:
        mem=lambda wildcards, attempt: attempt * config["variant_call"]["samtools"]["mem"],
        hrs=config["variant_call"]["samtools"]["hrs"],
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"samtools/{SAMTOOLS_VERSION}",
    shell:
        """
        samtools view -@ {threads} -h -b {input.bam} -o {output.chr_bam} {wildcards.chr} && samtools index -@ {threads} {output.chr_bam}
        """


rule clair3:
    input:
        ref=get_reference,
        chrom_bam=rules.split_bam_by_chrom.output.chr_bam,
        chrom_bam_bai=rules.split_bam_by_chrom.output.chr_bam_bai,
    output:
        chrom_vcf_gz=temp(
            "results/{tech}/{ref}/variant_call/clair3/{sample}/{chr}/merge_output.vcf.gz"
        )
    threads: config["variant_call"]["clair3"]["threads"]
    resources:
        mem=lambda wildcards, attempt: attempt * config["variant_call"]["clair3"]["mem"],
        hrs=config["variant_call"]["clair3"]["hrs"],
    container:
        CLAIR3_CNTR
    params:
        guppy_model = config["clair3"]["guppy_model"],
    shell:
        """
        /opt/bin/run_clair3.sh \
        --bam_fn={input.chrom_bam} \
        --ref_fn={input.ref} \
        --threads={threads} \
        --platform="{wildcards.tech}" \
        --model_path="/opt/models/{params.guppy_model}" \
        --output=$(dirname {output.chrom_vcf_gz}) 
        """


rule sniffles:
    input:
        ref=get_reference,
        chrom_bam=rules.split_bam_by_chrom.output.chr_bam,
        chrom_bam_bai=rules.split_bam_by_chrom.output.chr_bam_bai,
    output:
        vcf=temp(
            "results/{tech}/{ref}/variant_call/sniffles/{sample}/{chr}/{sample}_sniffles.vcf.gz"
        ),
        vcf_tbi=temp(
            "results/{tech}/{ref}/variant_call/sniffles/{sample}/{chr}/{sample}_sniffles.vcf.gz.tbi"
        ),
    params:
        trf = get_pipeline_resources(caller="sniffles", which_one="trf")
    threads: config["variant_call"]["sniffles"]["threads"]
    resources:
        mem=lambda wildcards, attempt: attempt * config["variant_call"]["sniffles"]["mem"],
        hrs=config["variant_call"]["sniffles"]["hrs"],
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"sniffles/{SNIFFLES_VERSION}",
    shell:
        """
        sniffles --threads {threads} -i {input.chrom_bam} {params.trf} --output-rnames -v {output.vcf}
        """


rule merge_chr_calls:
    input:
        chroms_vcf_gz=get_by_chrom,
    output:
        merged_vcf="results/{tech}/{ref}/variant_call/{caller}/{sample}/{sample}_{caller}.vcf.gz",
        merged_vcf_tbi="results/{tech}/{ref}/variant_call/{caller}/{sample}/{sample}_{caller}.vcf.gz.tbi"
    wildcard_constraints:
        caller="clair3|sniffles",
    threads: config["variant_call"]["bcftools"]["threads"]
    resources:
        mem=lambda wildcards, attempt: attempt * config["variant_call"]["bcftools"]["mem"],
        hrs=config["variant_call"]["bcftools"]["hrs"],
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"bcftools/{BCFTOOLS_VERSION}",
    shell:
        """
        bcftools concat {input.chroms_vcf_gz} -Oz -o {output.merged_vcf} && bcftools index --tbi --output {output.merged_vcf_tbi} {output.merged_vcf} 
        """


rule longphase:
    input:
        snv_vcf="results/{tech}/{ref}/variant_call/clair3/{sample}/{sample}_clair3.vcf.gz",
        bam="results/{tech}/{ref}/align/phased/non-trio/minimap2/{sample}/{sample}_sorted-linked.bam",
        bai="results/{tech}/{ref}/align/phased/non-trio/minimap2/{sample}/{sample}_sorted-linked.bam.bai",
        sv_vcf="results/{tech}/{ref}/variant_call/sniffles/{sample}/{sample}_sniffles.vcf.gz",
        ref=get_reference,
    output:
        phased_sv_vcf="results/{tech}/{ref}/variant_call/longphase/{sample}/{sample}_SV.vcf",
        phased_snp_vcf="results/{tech}/{ref}/variant_call/longphase/{sample}/{sample}.vcf",
    params:
        tech_arg = lambda wildcards: "pb" if wildcards.tech == "hifi" else wildcards.tech
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"longphase/{LONGPHASE_VERSION}",
    threads: config["variant_call"]["longphase"]["threads"]
    resources:
        mem=lambda wildcards, attempt: attempt * config["variant_call"]["longphase"]["mem"],
        hrs=config["variant_call"]["longphase"]["hrs"],
    shell:
        """
        longphase phase \
            -s {input.snv_vcf} \
            --sv-file {input.sv_vcf} \
            -b {input.bam} \
            -r {input.ref} \
            -t {threads} \
            --out-prefix $( echo {output.phased_snp_vcf} | sed 's/\.vcf//' ) \
            --{params.tech_arg}
        """


rule haplotag:
    input:
        snp_vcf=rules.longphase.output.phased_snp_vcf,
        sv_vcf=rules.longphase.output.phased_sv_vcf,
        bam="results/{tech}/{ref}/align/phased/non-trio/minimap2/{sample}/{sample}_sorted-linked.bam",
        bam_bai="results/{tech}/{ref}/align/phased/non-trio/minimap2/{sample}/{sample}_sorted-linked.bam.bai",
    output:
        haplotagged_bam="results/{tech}/{ref}/align/phased/non-trio/longphase/{sample}/{sample}_haplotagged_sorted-linked.bam",
        haplotagged_bam_bai="results/{tech}/{ref}/align/phased/non-trio/longphase/{sample}/{sample}_haplotagged_sorted-linked.bam.bai",
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"longphase/{LONGPHASE_VERSION}",
        f"samtools/{SAMTOOLS_VERSION}",
    threads: config["variant_call"]["longphase"]["threads"]
    resources:
        mem=lambda wildcards, attempt: attempt * config["variant_call"]["longphase"]["mem"],
        hrs=config["variant_call"]["longphase"]["hrs"],
    shell:
        """
        longphase haplotag \
            -s {input.snp_vcf} \
            --sv-file {input.sv_vcf} \
            -b {input.bam} \
            -t {threads} \
            --out-prefix $( echo {output.haplotagged_bam} | sed 's/\.bam//' ) \
        && samtools index -@ {threads} {output.haplotagged_bam}
        """
