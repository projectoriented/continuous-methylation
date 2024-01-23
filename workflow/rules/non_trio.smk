rule split_bam_by_chrom:
    input:
        bam="results/{tech}/{ref}/align/phased/non-trio/minimap2/{sample}/{sample}_sorted-5mC.bam",
        csi="results/{tech}/{ref}/align/phased/non-trio/minimap2/{sample}/{sample}_sorted-5mC.bam.csi",
    output:
        chr_bam=temp(
            "results/{tech}/{ref}/align/phased/non-trio/minimap2/{sample}/{sample}_sorted-5mC-{chr}.bam"
        ),
        chr_bam_bai=temp(
            "results/{tech}/{ref}/align/phased/non-trio/minimap2/{sample}/{sample}_sorted-5mC-{chr}.bam.bai"
        ),
    threads: 4
    resources:
        mem=lambda wildcards, attempt: attempt * 8,
        hrs=72,
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
    threads: 16
    resources:
        mem=calc_mem_gb,
        hrs=72,
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
    threads: 16
    resources:
        mem=calc_mem_gb,
        hrs=72,
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
    threads: 1
    resources:
        mem=calc_mem_gb,
        hrs=72,
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
        bam="results/{tech}/{ref}/align/phased/non-trio/minimap2/{sample}/{sample}_sorted-5mC.bam",
        csi="results/{tech}/{ref}/align/phased/non-trio/minimap2/{sample}/{sample}_sorted-5mC.bam.csi",
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
    threads: 16
    resources:
        mem=calc_mem_gb,
        hrs=72,
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
        bam="results/{tech}/{ref}/align/phased/non-trio/minimap2/{sample}/{sample}_sorted-5mC.bam",
        bam_csi="results/{tech}/{ref}/align/phased/non-trio/minimap2/{sample}/{sample}_sorted-5mC.bam.csi",
    output:
        haplotagged_bam="results/{tech}/{ref}/align/phased/non-trio/longphase/{sample}/{sample}_sorted-5mC-haplotagged.bam",
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"longphase/{LONGPHASE_VERSION}",
    threads: 16
    resources:
        mem=lambda wildcards, attempt: attempt * 4,
        hrs=72,
    shell:
        """
        longphase haplotag \
            -s {input.snp_vcf} \
            --sv-file {input.sv_vcf} \
            -b {input.bam} \
            -t {threads} \
            --out-prefix $( echo {output.haplotagged_bam} | sed 's/\.bam//' )
        """

rule index_haplotag:
    input:
        haplotagged_bam="results/{tech}/{ref}/align/phased/non-trio/longphase/{sample}/{sample}_sorted-5mC-haplotagged.bam",
    output:
        haplotagged_bam_csi="results/{tech}/{ref}/align/phased/non-trio/longphase/{sample}/{sample}_sorted-5mC-haplotagged.bam.csi",
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"samtools/{SAMTOOLS_VERSION}",
    threads: 16
    resources:
        mem=lambda wildcards, attempt: attempt * 4,
        hrs=72,
    shell:
        """
        samtools index -c -@ {threads} {input.haplotagged_bam}
        """

# rule haplotaggedness_bam:
#     input:
#         bam=rules.haplotag.output.haplotagged_bam,
#         bam_bai=rules.haplotag.output.haplotagged_bam_bai,
#     output:
#         hap1_bam=temp("results/{tech}/{ref}/align/phased/non-trio/longphase/{sample}/{sample}_hap1_haplotagged_sorted-linked.bam"),
#         hap2_bam=temp("results/{tech}/{ref}/align/phased/non-trio/longphase/{sample}/{sample}_hap2_haplotagged_sorted-linked.bam"),
#         hap_unknown_bam=temp("results/{tech}/{ref}/align/phased/non-trio/longphase/{sample}/{sample}_unknown_haplotagged_sorted-linked.bam"),
#         hap1_bam_bai=temp("results/{tech}/{ref}/align/phased/non-trio/longphase/{sample}/{sample}_hap1_haplotagged_sorted-linked.bam.bai"),
#         hap2_bam_bai=temp("results/{tech}/{ref}/align/phased/non-trio/longphase/{sample}/{sample}_hap2_haplotagged_sorted-linked.bam.bai"),
#         hap_unknown_bam_bai=temp("results/{tech}/{ref}/align/phased/non-trio/longphase/{sample}/{sample}_unknown_haplotagged_sorted-linked.bam.bai"),
#     envmodules:
#         "modules",
#         "modules-init",
#         "modules-gs/prod",
#         "modules-eichler/prod",
#         f"samtools/{SAMTOOLS_VERSION}",
#     threads: 16
#     resources:
#         mem=calc_mem_gb,
#         hrs=72,
#     shell:
#         """
#         samtools view -e '[HP]==1' -@ {threads} -bh {input.bam} > {output.hap1_bam} && samtools index -@ {threads} {output.hap1_bam}
#         samtools view -e '[HP]==2' -@ {threads} -bh {input.bam} > {output.hap2_bam} && samtools index -@ {threads} {output.hap2_bam}
#         samtools view -e '![HP]' -@ {threads} -bh {input.bam} > {output.hap_unknown_bam} && samtools index -@ {threads} {output.hap_unknown_bam}
#         """
