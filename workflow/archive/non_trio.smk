rule haplotaggedness_bam:
    input:
        bam=rules.haplotag.output.haplotagged_bam,
        bam_bai=rules.haplotag.output.haplotagged_bam_bai,
    output:
        hap1_bam=temp("results/{tech}/{ref}/align/phased/non-trio/longphase/{sample}/{sample}_hap1_haplotagged_sorted-linked.bam"),
        hap2_bam=temp("results/{tech}/{ref}/align/phased/non-trio/longphase/{sample}/{sample}_hap2_haplotagged_sorted-linked.bam"),
        hap_unknown_bam=temp("results/{tech}/{ref}/align/phased/non-trio/longphase/{sample}/{sample}_unknown_haplotagged_sorted-linked.bam"),
        hap1_bam_bai=temp("results/{tech}/{ref}/align/phased/non-trio/longphase/{sample}/{sample}_hap1_haplotagged_sorted-linked.bam.bai"),
        hap2_bam_bai=temp("results/{tech}/{ref}/align/phased/non-trio/longphase/{sample}/{sample}_hap2_haplotagged_sorted-linked.bam.bai"),
        hap_unknown_bam_bai=temp("results/{tech}/{ref}/align/phased/non-trio/longphase/{sample}/{sample}_unknown_haplotagged_sorted-linked.bam.bai"),
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"samtools/{SAMTOOLS_VERSION}",
    threads: 16
    resources:
        mem=calc_mem_gb,
        hrs=72,
    shell:
        """
        samtools view -e '[HP]==1' -@ {threads} -bh {input.bam} > {output.hap1_bam} && samtools index -@ {threads} {output.hap1_bam}
        samtools view -e '[HP]==2' -@ {threads} -bh {input.bam} > {output.hap2_bam} && samtools index -@ {threads} {output.hap2_bam}
        samtools view -e '![HP]' -@ {threads} -bh {input.bam} > {output.hap_unknown_bam} && samtools index -@ {threads} {output.hap_unknown_bam}
        """
