rule get_chrom_sizes:
    input:
        ref_fai=get_reference_fai,
    output:
        chrom_sizes="results/{tech}/{ref}/chrom.sizes"
    threads: 1
    resources:
        mem=calc_mem_gb,
        hrs=72,
    shell:
        """
        cut -f1,2 {input.ref_fai} > {output.chrom_sizes}
        """

rule remove_vc_unwanted:
    input:
        merged_vcf="results/{tech}/{ref}/variant_call/clair3/{sample}/{sample}_clair3.vcf.gz",
        merged_vcf_tbi="results/{tech}/{ref}/variant_call/clair3/{sample}/{sample}_clair3.vcf.gz.tbi"
    output:
        done = temp("results/{tech}/{ref}/variant_call/clair3/{sample}/.cleaned.txt")
    threads: 1
    resources:
        mem=calc_mem_gb,
        hrs=72,
    shell:
        """
        parent_dir=$(dirname {input.merged_vcf})
        find $parent_dir -type d -regex ".*chr.*" | xargs rm -rf && touch {output.done}
        """

# rule fastq:
#     input:
#         bam = get_bam,
#     output:
#         fastq_gz = temp("fastq/{sample}.fastq.gz")
#     envmodules:
#         "modules",
#         "modules-init",
#         "modules-gs/prod",
#         "modules-eichler/prod",
#         "samtools/1.14"
#     threads: 4
#     resources:
#         mem=calc_mem_gb,
#         hrs=72,
#     shell:
#         """
# 	    bgzip -@ {threads} < <(samtools fastq {input.bam} -@ {threads}) > {output.fastq_gz}
#         """
#
# rule split_fastq:
#     input:
#         fastq_gz = "fastq/{sample}.fastq.gz"
#     output:
#         fastq_gz_parts = temp(expand("fastq/parts/{{sample}}_{part}.fastq.gz", part=[f"{x:03d}" for x in range(1, fastq_split_parts)]))
#     params:
#         fastq_parts = fastq_split_parts
#     envmodules:
#         "modules",
#         "modules-init",
#         "modules-gs/prod",
#         "modules-eichler/prod",
#         "seqkit/2.6.1"
#     threads: 4
#     resources:
#         mem=calc_mem_gb,
#         hrs=72,
#     shell:
#         """
# 	    seqkit split2 {input.fastq_gz} --by-part {params.fastq_parts} --extension .gz --by-part-prefix {wildcards.sample}_ --out-dir $(dirname {output.fastq_gz_parts[0]})
#         """
