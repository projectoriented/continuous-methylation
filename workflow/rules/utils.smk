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
