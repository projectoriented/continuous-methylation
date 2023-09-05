rule get_chrom_sizes:
    input:
        ref_fai=get_reference + ".fai",
    output:
        chrom_sizes="results/{tech}/{ref}/chrom.sizes"
    threads: 1
    resources:
        mem=lambda wildcards, attempt: attempt * 4,
        hrs=72,
    shell:
        """
        cut -f1,2 {input.ref_fai} > {output.chrom_sizes}
        """