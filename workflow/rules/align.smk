rule align:
    input:
        fastq=get_fastq(which_one="other"),
        ref=get_reference,
    output:
        cell_bam=temp(
            "results/{tech}/{ref}/align/phased/{phase_type}/minimap2/{sample}/{sample}_{cell}_{suffix}.bam"
        ),
        cell_bam_bai=temp(
            "results/{tech}/{ref}/align/phased/{phase_type}/minimap2/{sample}/{sample}_{cell}_{suffix}.bam.bai"
        ),
    wildcard_constraints:
        suffix="hap1_sorted|hap2_sorted|non-binnable_sorted|sorted",
        cell="|".join(get_all_cell_names())
    params:
        tech_arg=lambda wildcards: f"map-{wildcards.tech}",
        mm2_params=MINIMAP2_PARAMS
    threads: 12
    resources:
        mem=calc_mem_gb,
        hrs=72,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"minimap2/{MINIMAP2_VERSION}",
        f"sambamba/{SAMBAMBA_VERSION}",
    shell:
        """
        # --MD tag is required by sniffles
        minimap2_args=""
        if [ {wildcards.phase_type} == "non-trio" ]
        then
            minimap2_args="--MD"
        else
            minimap2_args="-y --secondary=no"
        fi

        minimap2 \
            -t {threads} -Y --eqx \
            $minimap2_args \
            -ax {params.tech_arg} \
            {params.mm2_params} \
            {input.ref} {input.fastq} \
        | \
        sambamba view \
            --sam-input \
            --format bam \
            /dev/stdin \
        | \
        sambamba sort \
        --nthreads {threads} \
        --out {output.cell_bam} /dev/stdin \
        && sambamba index --nthreads {threads} {output.cell_bam}
        """


rule link_meth_tags:
    input:
        bam="results/{tech}/{ref}/align/phased/{phase_type}/minimap2/{sample}/{sample}_{cell}_{suffix}.bam",
        bam_bai="results/{tech}/{ref}/align/phased/{phase_type}/minimap2/{sample}/{sample}_{cell}_{suffix}.bam.bai",
        unmapped_bam=get_unmapped_bam
    output:
        linked_bam=temp("results/{tech}/{ref}/align/phased/{phase_type}/minimap2/{sample}/{sample}_{cell}_{suffix}-5mC.bam"),
        linked_bam_bai=temp("results/{tech}/{ref}/align/phased/{phase_type}/minimap2/{sample}/{sample}_{cell}_{suffix}-5mC.bam.bai"),
    wildcard_constraints:
        suffix="hap1_sorted|hap2_sorted|non-binnable_sorted|sorted",
        cell="|".join(get_all_cell_names())
    threads: 16
    resources:
        mem=calc_mem_gb,
        hrs=72,
        disk="250G"
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
          --tmp {resources.tmpdir} \
          --aln {input.bam} \
          --sample {wildcards.sample}_{wildcards.cell}_{wildcards.suffix} \
          --methyl_bams "$(echo {input.unmapped_bam})" \
          --output {output.linked_bam}
        """

rule haplotag_trio:
    input:
        bam="results/{tech}/{ref}/align/phased/trio/minimap2/{sample}/{sample}_{cell}_{suffix}-5mC.bam",
        bam_bai="results/{tech}/{ref}/align/phased/trio/minimap2/{sample}/{sample}_{cell}_{suffix}-5mC.bam.bai",
    output:
        haplotagged_bam=temp("results/{tech}/{ref}/align/phased/{phase_type}/minimap2/{sample}/{sample}_{cell}_{suffix}-5mC-haplotagged.bam"),
    wildcard_constraints:
        suffix="hap1_sorted|hap2_sorted|non-binnable_sorted|sorted",
        cell="|".join(get_all_cell_names())
    threads: 1
    resources:
        mem=calc_mem_gb,
        hrs=72,
    run:
        import pysam

        tag_dict = {
            "hap1_sorted": 1,
            "hap2_sorted": 2,
        }
        samfile = pysam.AlignmentFile(input.bam,check_sq=False)
        haplotag = pysam.AlignmentFile(output.haplotagged_bam,"wb",template=samfile)
        for read in samfile.fetch(until_eof=True):
            read.tags = read.tags + [("HP", tag_dict[wildcards.suffix])]
            haplotag.write(read)

        haplotag.close()
        # pysam.index(output.haplotagged_bam)
        samfile.close()

rule merge_align:
    input:
        bams=gather_tech_bams,
    output:
        bam="results/{tech}/{ref}/align/phased/{phase_type}/minimap2/{sample}/{sample}_sorted-{suffix}.bam",
    wildcard_constraints:
        suffix="5mC-haplotagged|5mC"
    threads: 8
    resources:
        mem=lambda wildcards, attempt: attempt * 4,
        hrs=72,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"samtools/{SAMTOOLS_VERSION}"
    shell:
        """
        if [[ $( echo "{input.bams}" | tr ' ' '\\n' | wc -l ) -eq 1 ]]; then
            mv {input.bams} {output.bam}
        else
            samtools merge -fr -@ {threads} -o {output.bam} {input.bams} 
        fi
        """

rule index_merge_align:
    input:
        bam="results/{tech}/{ref}/align/phased/{phase_type}/minimap2/{sample}/{sample}_sorted-{suffix}.bam",
    output:
        bam_csi="results/{tech}/{ref}/align/phased/{phase_type}/minimap2/{sample}/{sample}_sorted-{suffix}.bam.csi",
    wildcard_constraints:
        suffix="5mC-haplotagged|5mC"
    threads: 8
    resources:
        mem=lambda wildcards, attempt: attempt * 4,
        hrs=72,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"samtools/{SAMTOOLS_VERSION}"
    shell:
        """
        samtools index -c -@ {threads} {input.bam}
        """