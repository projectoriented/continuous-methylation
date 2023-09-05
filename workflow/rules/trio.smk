if config["tech"] == "ont":
    rule canu_trio:
        input:
            unpack(get_assembly_trio_inputs(which_one="parental_illumina")),
            ont_fastq=get_assembly_trio_inputs(which_one="fastq"),
        output:
            hap1="results/ont/assemblies/canu/trio/{sample}/haplotype/haplotype-1.fasta.gz",
            hap2="results/ont/assemblies/canu/trio/{sample}/haplotype/haplotype-2.fasta.gz"
        wildcard_constraints:
            phase_type="trio"
        params:
            canu_max_mem = int(config["assembly"]["canu"]["mem"]) * int(config["assembly"]["canu"]["threads"])
        threads: config["assembly"]["canu"]["threads"]
        resources:
            mem=lambda wildcards, attempt: attempt * config["assembly"]["canu"]["mem"],
            hrs=config["assembly"]["canu"]["hrs"],
        envmodules:
            "modules",
            "modules-init",
            "modules-gs/prod",
            "modules-eichler/prod",
            f"canu/{CANU_VERSION}",
        shell:
            """
            canu \
                -p {wildcards.sample} \
                -d $(echo {output.hap1} | sed 's/\/haplotype.*//') \
                -haplotype \
                genomeSize=3.1g \
                maxInputCoverage=150 \
                -haplotype1 $(cat {input.pat}) \
                -haplotype2 $(cat {input.mat}) \
                -nanopore $(cat {input.ont_fastq}) \
                maxThreads={threads} \
                maxMemory={params.canu_max_mem}g \
                useGrid=false
            """

    rule rename_canu_output:
        input:
            hap1 = "results/ont/assemblies/canu/trio/{sample}/haplotype/haplotype-1.fasta.gz",
            hap2 = "results/ont/assemblies/canu/trio/{sample}/haplotype/haplotype-2.fasta.gz"
        output:
            hap1 = "results/ont/assemblies/canu/trio/{sample}/haplotype/{sample}_hap1.fasta.gz",
            hap2 = "results/ont/assemblies/canu/trio/{sample}/haplotype/{sample}_hap2.fasta.gz"
        threads: config["default"]["threads"]
        resources:
            mem=lambda wildcards, attempt: attempt * config["default"]["mem"],
            hrs=config["default"]["hrs"],
        shell:
            """
            mv {input.hap1} {output.hap1}
            mv {input.hap2} {output.hap2}
            """
elif config["tech"] == "hifi":
    rule yak_parents:
        input:
            parental_illumina=get_yak_input,
        output:
            parental_yak=temp("results/hifi/yak/parents/{family}/{parental}.yak"),
        threads: config["assembly"]["yak"]["threads"]
        resources:
            mem=lambda wildcards, attempt: attempt * config["assembly"]["yak"]["mem"],
            hrs=config["assembly"]["yak"]["hrs"],
        envmodules:
            "modules",
            "modules-init",
            "modules-gs/prod",
            "modules-eichler/prod",
            f"yak/{YAK_VERSION}",
        shell:
            """
            yak count \
                -k 31 -b 37 \
                -t {threads} \
                -o {output.parental_yak} \
                <(cat $(cat {input.parental_illumina}) ) <(cat $(cat {input.parental_illumina}) )
            """

    rule hifiasm_dual:
        input:
            fastq=get_assembly_trio_inputs(which_one="fastq"),
        output:
            asm_hap1=temp("results/hifi/assemblies/hifiasm/trio/{sample}.hifiasm.bp.hap1.p_ctg.gfa"),
            asm_hap2=temp("results/hifi/assemblies/hifiasm/trio/{sample}.hifiasm.bp.hap2.p_ctg.gfa")
        threads: config["assembly"]["hifiasm"]["threads"]
        resources:
            mem=lambda wildcards, attempt: attempt * config["assembly"]["hifiasm"]["mem"],
            hrs=config["assembly"]["hifiasm"]["hrs"],
        log:
            "results/hifi/assemblies/hifiasm/trio/{sample}-primary.log"
        envmodules:
            "modules",
            "modules-init",
            "modules-gs/prod",
            "modules-eichler/prod",
            f"hifiasm/{HIFIASM_VERSION}",
        shell:
            """
            hifiasm \
                -o results/hifi/assemblies/hifiasm/trio/{wildcards.sample}.hifiasm \
                -t {threads} \
                $( cat {input.fastq} ) 2> {log}
            """


    rule hifiasm_trio:
        input:
            unpack(get_parental_yak),
            asm=rules.hifiasm_dual.output.asm_hap1
        output:
            asm_pat="results/hifi/assemblies/hifiasm/trio/{sample}.hifiasm.dip.hap1.p_ctg.gfa",
            asm_mat="results/hifi/assemblies/hifiasm/trio/{sample}.hifiasm.dip.hap2.p_ctg.gfa",
        threads: config["assembly"]["hifiasm"]["threads"]
        resources:
            mem=lambda wildcards, attempt: attempt * config["assembly"]["hifiasm"]["mem"],
            hrs=config["assembly"]["hifiasm"]["hrs"],
        log:
            "results/hifi/assemblies/hifiasm/trio/{sample}-trio.log"
        envmodules:
            "modules",
            "modules-init",
            "modules-gs/prod",
            "modules-eichler/prod",
            f"hifiasm/{HIFIASM_VERSION}",
        shell:
            """
            hifiasm \
                -t {threads} \
                -o results/hifi/assemblies/hifiasm/trio/{wildcards.sample}.hifiasm \
                -1 {input.pat} \
                -2 {input.mat} \
                /dev/null 2> {log}
            """

rule extract_hap_reads:
    input:
        hap_fa=lambda wildcards: "results/{tech}/assemblies/canu/trio/{sample}/haplotype/{sample}_{hap}.fasta.gz" if wildcards.tech == "ont" else "results/{tech}/assemblies/hifiasm/trio/{sample}.hifiasm.dip.{hap}.p_ctg.gfa",
    output:
        hap_read_names=temp(
            "results/{tech}/{ref}/align/phased/trio/{sample}/fastq/{sample}_{hap}.txt"
        ),
    threads: config["default"]["threads"]
    resources:
        mem=lambda wildcards, attempt: attempt * config["default"]["mem"],
        hrs=config["default"]["hrs"],
    shell:
        """
        if [[ {wildcards.tech} == "ont" ]]; then
            zcat {input.hap_fa} | grep -E "^>" | sed -r 's/(>)(.*)( runid.*)/\\2/g' > {output.hap_read_names}
        else
            tail -n +2 {input.hap_fa} | cut -f5 > {output.hap_read_names}
        fi
        """


rule correspond_fastq_reads_to_hap:
    input:
        cell_fastq=get_tech_files(which_one="fofn"),
        hap_read_names=rules.extract_hap_reads.output.hap_read_names,
    output:
        cell_hap_fastq=temp(
            "results/{tech}/{ref}/align/phased/{phase_type}/{sample}/fastq/{cell}_{hap}.fastq.gz"
        ),
    wildcard_constraints:
        phase_type="trio"
    threads: config["default"]["threads"]
    resources:
        mem=lambda wildcards, attempt: attempt * (config["default"]["mem"] * 2),
        hrs=config["default"]["hrs"],
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"seqtk/{SEQTK_VERSION}",
    shell:
        """
        seqtk subseq {input.cell_fastq} {input.hap_read_names} | gzip -c > {output.cell_hap_fastq}
        """


rule align_hap_specific_reads:
    input:
        cell_hap_fastq=rules.correspond_fastq_reads_to_hap.output.cell_hap_fastq,
        ref=get_reference,
    output:
        cell_hap_bam=temp(
            "results/{tech}/{ref}/align/phased/{phase_type}/{sample}/{sample}_{cell}_{hap}_sorted.bam"
        ),
        cell_hap_bam_bai=temp(
            "results/{tech}/{ref}/align/phased/{phase_type}/{sample}/{sample}_{cell}_{hap}_sorted.bam.bai"
        ),
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
        minimap2 -t {threads} -I 10G -Y -y --secondary=no --eqx -a -x map-{wildcards.tech} {input.ref} {input.cell_hap_fastq} | sambamba view --sam-input --format bam /dev/stdin | sambamba sort --nthreads {threads} --out {output.cell_hap_bam} /dev/stdin && sambamba index --nthreads {threads} {output.cell_hap_bam}
        """


rule link_meth_tags:
    input:
        cell_hap_bam=rules.align_hap_specific_reads.output.cell_hap_bam,
        cell_hap_bam_bai=rules.align_hap_specific_reads.output.cell_hap_bam_bai,
        unmapped_bam=get_tech_files(which_one="unmapped_bam_fofn")
    output:
        cell_hap_linked_bam=temp("results/{tech}/{ref}/align/phased/{phase_type}/{sample}/{sample}_{cell}_{hap}_sorted-linked.bam"),
        cell_hap_linked_bam_bai=temp("results/{tech}/{ref}/align/phased/{phase_type}/{sample}/{sample}_{cell}_{hap}_sorted-linked.bam.bai")
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
          --aln {input.cell_hap_bam} \
          --sample {wildcards.sample}_{wildcards.cell}_{wildcards.hap} \
          --methyl_bams "$(echo {input.unmapped_bam})" \
          --output {output.cell_hap_linked_bam}
        """


rule combine_bam_by_hap:
    input:
        hap_bams=gather_tech_bams(which_one="hap_specific"),
    output:
        merged_hap_bam=temp("results/{tech}/{ref}/align/phased/{phase_type}/{sample}/{sample}_{hap}_sorted-linked.bam"),
        merged_hap_bam_bai=temp("results/{tech}/{ref}/align/phased/{phase_type}/{sample}/{sample}_{hap}_sorted-linked.bam.bai"),
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
        if [[ $( echo "{input.hap_bams}" | tr ' ' '\\n' | wc -l ) -eq 1 ]]; then
            mv {input.hap_bams} {output.merged_hap_bam} && sambamba index --nthreads {threads} {output.merged_hap_bam}
        else
            sambamba merge --nthreads {threads} {output.merged_hap_bam} {input.hap_bams} && sambamba index --nthreads {threads} {output.merged_hap_bam}
        fi
        """


rule merge_hap_bams:
    input:
        hap_bams = expand("results/{{tech}}/{{ref}}/align/phased/{{phase_type}}/{{sample}}/{{sample}}_{hap}_sorted-linked.bam", hap=HAPS)
    output:
        merged_bam=temp("results/{tech}/{ref}/align/phased/{phase_type}/{sample}/{sample}_sorted-linked.bam"),
        merged_bam_bai=temp("results/{tech}/{ref}/align/phased/{phase_type}/{sample}/{sample}_sorted-linked.bam.bai"),
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
        sambamba merge --nthreads {threads} {output.merged_bam} {input.hap_bams} && sambamba index --nthreads {threads} {output.merged_bam}
        """
