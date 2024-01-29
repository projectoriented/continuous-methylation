if config["tech"] == "ont":
    rule canu_trio:
        input:
            unpack(get_assembly_trio_inputs(which_one="parental_illumina")),
            ont_fastq=get_assembly_trio_inputs(which_one="fastq"),
        output:
            hap1="results/ont/assemblies/canu/trio/{sample}/haplotype/haplotype-1.fasta.gz",
            hap2="results/ont/assemblies/canu/trio/{sample}/haplotype/haplotype-2.fasta.gz",
            unknown_hap="results/ont/assemblies/canu/trio/{sample}/haplotype/haplotype-unknown.fasta.gz"
        wildcard_constraints:
            phase_type="trio"
        params:
            canu_max_mem = lambda wildcards, threads, resources: int(resources["mem"]) * int(threads)
        threads: 16
        resources:
            mem=lambda wildcards, attempt, input: attempt * 10,
            hrs=72,
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
            hap2 = "results/ont/assemblies/canu/trio/{sample}/haplotype/haplotype-2.fasta.gz",
            unknown_hap = "results/ont/assemblies/canu/trio/{sample}/haplotype/haplotype-unknown.fasta.gz"
        output:
            hap1 = "results/ont/assemblies/canu/trio/{sample}/haplotype/{sample}_hap1.fasta.gz",
            hap2 = "results/ont/assemblies/canu/trio/{sample}/haplotype/{sample}_hap2.fasta.gz",
            unknown_hap= "results/ont/assemblies/canu/trio/{sample}/haplotype/{sample}_unknown.fasta.gz"
        threads: 1
        resources:
            mem=1,
            hrs=72,
        shell:
            """
            mv {input.hap1} {output.hap1}
            mv {input.hap2} {output.hap2}
            mv {input.unknown_hap} {output.unknown_hap}
            """
elif config["tech"] == "hifi":
    rule yak_parents:
        input:
            parental_illumina=get_yak_input,
        output:
            parental_yak=temp("results/hifi/yak/parents/{family}/{parental}.yak"),
        threads: 8
        resources:
            mem=lambda wildcards, attempt, input: attempt * 10,
            hrs=72,
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
        threads: 16
        resources:
            mem=lambda wildcards, attempt: attempt * 6,
            hrs=72,
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
            asm_all="results/hifi/assemblies/hifiasm/trio/{sample}.hifiasm.dip.r_utg.noseq.gfa",
        threads: 16
        resources:
            mem=lambda wildcards, attempt: attempt * 4,
            hrs=72,
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
    threads: 1
    resources:
        mem=calc_mem_gb,
        hrs=72,
    shell:
        """
        if [[ {wildcards.tech} == "ont" ]]; then
            zcat {input.hap_fa} | grep -E "^>" | sed -r 's/(>)(.*)( runid.*)/\\2/g' > {output.hap_read_names}
        else
            tail -n +2 {input.hap_fa} | cut -f5 > {output.hap_read_names}
        fi
        """


rule extract_non_hap_reads:
    input:
        asm=lambda wildcards: "results/ont/assemblies/canu/trio/{sample}/haplotype/{sample}_unknown.fasta.gz" if wildcards.tech == "ont" else "results/hifi/assemblies/hifiasm/trio/{sample}.hifiasm.dip.r_utg.noseq.gfa",
    output:
        read_names=temp(
            "results/{tech}/{ref}/align/phased/trio/{sample}/fastq/{sample}_non-binnable.txt"
        ),
    threads: 1
    resources:
        mem=calc_mem_gb,
        hrs=72,
    shell:
        """
        if [[ {wildcards.tech} == "ont" ]]; then
            zcat {input.asm} | grep -E "^>" | sed -r 's/(>)(.*)( runid.*)/\\2/g' > {output.read_names}
        else
            grep -w "HG:A:a" {input.asm} | cut -f5 > {output.read_names}
        fi
        """

rule correspond_reads:
    input:
        cell_fastq=get_fastq(which_one="correspond_reads"),
        read_names="results/{tech}/{ref}/align/phased/trio/{sample}/fastq/{sample}_{suffix}.txt",
    output:
        fastq=temp(
            "results/{tech}/{ref}/align/phased/trio/{sample}/fastq/{sample}_{suffix}_{cell}.fastq.gz"
        ),
    wildcard_constraints:
        cell="|".join(get_all_cell_names())
    threads: 1
    resources:
        mem=calc_mem_gb,
        hrs=72,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"seqtk/{SEQTK_VERSION}",
    shell:
        """
        seqtk subseq {input.cell_fastq} {input.read_names} | gzip -c > {output.fastq}
        """

