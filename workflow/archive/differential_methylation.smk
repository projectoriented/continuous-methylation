# -------- Rule-specific helper functions -------- #

def get_dss_summaries_by_chrX(wildcards):
    targets = []

    subset_df = dss_df.loc[dss_df["group_name"] == wildcards.group_name]

    for row in subset_df.itertuples():

        if pd.isnull(row.maternal_illumina_fofn):
            phase_type = "non-trio"
        else:
            phase_type = "trio"

            targets.extend(
                [
                    f"results/{wildcards.tech}/analysis/methylation/{{ref}}/dss/{wildcards.group_name}/no-intersection/tmp/haplotype/{{sample}}_{row.sample}_{hap}_trio_chrX-{{analysis_type}}_{{dss_category}}-summary.tsv"
                    for hap in HAPS
                ]
            )

        targets.append(f"results/{wildcards.tech}/analysis/methylation/{{ref}}/dss/{wildcards.group_name}/no-intersection/tmp/{{sample}}_{row.sample}_{phase_type}_chrX-{{analysis_type}}_{{dss_category}}-summary.tsv")

    return targets


def get_dss_summaries_by_chrom(wildcards):
    targets = []

    subset_df = dss_df.loc[dss_df["group_name"] == wildcards.group_name]

    for row in subset_df.itertuples():
        if pd.isnull(row.maternal_illumina_fofn):
            phase_type = "non-trio"
        else:
            phase_type = "trio"

            targets.extend(
                [
                    f"results/{wildcards.tech}/analysis/methylation/{{ref}}/dss/{wildcards.group_name}/intersection/{{groupA}}/haplotype/{row.sample}_{hap}_trio_{{chr}}-{{analysis_type}}_{{dss_category}}-summary.tsv"
                    for hap in HAPS
                ]
            )

        targets.append(f"results/{wildcards.tech}/analysis/methylation/{{ref}}/dss/{wildcards.group_name}/intersection/{{groupA}}/{row.sample}_{phase_type}_{{chr}}-{{analysis_type}}_{{dss_category}}-summary.tsv")

    return targets


def get_dss_inputs(wildcards):

    reference = wildcards.ref.split("-")[0]

    file_names = []

    if wildcards.chr != "chrX":
        subset_df = dss_df.loc[dss_df["sample"].str.contains(fr"{wildcards.groupA}|{wildcards.groupB}")]
        file_pattern = [f"results/{wildcards.tech}/analysis/methylation/{{ref}}/dss/txt/{wildcards.group_name}/{{sample}}_cpg-pileup_{{phase_type}}_{{{{chr}}}}.txt"]
    else:
        subset_df = dss_df.query(fr"sample == '{wildcards.sample}'")
        file_pattern = [ f"results/{wildcards.tech}/analysis/methylation/{{ref}}/dss/txt/{wildcards.group_name}/{{sample}}_{suffix}_{{phase_type}}_{{{{chr}}}}.txt" for suffix in ['hap1_cpg-pileup', 'hap2_cpg-pileup']]

    for row in subset_df.itertuples():
        if pd.isnull(row.maternal_illumina_fofn):
            phase_type = "non-trio"
        else:
            phase_type = "trio"

        if reference in row.reference_name:
            file_names.extend(
                [
                    f.format(ref=row.reference_name, sample=row.sample, phase_type=phase_type) for f in file_pattern
                ]
            )
        else:
            pass

    return file_names

def get_dss_params(wildcards):

    if wildcards.chr != "chrX":
        param_dict = {
            "groupA": wildcards.groupA,
            "groupB": wildcards.groupB,
            "sample_names": dss_df.loc[dss_df["sample"].str.contains(fr"{wildcards.groupA}|{wildcards.groupB}"), "sample"].tolist(),
        }
    else:
        sn = [f"{wildcards.sample}_{x}" for x in ["pat", "mat"]]
        param_dict = {
            "groupA": sn[0],
            "groupB": sn[1],
            "sample_names": sn,
        }

    return param_dict


def get_dss_groups(wildcards):

    path_pattern = "results/{{tech}}/analysis/methylation/{{ref}}/dss/{group_name}/original/{{analysis_type}}/{{chr}}_{groupA}_vs_{groupB}_DMR.tsv"

    subset_df = dss_df.query(fr"group_name == '{wildcards.group_name}'")

    group_a_name = subset_df.query("group == 'A'")["sample"]

    group_b = subset_df.query("group == 'B'")

    priority_b = group_b[group_b["maternal_illumina_fofn"].notnull()]["sample"]
    group_b = group_b[group_b["maternal_illumina_fofn"].isnull()]["sample"]

    b_list = [path_pattern.format(groupA=a, groupB=b, group_name=wildcards.group_name) for a in group_a_name for b in group_b]

    # Make sure A is proband versus sibling or A is the trio-phased sample declared as A
    a_list = [path_pattern.format(groupA=a, groupB=b, group_name=wildcards.group_name) for a in group_a_name for b in priority_b]

    return {
        "a": a_list,
        "b": b_list
    }

def get_dss_prepare_txt_output_chrX(which_one):
    def inner(wildcards):

        reference = manifest_df.loc[(manifest_df["sample"] == wildcards.other_sample) & (manifest_df["reference_name"].str.contains(str(wildcards.ref))), "reference_name"][0]

        if which_one == "haplotype":
            return f"results/{wildcards.tech}/analysis/methylation/{reference}/dss/txt/{{group_name}}/{wildcards.other_sample}_{{hap}}_cpg-pileup_{{phase_type}}_chrX.txt"
        elif which_one == "non-haplotype":
            return f"results/{wildcards.tech}/analysis/methylation/{reference}/dss/txt/{{group_name}}/{wildcards.other_sample}_cpg-pileup_{{phase_type}}_chrX.txt"
        else:
            sys.exit(f"Unsupported param: {which_one}")

    return inner


def get_dss_prepare_txt_output(which_one):
    def inner(wildcards):

        reference = manifest_df.loc[(manifest_df["sample"] == wildcards.sample) & (
            manifest_df["reference_name"].str.contains(str(wildcards.ref))), "reference_name"][0]

        if which_one == "haplotype":
            return f"results/{wildcards.tech}/analysis/methylation/{reference}/dss/txt/{{group_name}}/{{sample}}_{{hap}}_cpg-pileup_{{phase_type}}_{{chr}}.txt"
        elif which_one == "non-haplotype":
            return f"results/{wildcards.tech}/analysis/methylation/{reference}/dss/txt/{{group_name}}/{{sample}}_cpg-pileup_{{phase_type}}_{{chr}}.txt"
        else:
            sys.exit(f"Unsupported param: {which_one}")

    return inner

def get_anno_dict(wildcards):
    return config["annotations"][wildcards.ref]

def get_tech_specific_bed(wildcards):
    if wildcards.tech == "hifi":
        return "results/hifi/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_{suffix}.combined.bed.gz"
    else:
        return "results/{tech}/{ref}/methylation/phased/{phase_type}/{sample}/{sample}_{suffix}.bed.gz"

will_exclude_regions = False

if isinstance(config["dss_exclude_regions"], str):
    if os.path.exists(config["dss_exclude_regions"]):
        will_exclude_regions = True
        rule exclude_regions:
            input:
                bed = get_tech_specific_bed,
                regions_to_exclude = config["dss_exclude_regions"]
            output:
                excluded_regions = temp("results/{tech}/{ref}/methylation/phased/{phase_type}/{sample}/filtered/{sample}_{suffix}.bed.gz")
            wildcard_constraints:
                suffix="cpg-pileup|hap1_cpg-pileup|hap2_cpg-pileup"
            threads: 1
            resources:
                mem=lambda wildcards, attempt: attempt * 4,
                hrs=72
            envmodules:
                "modules",
                "modules-init",
                "modules-gs/prod",
                "modules-eichler/prod",
                f"bedtools/{BEDTOOLS_VERSION}"
            shell:
                """
                bedtools subtract -a {input.bed} -b {input.regions_to_exclude} | gzip -c > {output.excluded_regions}
                """

rule dss_prepare_in_txt:
    input:
        bed=get_tech_specific_bed if not will_exclude_regions else rules.exclude_regions.output.excluded_regions,
    output:
        bed_by_chrom=temp("results/{tech}/analysis/methylation/{ref}/dss/txt/{group_name}/{sample}_{suffix}_{phase_type}_{chr}.txt")
    wildcard_constraints:
        suffix="cpg-pileup|hap1_cpg-pileup|hap2_cpg-pileup"
    threads: 1
    resources:
        mem=calc_mem_gb,
        hrs=72
    params:
        min_mod = config.get("dss_min_mod", 0)
    shell:
        """
        if [ {wildcards.tech} == "hifi" ]; then
            # Columns grabbed are based on this documentation: https://github.com/PacificBiosciences/pb-CpG-tools#bed-file-format
            
            # chrom, start, n_valid, n_mod
            zgrep -w {wildcards.chr} {input.bed} | awk '$7 >= {params.min_mod} {{print $1,$2,$6,$7}}' FS='\\t' OFS='\\t' > {output.bed_by_chrom}
        elif [ {wildcards.tech} == "ont" ]; then
            # Columns grabbed are based on this documentation: https://github.com/nanoporetech/modkit/#bedmethyl-column-descriptions
            
            # chrom, start, n_valid, n_mod
            zgrep -w {wildcards.chr} {input.bed} | awk '$12 >= {params.min_mod} {{print $1,$2,$10,$12}}' FS='\\t' OFS='\\t' > {output.bed_by_chrom}
        else
            echo "Invalid tech wildcard: {wildcards.tech}" 1>&2; exit 1 
        fi
        """

rule dss_autosomes:
    input:
        file_names = get_dss_inputs
    output:
        dmr="results/{tech}/analysis/methylation/{ref}/dss/{group_name}/original/{analysis_type}/{chr}_{groupA}_vs_{groupB}_DMR.tsv",
        dml="results/{tech}/analysis/methylation/{ref}/dss/{group_name}/original/{analysis_type}/{chr}_{groupA}_vs_{groupB}_DML.tsv"
    wildcard_constraints:
        chr="chr[0-9]+"
    params:
        get_dss_params,
        output_prefix="results/{tech}/analysis/methylation/{ref}/dss/{group_name}/original/{analysis_type}/{chr}_{groupA}_vs_{groupB}"
    threads: 8
    resources:
        mem=calc_mem_gb,
        hrs=72
    log: "results/{tech}/analysis/methylation/{ref}/dss/{group_name}/original/{analysis_type}/log/{chr}_{groupA}_vs_{groupB}.log"
    container:
        DSS_CNTR
    script:
        "../scripts/DSS.R"

rule dss_chrX:
    input:
        file_names = get_dss_inputs
    output:
        dmr="results/{tech}/analysis/methylation/{ref}/dss/{group_name}/original/{chr}/{analysis_type}/{sample}_{groupA}_vs_{groupB}_DMR.tsv",
        dml="results/{tech}/analysis/methylation/{ref}/dss/{group_name}/original/{chr}/{analysis_type}/{sample}_{groupA}_vs_{groupB}_DML.tsv"
    wildcard_constraints:
        chr="chrX",
        groupA="pat",
        groupB="mat"
    params:
        get_dss_params,
        output_prefix="results/{tech}/analysis/methylation/{ref}/dss/{group_name}/original/{chr}/{analysis_type}/{sample}_{groupA}_vs_{groupB}"
    threads: 8
    resources:
        mem=calc_mem_gb,
        hrs=72
    log: "results/{tech}/analysis/methylation/{ref}/dss/{group_name}/original/{chr}/{analysis_type}/log/{sample}_{groupA}_vs_{groupB}.log"
    container:
        DSS_CNTR
    script:
        "../scripts/DSS.R"


rule softlink_chrX:
    input:
        dss_out = "results/{tech}/analysis/methylation/{ref}/dss/{group_name}/original/{chr}/{analysis_type}/{sample}_pat_vs_mat_DMR.tsv"
    output:
        softlinked_dmr = "results/{tech}/analysis/methylation/{ref}/dss/{group_name}/no-intersection/{chr}_{sample}_{analysis_type}_DMR.tsv"
    wildcard_constraints:
        chr="chrX"
    threads: 1
    resources:
        mem=calc_mem_gb,
        hrs=72
    shell:
        """
        ln -sf $(readlink -f {input.dss_out}) {output.softlinked_dmr}
        """

rule dss_overlap:
    input:
        unpack(get_dss_groups)
    output:
        dmr="results/{tech}/analysis/methylation/{ref}/dss/{group_name}/intersection/{groupA}/{chr}_{analysis_type}_DMR.tsv"
    threads: 1
    resources:
        mem=calc_mem_gb,
        hrs=72
    log: "results/{tech}/analysis/methylation/{ref}/dss/{group_name}/intersection/{groupA}/{chr}_{analysis_type}_DMR.log"
    run:
        from pybedtools import BedTool

        # logging
        sys.stderr = open(log[0],"w")

        a = input.a[0] if isinstance(input.a, list) else input.a
        b = input.b

        # Doing this only to get the header name.
        a_df = pd.read_table(a, header=0, dtype={"start": int, "end": int})
        target_cols = a_df.columns

        # BedTool
        a_bed = BedTool.from_dataframe(a_df)
        b_dict = {idx: BedTool.from_dataframe(pd.read_table(fp, header=0, dtype={"start": int, "end": int})) for idx, fp in enumerate(b)}

        # intersection_df = ( a_bed + list(b_dict.values()) ).to_dataframe(disable_auto_names=True, names=target_cols,header=None)
        intersection_df = a_bed.intersect(b_dict.get(0) + b_dict.get(1), wa=False).to_dataframe(disable_auto_names=True,names=target_cols,header=None)

        if intersection_df.empty:
            print(f"No intersection happened for {wildcards.groupA} in {wildcards.chr}", file=sys.stderr)
            intersection_df = pd.concat([ pd.read_table(x, header=0, dtype={"start": int, "end": int}) for x in [a] + b])

        intersection_df.to_csv(output.dmr, sep='\t', header=True, index=False)



rule dss_summary_table_chrX:
    input:
        dss_out = "results/{tech}/analysis/methylation/{ref}/dss/{group_name}/no-intersection/chrX_{sample}_{analysis_type}_DMR.tsv",
        sample_bed = get_dss_prepare_txt_output_chrX(which_one="non-haplotype")
    output:
        sample_summary_table = temp("results/{tech}/analysis/methylation/{ref}/dss/{group_name}/no-intersection/tmp/{sample}_{other_sample}_{phase_type}_chrX-{analysis_type}_{dss_category}-summary.tsv")
    threads: 1
    resources:
        mem= calc_mem_gb,
        hrs=72
    script:
        "../scripts/calculate_percent_meth.py"


rule dss_summary_table_chrX_haplotype:
    input:
        dss_out = "results/{tech}/analysis/methylation/{ref}/dss/{group_name}/no-intersection/chrX_{sample}_{analysis_type}_DMR.tsv",
        sample_bed = get_dss_prepare_txt_output_chrX(which_one="haplotype")
    output:
        sample_summary_table = temp("results/{tech}/analysis/methylation/{ref}/dss/{group_name}/no-intersection/tmp/haplotype/{sample}_{other_sample}_{hap}_{phase_type}_chrX-{analysis_type}_{dss_category}-summary.tsv")
    threads: 1
    resources:
        mem=calc_mem_gb,
        hrs=72
    script:
        "../scripts/calculate_percent_meth.py"


rule dss_summary_table:
    input:
        dss_out="results/{tech}/analysis/methylation/{ref}/dss/{group_name}/intersection/{groupA}/{chr}_{analysis_type}_{dss_category}.tsv",
        sample_bed=get_dss_prepare_txt_output(which_one="non-haplotype"),
    output:
        sample_summary_table=temp("results/{tech}/analysis/methylation/{ref}/dss/{group_name}/intersection/{groupA}/{sample}_{phase_type}_{chr}-{analysis_type}_{dss_category}-summary.tsv")
    threads: 1
    resources:
        mem=calc_mem_gb,
        hrs=72
    script:
        "../scripts/calculate_percent_meth.py"


rule dss_summary_table_by_haplotype:
    input:
        dss_out="results/{tech}/analysis/methylation/{ref}/dss/{group_name}/intersection/{groupA}/{chr}_{analysis_type}_{dss_category}.tsv",
        sample_bed=get_dss_prepare_txt_output(which_one="haplotype"),
    output:
        sample_summary_table=temp("results/{tech}/analysis/methylation/{ref}/dss/{group_name}/intersection/{groupA}/haplotype/{sample}_{hap}_{phase_type}_{chr}-{analysis_type}_{dss_category}-summary.tsv")
    threads: 1
    resources:
        mem=calc_mem_gb,
        hrs=72
    script:
        "../scripts/calculate_percent_meth.py"


rule merge_dss_summary_table_by_group:
    input:
        all_samples=get_dss_summaries_by_chrom,
    output:
        chrom_summary_table="results/{tech}/analysis/methylation/{ref}/dss/{group_name}/intersection/{groupA}/{chr}_{analysis_type}_{dss_category}-summary.tsv.gz"
    threads: 1
    resources:
        mem=calc_mem_gb,
        hrs=72
    script:
        "../scripts/merge_dss_summary_tables.py"


rule merge_dss_summary_table_by_sample:
    input:
        all_samples=get_dss_summaries_by_chrX,
    output:
        chrom_summary_table="results/{tech}/analysis/methylation/{ref}/dss/{group_name}/no-intersection/{sample}_chrX_{analysis_type}_{dss_category}-summary.tsv.gz"
    threads: 1
    resources:
        mem=calc_mem_gb,
        hrs=72
    script:
        "../scripts/merge_dss_summary_tables.py"

rule add_annotation:
    input:
        unpack(get_anno_dict),
        merged_summary = "results/{tech}/analysis/methylation/{ref}/dss/{group_name}/{intersect_strategy}/{prefix}_{analysis_type}_{dss_category}-summary.tsv.gz",
    output:
        added_annotation = "results/{tech}/analysis/methylation/{ref}/dss/{group_name}/{intersect_strategy}/{prefix}_{analysis_type}_{dss_category}-annotated-summary.tsv.gz"
    threads: 1
    resources:
        mem=calc_mem_gb,
        hrs=72
    run:
        from pybedtools import BedTool
        import pandas as pd

        df = pd.read_table(input.merged_summary,header=0)
        df["id"] = df.apply(lambda row: f"{row.chr}-{row.start}-{row.end}-{row.length}",axis=1)
        mod_df = df.loc[:, ['chr', 'start', 'end', 'id']].sort_values(by=["chr", "start", "end"])

        starting_bed = BedTool.from_dataframe(mod_df)

        keys = [k for k in input.keys() if k != "merged_summary"]

        for key in keys:

            anno_path = input.get(key)
            anno_df = pd.read_table(anno_path,header=None,names=["chr", "pos", "end", key]).sort_values(by=["chr", "pos", "end"])
            anno_bed = BedTool.from_dataframe(anno_df)

            if "gene" in key.lower():
                closest_genes_df = starting_bed.closest(anno_bed,D="ref").groupby(g=[4,9], c=8, o="collapse").to_dataframe(disable_auto_names=True, names=[
                    'id', 'distance(bp)', key],header=None)

                for idx, row in closest_genes_df.iterrows():
                    if row[key] == ".":
                        closest_genes_df.loc[idx, key] = "N/A"
                        closest_genes_df.loc[idx, "distance(bp)"] = "N/A"
                    else:
                        closest_genes_df.loc[idx, key] = ",".join(pd.Series(row[key].split(",")).drop_duplicates().to_list())

                closest_genes_df = closest_genes_df[['id', key, 'distance(bp)']]

                df = df.merge(closest_genes_df,on="id",how="left")

                del closest_genes_df

            else:
                hits = starting_bed.intersect(anno_bed,wa=True,wb=True)

                try:
                    len(hits[0])

                    annotated_df = hits.groupby(g=[4], c=8, o="collapse").to_dataframe(disable_auto_names=True, names=['id', key], header=None)
                    annotated_df[key] = annotated_df.apply(lambda x: ",".join(pd.Series(x[key].split(",")).drop_duplicates().to_list()),axis=1)
                    annotated_df = annotated_df[['id', key]]

                    df = df.merge(annotated_df,on="id",how="left")
                    del annotated_df
                except IndexError:
                    df[key] = "N/A"

                del hits

        df.drop(columns=["id"],inplace=True)
        df.fillna("N/A", inplace=True)

        df.to_csv(output.added_annotation, sep='\t', header=True, index=False)
