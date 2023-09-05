import sys
import os

import pandas as pd

# --------  Switches -------- #
methylation_analysis = config.get("differential_methylation", None)

# --------  Load sample sheet -------- #
manifest_df = pd.read_table(config["manifest"], dtype=str, header=0)
manifest_df.set_index(["sample", "reference_name"], drop=False, inplace=True)

if methylation_analysis:
    dss_df = pd.read_table(config["dss_manifest"], dtype=str, header=0)
    dss_df = dss_df.merge(
        manifest_df.reset_index(drop=True),
        how="inner"
    ).set_index(["sample"],drop=False)

# -------- Wildcard constraints -------- #
wildcard_constraints:
    tech=TECH,
    sample="|".join(manifest_df.index.get_level_values("sample")),
    parental="pat|mat",
    hap="|".join(HAPS),
    phase_type="trio|non-trio",
    intersect_strategy="no-intersection|intersection",
    analysis_type="two_group|model_based",
    dss_category="DMR|DML",
    group_name="|".join(manifest_df["group_name"].unique())

# -------- Helper functions -------- #
def get_final_output(wildcards):

    final_output = []

    final_output.extend(get_methyl_targets())

    if methylation_analysis:
        final_output.extend(get_final_dss_targets())

    return final_output

def get_cpg_bams(wildcards):
    if wildcards.phase_type == "trio":
        if "hap" in wildcards.suffix:
            current_hap = wildcards.suffix.split("_")[0]
            return {
                "bam": f"results/{TECH}/{wildcards.ref}/align/phased/trio/{wildcards.sample}/{wildcards.sample}_{current_hap}_sorted-linked.bam",
                "bai": f"results/{TECH}/{wildcards.ref}/align/phased/trio/{wildcards.sample}/{wildcards.sample}_{current_hap}_sorted-linked.bam.bai"
            }
        else:
            return {
                "bam": f"results/{TECH}/{wildcards.ref}/align/phased/trio/{wildcards.sample}/{wildcards.sample}_sorted-linked.bam",
                "bai": f"results/{TECH}/{wildcards.ref}/align/phased/trio/{wildcards.sample}/{wildcards.sample}_sorted-linked.bam.bai"
            }
    elif wildcards.phase_type == "non-trio":
        if "hap" in wildcards.suffix or "unknown" in wildcards.suffix and TECH == "ont":
            current_hap = wildcards.suffix.split("_")[0]
            return {
                "bam": f"results/{TECH}/{wildcards.ref}/align/phased/non-trio/longphase/{wildcards.sample}/{wildcards.sample}_{current_hap}_haplotagged_sorted-linked.bam",
                "bai": f"results/{TECH}/{wildcards.ref}/align/phased/non-trio/longphase/{wildcards.sample}/{wildcards.sample}_{current_hap}_haplotagged_sorted-linked.bam.bai",
            }
        else:
            return {
                "bam": f"results/{TECH}/{wildcards.ref}/align/phased/non-trio/longphase/{wildcards.sample}/{wildcards.sample}_haplotagged_sorted-linked.bam",
                "bai": f"results/{TECH}/{wildcards.ref}/align/phased/non-trio/longphase/{wildcards.sample}/{wildcards.sample}_haplotagged_sorted-linked.bam.bai",
            }

def get_methyl_targets():
    methyl_files = []
    for row in manifest_df.itertuples():
        if pd.notnull(row.maternal_illumina_fofn) and pd.notnull(row.paternal_illumina_fofn):
            phase_type = "trio"
            if TECH == "ont":
                methyl_files.extend(
                    [
                        f"results/{TECH}/{row.reference_name}/methylation/phased/{phase_type}/{row.sample}/{row.sample}_{hap}_cpg-pileup.{ext}"
                        for hap in HAPS for ext in ["bed.gz", "bw"]
                    ],
                )
                methyl_files.extend(
                    [f"results/{TECH}/{row.reference_name}/methylation/phased/{phase_type}/{row.sample}/{row.sample}_cpg-pileup.{ext}" for ext in ["bed.gz", "bw"]]
                )
            elif TECH == "hifi":
                methyl_files.extend(
                    [
                        f"results/{TECH}/{row.reference_name}/methylation/phased/{phase_type}/{row.sample}/{row.sample}_{hap}_cpg-pileup.combined.{ext}"
                        for hap in HAPS for ext in ["bed.gz", "bw"]
                    ]
                )
                methyl_files.extend(
                    [
                        f"results/{TECH}/{row.reference_name}/methylation/phased/{phase_type}/{row.sample}/{row.sample}_cpg-pileup.combined.{ext}" for ext in ["bed.gz", "bw"]
                    ]
                )
        else:
            phase_type = "non-trio"
            if TECH == "ont":
                methyl_files.extend(
                    [
                        f"results/{TECH}/{row.reference_name}/methylation/phased/{phase_type}/{row.sample}/{row.sample}_{hap}_cpg-pileup.{ext}"
                        for hap in HAPS+["unknown"] for ext in ["bed.gz", "bw"]
                    ],
                )

                methyl_files.extend(
                    [
                        f"results/{TECH}/{row.reference_name}/methylation/phased/{phase_type}/{row.sample}/{row.sample}_cpg-pileup.{ext}"
                        for ext in ["bed.gz", "bw"]
                    ]

                )
            elif TECH == "hifi":
                methyl_files.extend(
                    [
                        f"results/{TECH}/{row.reference_name}/methylation/phased/{phase_type}/{row.sample}/{row.sample}_cpg-pileup.combined.{ext}"
                        for ext in ["bed.gz", "bw"]
                    ]
                )

    return methyl_files

def get_pipeline_resources(which_one, caller):
    def inner(wildcards):
        modify_reference_name = wildcards.ref.split("-")[0]
        if which_one == "trf":
            try:
                trf_file = config["resources"][modify_reference_name]["trf"]
            except KeyError:
                trf_file = None
            if caller == "sniffles" and config["sniffles"]["trf"]:
                if trf_file:
                    return f"--tandem-repeats {trf_file}"
                else:
                    return ""
            return trf_file
        else:
            sys.exit(f"Unsupported param in get_pipeline_resources: {which_one}")

    return inner


def get_reference(wildcards):
    try:
        reference_path = config["reference"][wildcards.ref]
    except KeyError:
        reference_path = manifest_df.at[(wildcards.sample, wildcards.ref), "reference_path"]

    return reference_path


def get_reference_fai(wildcards):
    reference = get_reference(wildcards)
    return reference + ".fai"

def get_assembly_trio_inputs(which_one):
    def inner(wildcards):
        if which_one == "fastq":
            return manifest_df.loc[manifest_df["sample"] == wildcards.sample, "fofn"][0]
        elif which_one == "parental_illumina":
            return {
                "pat": manifest_df.loc[manifest_df["sample"] == wildcards.sample, "paternal_illumina_fofn"][0],
                "mat": manifest_df.loc[manifest_df["sample"] == wildcards.sample, "maternal_illumina_fofn"][0]
            }
        else:
            sys.exit(f"Unsupported param in get_assembly_trio_inputs: {which_one}")
    return inner


def get_yak_input(wildcards):
    yak_df = manifest_df.copy()
    yak_df = yak_df.reset_index(drop=True).set_index(["group_name"], drop=False)
    mat_illumina = yak_df.loc[yak_df.at[wildcards.family, "maternal_illumina_fofn"].notnull(), "maternal_illumina_fofn"][0]
    pat_illumina = yak_df.loc[yak_df.at[wildcards.family,  "paternal_illumina_fofn"].notnull(), "paternal_illumina_fofn"][0]
    if wildcards.parental == "mat":
        return mat_illumina
    elif wildcards.parental == "pat":
        return pat_illumina
    else:
        sys.exit(f"Unsupported param: {wildcards.parental}")


def get_parental_yak(wildcards):
    family = manifest_df.query(fr"sample == '{wildcards.sample}'")["group_name"][0]
    return {
        "mat": f"results/hifi/yak/parents/{family}/mat.yak",
        "pat": f"results/hifi/yak/parents/{family}/pat.yak",
    }


def get_tech_files(which_one):

    def inner(wildcards):

        fp = manifest_df.loc[(manifest_df["sample"] == wildcards.sample) & (manifest_df["reference_name"] == wildcards.ref), which_one][0]
        fofn_df = pd.read_table(fp,header=None,names=["file_path"])
        fofn_df["basename"] = fofn_df.apply(
            lambda row: os.path.basename(row.file_path),axis=1
        )
        fofn_df["n"] = fofn_df.index
        if which_one == "fofn":

            fofn_df["cell"] = fofn_df.basename.str.extract(r"(.*)\.fastq.*")

            # If no duplicates are present in file_path but exists in cell, then add unique identifier to cell
            if fofn_df["cell"].duplicated().any() and not fofn_df["file_path"].duplicated().any():
                fofn_df["cell"] = fofn_df.apply(lambda x: f"{x.cell}-{x.n}", axis=1)

            fofn_df.set_index(["cell"],inplace=True)
            return fofn_df.at[wildcards.cell, "file_path"]

        elif which_one == "unmapped_bam_fofn":

            fofn_df["cell"] = fofn_df.basename.str.extract(r"(.*)\.bam.*")

            # If no duplicates are present in file_path but exists in cell, then add unique identifier to cell
            if fofn_df["cell"].duplicated().any() and not fofn_df["file_path"].duplicated().any():
                fofn_df["cell"] = fofn_df.apply(lambda x: f"{x.cell}-{x.n}", axis=1)

            fofn_df.set_index(["cell"],inplace=True)
            return fofn_df.at[wildcards.cell, "file_path"]

    return inner

def gather_tech_bams(which_one):
    def inner(wildcards):
        fp = manifest_df.loc[(manifest_df["sample"] == wildcards.sample) & (manifest_df["reference_name"] == wildcards.ref), "fofn"][0]
        fofn_df = pd.read_table(fp,header=None,names=["file_path"])

        fofn_df["basename"] = fofn_df.apply(
            lambda row: os.path.basename(row.file_path),axis=1
        )
        fofn_df["n"] = fofn_df.index

        fofn_df["cell"] = fofn_df.basename.str.extract(r"(.*)\.fastq.*")

        # If no duplicates are present in file_path but exists in cell, then add unique identifier to cell
        if fofn_df["cell"].duplicated().any() and not fofn_df["file_path"].duplicated().any():
            fofn_df["cell"] = fofn_df.apply(lambda x: f"{x.cell}-{x.n}",axis=1)

        if which_one == "hap_specific":
            return expand("results/{tech}/{{ref}}/align/phased/{{phase_type}}/{{sample}}/{{sample}}_{cell}_{{hap}}_sorted-linked.bam",cell=fofn_df.cell.tolist(), tech=[TECH])
        elif which_one == "not_hap_specific":
            return expand("results/{tech}/{{ref}}/align/phased/{{phase_type}}/minimap2/{{sample}}/{{sample}}_{cell}_sorted-linked.bam",cell=fofn_df.cell.tolist(), tech=[TECH])

    return inner


def get_by_chrom(wildcards):
    reference_fai = get_reference_fai(wildcards)
    fai_df = pd.read_table(
        reference_fai,
        header=None,
        names=["name", "length", "offset", "linebases", "linewidth"],
    )["name"]

    if wildcards.caller == "sniffles":
        return expand(
            "results/{tech}/{{ref}}/variant_call/sniffles/{{sample}}/{chr}/{{sample}}_sniffles.vcf.gz",
            chr=fai_df.tolist(),
            tech=[TECH],
        )
    elif wildcards.caller == "clair3":
        return expand(
            "results/{tech}/{{ref}}/variant_call/clair3/{{sample}}/{chr}/merge_output.vcf.gz",
            chr=fai_df.tolist(),
            tech=[TECH],
        )

def get_final_dss_targets():

    subset_df = dss_df[dss_df["group_name"].isin(manifest_df["group_name"])]

    reference_name_list =  set([x.split("-")[0] for x in dss_df.reference_name.unique()])
    # chrom_set = config["dm_chrom"].split(",")

    group_name_list = subset_df["group_name"].unique().tolist()
    group_a_list = subset_df.loc[dss_df["group"] == "A", "sample"]

    autosomes_set = ['chr{}'.format(x) for x in list(range(1, 23))]

    targets = []

    for ref in reference_name_list:
        current_ref_targets = []

        autosomes = expand(
            "results/{tech}/analysis/methylation/{{ref}}/dss/{group_name}/intersection/{groupA}/{{chrom}}_two_group_DMR-annotated-summary.tsv.gz",
            # chrom = autosomes_set,
            zip,
            group_name = group_name_list,
            groupA = group_a_list,
            tech=[TECH] * len(group_name_list),
        )

        # Perform chrX analysis on female trio samples only
        subset_df = subset_df.loc[~subset_df.maternal_illumina_fofn.isnull()].query("sex == 'F'")
        if not subset_df.empty:
            group_name_list = subset_df["group_name"].tolist()
            sample_list = subset_df["sample"].tolist()

            sex = expand(
                "results/{tech}/analysis/methylation/{{ref}}/dss/{group_name}/no-intersection/{sample}_chrX_two_group_DMR-annotated-summary.tsv.gz",
                zip,
                group_name=group_name_list,
                sample=sample_list,
                tech=[TECH] * len(group_name_list),
            )

            current_ref_targets.extend(sex)

        current_ref_targets.extend(autosomes)

        targets.extend([x.format(ref=ref, chrom=c) for x in current_ref_targets for c in autosomes_set])

    return targets
