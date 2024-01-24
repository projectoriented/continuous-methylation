import sys
import os

import pandas as pd

# --------  Load sample sheet -------- #
manifest_df = pd.read_table(config["manifest"], dtype=str, header=0)
manifest_df.set_index(["sample", "reference_name"], drop=False, inplace=True)

# -------- Wildcard constraints -------- #
wildcard_constraints:
    tech=TECH,
    sample="|".join(manifest_df.index.get_level_values("sample")),
    parental="pat|mat",
    hap="|".join(HAPS),
    phase_type="trio|non-trio",
    analysis_type="two_group|model_based",
    dss_category="DMR|DML",
    group_name="|".join(manifest_df["group_name"].unique())

# -------- Helper functions -------- #
def get_final_output(wildcards):

    final_output = []

    final_output.extend(get_methyl_targets())
    final_output.extend(cleanup_things())

    return final_output

def get_methyl_targets():
    methyl_files = []
    for row in manifest_df.itertuples():
        if pd.notnull(row.maternal_illumina_fofn) and pd.notnull(row.paternal_illumina_fofn):
            phase_type = "trio"
        else:
            phase_type = "non-trio"

        if TECH == "ont":
            methyl_files.extend(
                [
                    f"results/{TECH}/{row.reference_name}/methylation/phased/{phase_type}/{row.sample}/{row.sample}_cpg-pileup_{suffix}.{ext}"
                    for suffix in ["combined", "unknown"] + HAPS for ext in ["bed.gz", "bw"]
                ],
            )
        elif TECH == "hifi":
            methyl_files.extend(
                [
                    f"results/{TECH}/{row.reference_name}/methylation/phased/{phase_type}/{row.sample}/{row.sample}_cpg-pileup.{suffix}.{ext}"
                    for suffix in ["combined"] + HAPS for ext in ["bed.gz", "bw"]
                ]
            )
        else:
            raise ValueError(f"Unaccounted for {row}")

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

def calc_mem_gb(wildcards, input, attempt, threads):
    mb = max(1.5 * input.size_mb, 1000)
    gb = int(mb / 1000)

    if threads != 1:
        gb = int(max(gb / threads, 2))

    return gb * attempt

def get_reference_fai(wildcards):
    reference = get_reference(wildcards)
    return reference + ".fai"

def get_assembly_trio_inputs(which_one):
    def inner(wildcards):
        if which_one == "fastq":
            return manifest_df.loc[manifest_df["sample"] == wildcards.sample, "fofn"].values[0]
        elif which_one == "parental_illumina":
            return {
                "pat": manifest_df.loc[manifest_df["sample"] == wildcards.sample, "paternal_illumina_fofn"].values[0],
                "mat": manifest_df.loc[manifest_df["sample"] == wildcards.sample, "maternal_illumina_fofn"].values[0]
            }
        else:
            sys.exit(f"Unsupported param in get_assembly_trio_inputs: {which_one}")
    return inner


def get_yak_input(wildcards):
    yak_df = manifest_df.copy()
    yak_df = yak_df.loc[
        yak_df["maternal_illumina_fofn"].notnull(),
        ["group_name", "maternal_illumina_fofn","paternal_illumina_fofn"]
    ].drop_duplicates()
    yak_df = yak_df.reset_index(drop=True).set_index(["group_name"])
    mat_illumina = yak_df.at[wildcards.family, "maternal_illumina_fofn"]
    pat_illumina = yak_df.at[wildcards.family,  "paternal_illumina_fofn"]
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

def get_all_cell_names():
    all_cell_names = []
    sample_list = manifest_df["sample"].tolist()
    ref_name = manifest_df["reference_name"].tolist()

    for x, r in zip(sample_list,ref_name):
        all_cell_names.extend(get_fofn_df(which_one="fofn",sample_name=x,ref_name=ref_name).index.tolist())

    return all_cell_names

def get_fofn_df(which_one, sample_name, ref_name):
    fp = manifest_df.loc[
        (manifest_df["sample"] == sample_name) & (manifest_df["reference_name"] == ref_name), which_one].values[0]
    fofn_df = pd.read_table(fp,header=None,names=["file_path"])
    fofn_df["basename"] = fofn_df.apply(
        lambda row: os.path.basename(row.file_path),axis=1
    )
    fofn_df["n"] = fofn_df.index
    if which_one == "fofn":
        fofn_df["cell"] = fofn_df.basename.str.extract(r"(.*)\.fastq.*")
    elif which_one == "unmapped_bam_fofn":
        fofn_df["cell"] = fofn_df.basename.str.extract(r"(.*)\.bam.*")
    else:
        raise ValueError(f"Invalid argument for {which_one}")

    # If no duplicates are present in file_path but exists in cell, then add unique identifier to cell
    if fofn_df["cell"].duplicated().any() and not fofn_df["file_path"].duplicated().any():
        fofn_df["cell"] = fofn_df.apply(lambda x: f"{x.cell}-{x.n}",axis=1)

    fofn_df.set_index(["cell"],inplace=True)
    return fofn_df

def get_5mC_bams(wildcards):
    if wildcards.phase_type == "trio":
        phase_type = "trio"
        aligner = "minimap2"
    else:
        phase_type = "non-trio"
        aligner = "longphase"

    return {
        "bam": f"results/{TECH}/{wildcards.ref}/align/phased/{phase_type}/{aligner}/{wildcards.sample}/{wildcards.sample}_sorted-5mC-haplotagged.bam",
        "csi": f"results/{TECH}/{wildcards.ref}/align/phased/{phase_type}/{aligner}/{wildcards.sample}/{wildcards.sample}_sorted-5mC-haplotagged.bam.csi"
    }

def get_unmapped_bam(wildcards):
    fofn_df = get_fofn_df(which_one="unmapped_bam_fofn",sample_name=wildcards.sample,ref_name=wildcards.ref)

    return fofn_df.at[wildcards.cell, "file_path"]

def get_fastq(which_one):
    def inner(wildcards):
        fofn_df = get_fofn_df(which_one="fofn",sample_name=wildcards.sample,ref_name=wildcards.ref)
        if which_one == "correspond_reads":
            return fofn_df.at[wildcards.cell, "file_path"]
        else:
            if wildcards.phase_type == "trio":
                binned_type = wildcards.suffix.split("_")[0]
                return f"results/{{tech}}/{{ref}}/align/phased/trio/{{sample}}/fastq/{{sample}}_{binned_type}_{{cell}}.fastq.gz"
            else:
                return fofn_df.at[wildcards.cell, "file_path"]
    return inner

def gather_tech_bams(wildcards):
    # Just get cell/movie names- so either fofn or unmapped_bam_fofn works
    fofn_df = get_fofn_df(which_one="unmapped_bam_fofn",sample_name=wildcards.sample,ref_name=wildcards.ref)

    if wildcards.phase_type == "trio":
        actual_suffix = ["hap1_sorted-5mC-haplotagged", "hap2_sorted-5mC-haplotagged", "non-binnable_sorted-5mC"]
    else:
        actual_suffix = "sorted-5mC"

    return expand(
        "results/{tech}/{{ref}}/align/phased/{{phase_type}}/minimap2/{{sample}}/{{sample}}_{cell}_{actual_suffix}.bam",
        cell=fofn_df.index.tolist(), tech=[TECH],
        actual_suffix=actual_suffix
    )

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

def cleanup_things():
    cleanup = []

    vc_pattern = "results/{tech}/{ref}/variant_call/clair3/{sample}/.cleaned.txt"
    for row in manifest_df.itertuples():
        if not (pd.notnull(row.maternal_illumina_fofn) and pd.notnull(row.paternal_illumina_fofn)):
            cleanup.append(
                vc_pattern.format(
                    tech=TECH,
                    ref=row.reference_name,
                    sample=row.sample
                )
            )
    return cleanup