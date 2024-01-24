# --------  Switches -------- #
methylation_analysis = config.get("differential_methylation", None)

if methylation_analysis:
    dss_df = pd.read_table(config["dss_manifest"], dtype=str, header=0)
    dss_df = dss_df.merge(
        manifest_df.reset_index(drop=True),
        how="inner"
    ).set_index(["sample"],drop=False)

def get_final_dss_targets():

    subset_df = dss_df[dss_df["group_name"].isin(manifest_df["group_name"])]

    reference_name_list =  set([x.split("-")[0] for x in dss_df.reference_name.unique()])
    # chrom_set = config["dm_chrom"].split(",")

    group_name_list = subset_df["group_name"].unique().tolist()
    group_a_list = subset_df.loc[dss_df["group"] == "A", "sample"].tolist()

    # Sort & make sure they are same length.
    group_name_list.sort()
    group_a_list.sort()
    assert len(group_name_list) == len(group_a_list), "dss: group_name and groupA has differing lengths."

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