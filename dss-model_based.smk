import pandas as pd
import sys

# --------- Globals --------- #
DSS_CNTR="/net/eichler/vol26/7200/software/containers/R/DSS/2.46.0/bioconductor-dss_2.46.0--r42hc0cfd56_0.sif"
AUTOSOMES=['chr{}'.format(x) for x in list(range(1, 23))]
SEX=["chrX"]
HAP=["hap1", "hap2"]
MODEL_ONE_COVARIATES=["case", "familyid", "caseasd", "casemother", "casesibling"]
MODEL_TWO_COVARIATES=["case", "familyid", "allele"]

# --------- Load files --------- #
configfile: "resources.yaml"
configfile: "config.yaml"

df = pd.read_table(config["manifest"], dtype=str, header=0)

# # Make sure proband is always first
# df["sort_order"] = df.apply(lambda row: "a" if "_p" in row["sample"] else "b", axis=1)
# df.sort_values("sort_order", inplace=True)

# --------- Constraints --------- #
wildcard_constraints:
    model="model_one|model_two",
    chr="|".join(AUTOSOMES+SEX),
    sample="|".join(df["sample"]),
    suffix="cpg|hap1_cpg|hap2_cpg"

# --------- Input functions --------- #
def get_final_dss_targets(wildcards):

    model_one = expand(
        # "results/model_one/{chr}_{covariates}_DMR.tsv",
        "results/model_one/intersection/annotation/wgs_DMR-annotated_summary.tsv.gz",
        # covariates=["case", "familyid", "case-coef1", "case-coef2"]
    )

    model_two = expand(
        # "results/model_two/{chr}_{covariates}_DMR.tsv",
        "results/model_two/intersection/annotation/wgs_DMR-annotated_summary.tsv.gz",
        # covariates=["case", "familyid", "allele"]
    )

    return model_one+model_two

def get_dss_inputs(model="model_one"):
    def inner(wildcards):

        if model == "model_one":
            sample_names = df.query("sex == 'F'").groupby("family_id")["sample"].filter(lambda x: x.count() == 3).tolist()
            file_names = [f"results/mCG/{s}_cpg.txt" for s in sample_names]
        elif model == "model_two":
            sample_names = df.loc[~df["sample"].str.contains("mo|fa"), "sample"].tolist()
            file_names = [f"results/mCG/{s}_{h}_cpg.txt" for s in sample_names for h in HAP]
        else:
            raise ValueError("Invalid param in get_dss_inputs")
        return file_names
    return inner

def get_dss_params(model="model_one"):

    def inner(wildcards):
        if model == "model_one":

            case_converter = {
                "_p": "asd",
                "_m": "mother",
                "_s": "sibling"
            }

            sample_names = df.query("sex == 'F'").groupby("family_id")["sample"].filter(lambda x: x.count() == 3).tolist()
            allele_names = [""]
            case_names = [case_converter.get(x[x.index("_"):-1]) for x in sample_names]

        elif model == "model_two":

            case_converter = {
                "_p" : "asd",
                "_s" : "control"
            }

            hap_converter = {
                "hap1": "paternal",
                "hap2": "maternal"
            }

            sample_names = [f"{x}_{h}" for x in df.loc[~df["sample"].str.contains("mo|fa"), "sample"].tolist() for h in HAP]
            allele_names = [hap_converter.get(x[x.index("hap"):]) for x in sample_names]
            case_names = [case_converter.get(x[x.index("_"):-6]) for x in sample_names]
        else:
            pass

        param_dict = {
            "sample_names": sample_names,
            "case": case_names,
            "family_names": [x.split("_")[0] for x in sample_names],
            "allele_names": allele_names,
            "model": model,
        }

        return param_dict
    return inner

def get_dss_groups(wildcards):
    if wildcards.model == "model_one":
        return {
            "a": [f"results/model_one/wgs_case_DMR.tsv"],
            "b": [f"results/model_one/wgs_{c}_DMR.tsv" for c in MODEL_ONE_COVARIATES[1:]]
        }
    elif wildcards.model == "model_two":
        return {
            "a": [f"results/model_two/wgs_allele_DMR.tsv"],
            "b": [f"results/model_two/wgs_{c}_DMR.tsv" for c in MODEL_TWO_COVARIATES if "allele" not in c]
        }

def get_dss_summaries_by_chrom(wildcards):
    targets = []

    file_pattern = "results/{{model}}/intersection/temp/wgs_DMR-{sample}_{suffix}.tsv"

    for row in df.itertuples():
        if "_mo" in row.sample or "_fa" in row.sample:
            targets.append(file_pattern.format(sample=row.sample, suffix="cpg"))
        else:
            for h in HAP:
                hap_suffix = f"{h}_cpg"
                targets.append(file_pattern.format(sample=row.sample,suffix=hap_suffix))
            targets.append(file_pattern.format(sample=row.sample,suffix="cpg"))

    return targets

def get_anno_dict(wildcards):
    return config["annotations"]["GRCh38"]

# --------- Begin --------- #
rule all:
    input: get_final_dss_targets

rule dss_prepare_in_txt:
    input:
        bed="bed-pileups/{sample}_{suffix}-pileup.bed.gz",
    output:
        bed_by_chrom=temp("results/mCG/{sample}_{suffix}.txt")
    threads: config["default"]["threads"]
    resources:
        mem=lambda wildcards, attempt: attempt * config["default"]["mem"],
        hrs=config["default"]["hrs"]
    shell:
        """
        # Columns grabbed are based on this documentation: https://github.com/nanoporetech/modkit/#bedmethyl-column-descriptions

        # chrom, start, n_valid, n_mod
        zcat {input.bed} | awk '{{print $1,$2,$10,$12}}' FS='\\t' OFS='\\t' > {output.bed_by_chrom}
        """

rule dss_model_one:
    input:
        file_names = get_dss_inputs(model="model_one")
    output:
        case_dmr = "results/model_one/wgs_case_DMR.tsv",
        family_dmr = "results/model_one/wgs_familyid_DMR.tsv",
        case_coef1_dmr = "results/model_one/wgs_casesibling_DMR.tsv",
        case_other_dmr = "results/model_one/wgs_casemother_DMR.tsv",
        case_coef2_dmr = "results/model_one/wgs_caseasd_DMR.tsv",
        case_dml = "results/model_one/wgs_case_DML.tsv",
        family_dml = "results/model_one/wgs_familyid_DML.tsv",
        case_coef1_dml = "results/model_one/wgs_casesibling_DML.tsv",
        case_other_dml = "results/model_one/wgs_casemother_DML.tsv",
        case_coef2_dml = "results/model_one/wgs_caseasd_DML.tsv",
    wildcard_constraints:
        chr="chr[0-9]+|chrX",
    params:
        get_dss_params(model="model_one"),
        output_prefix="results/model_one/wgs",
    threads: config["analysis"]["dss"]["threads"]
    resources:
        mem=lambda wildcards, attempt: attempt * config["analysis"]["dss"]["mem"],
        hrs=config["analysis"]["dss"]["hrs"]
    log: "results/model_one/log/wgs.log"
    container:
        DSS_CNTR
    script:
        "scripts/DSS.R"

rule dss_model_two:
    input:
        file_names = get_dss_inputs(model="model_two")
    output:
        case_dmr = "results/model_two/wgs_case_DMR.tsv",
        family_dmr = "results/model_two/wgs_familyid_DMR.tsv",
        allele_dmr = "results/model_two/wgs_allele_DMR.tsv",
        case_dml = "results/model_two/wgs_case_DML.tsv",
        family_dml = "results/model_two/wgs_familyid_DML.tsv",
        allele_dml = "results/model_two/wgs_allele_DML.tsv",
    params:
        get_dss_params(model="model_two"),
        output_prefix="results/model_two/wgs"
    threads: config["analysis"]["dss"]["threads"]
    resources:
        mem=lambda wildcards, attempt: attempt * config["analysis"]["dss"]["mem"],
        hrs=config["analysis"]["dss"]["hrs"]
    log: "results/model_two/log/wgs.log"
    container:
        DSS_CNTR
    script:
        "scripts/DSS.R"


rule dss_intersect:
    input:
        unpack(get_dss_groups)
    output:
        dmr="results/{model}/intersection/wgs_DMR.tsv"
    threads: config["default"]["threads"]
    resources:
        mem=lambda wildcards, attempt: attempt * config["default"]["mem"],
        hrs=config["default"]["hrs"]
    log: "results/{model}/intersection/log/wgs_DMR.log"
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

        intersection_df = ( a_bed + list(b_dict.values()) ).to_dataframe(disable_auto_names=True, names=target_cols,header=None)

        if intersection_df.empty:
            print(f"No intersection happened for {wildcards.model}", file=sys.stderr)
            print(f"Combining both files in output.\nA: {a}\nB: {b}", file=sys.stderr)
            intersection_df = pd.concat([ pd.read_table(x, header=0, dtype={"start": int, "end": int}) for x in [a] + b])
        else:
            print(f"Before intersection {a_df.shape} (row, columns);\nA: {a};\nB: {b};", file=sys.stderr)
            print(f"After intersection {intersection_df.shape} (row, columns);\nA: {a};\nB: {b};", file=sys.stderr)

        intersection_df.to_csv(output.dmr, sep='\t', header=True, index=False)


rule dss_summary_table:
    input:
        dss_out = "results/{model}/intersection/wgs_DMR.tsv",
        sample_bed = "results/mCG/{sample}_{suffix}.txt"
    output:
        sample_summary_table = temp("results/{model}/intersection/temp/wgs_DMR-{sample}_{suffix}.tsv")
    threads: config["default"]["threads"]
    resources:
        mem= lambda wildcards,attempt: attempt * config["default"]["mem"],
        hrs=config["default"]["hrs"]
    script:
        "scripts/calculate_percent_meth.py"


rule merge_dss_summary_table:
    input:
        all_samples=get_dss_summaries_by_chrom,
    output:
        chrom_summary_table=temp("results/{model}/intersection/annotation/wgs_DMR-prcnt_summary.tsv.gz")
    threads: config["default"]["threads"]
    resources:
        mem= lambda wildcards,attempt: attempt * config["default"]["mem"],
        hrs=config["default"]["hrs"]
    script:
        "scripts/merge_dss_summary_tables.py"


rule add_annotation:
    input:
        unpack(get_anno_dict),
        merged_summary = "results/{model}/intersection/annotation/wgs_DMR-prcnt_summary.tsv.gz",
    output:
        added_annotation = "results/{model}/intersection/annotation/wgs_DMR-annotated_summary.tsv.gz"
    threads: config["default"]["threads"]
    resources:
        mem= lambda wildcards,attempt: attempt * config["default"]["mem"],
        hrs=config["default"]["hrs"]
    run:
        from pybedtools import BedTool

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
                closest_genes_df.rename(columns={"distance(bp)": f"distance(bp)_{key}"}, inplace=True)

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
