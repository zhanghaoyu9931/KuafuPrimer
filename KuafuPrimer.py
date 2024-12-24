import os
import pandas as pd
import argparse
import shutil
import multiprocessing
import datetime

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from tqdm import tqdm
from KuafuPrimerFuncs import *

__author__ = "Haoyu Zhang"
__copyright__ = "Copyright (c) 2024 Zhulab"
__license__ = "The MIT License (MIT)"
__version__ = "1.0"


def blast_to_silvaRef(input_fna, output_txt):
    # [deplicated] use the abundance table as the input now.
    def mothur_class(input_fna):
        os.system(f"bash classify_seqs.sh {input_fna}")
        taxa_level_dict = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
        meta_res = input_fna.replace(".fasta", ".nr_v138_1.wang.taxonomy")
        df_metagenome = pd.read_csv(meta_res, sep="\t", header=None)
        df_metagenome.columns = ["reads_nm", "taxonomic_info"]

        # accuracy of different taxa level
        for i, taxa_ln in enumerate(taxa_level_dict):
            df_metagenome[taxa_ln] = df_metagenome["taxonomic_info"].apply(
                lambda x: x.split(";")[i].split("(")[0]
                if (len(x.split(";")) > i)
                else ""
            )

        df_metagenome = pd.DataFrame(
            df_metagenome, columns=["reads_nm"] + taxa_level_dict
        )
        rm_files = input_fna.replace(".fasta", ".nr_v138_1.wang*")
        os.system(f"rm {rm_files}")
        return df_metagenome

    df_metagenome = mothur_class(input_fna)
    # 1107: drop those seqs not belong to Bacteria
    df_metagenome = df_metagenome[df_metagenome["Kingdom"] == "Bacteria"].reset_index(
        drop=True
    )

    # return the genus table
    srr_id = os.path.basename(output_txt).split("_")[0]
    df_metagenome["Genus"] = df_metagenome["Genus"].apply(
        lambda x: "unclassified" if "unclassified" in x else x
    )
    all_reads = df_metagenome.shape[0]
    df_metagenome = pd.DataFrame(df_metagenome["Genus"].value_counts())
    # 根据pandas版本不同，这里可能会报错：Genus改成count即可
    df_metagenome.rename(columns={"count": srr_id}, inplace=True)
    df_metagenome["Genus"] = list(df_metagenome.index)
    # reset index
    df_metagenome.reset_index(drop=True, inplace=True)
    df_metagenome = df_metagenome[df_metagenome[srr_id] >= 10].reset_index(drop=True)
    # get the relative abundance
    df_metagenome[srr_id] = df_metagenome[srr_id] / all_reads
    return df_metagenome


### main func: user interface
def design_pri_vregion(
    vv,
    core_micro,
    output,
    num_every_spe,
    primer_lens_list,
    taxo_df,
    SILVA_set_pred_16sDeepSeg,
    extend_bp_num,
    rand_seed,
    rm_tmp_files,
    step_search,
):
    # fastPrimer to design primers
    if ("v9" in vv) or ("v1" in vv):
        mismatch_cutoff = 0.75  # release for v9
    else:
        mismatch_cutoff = 0.9  # this threshold could be high
        
    design_primer(
        microbiota_target="designed_primers",
        core_microbiota=core_micro,
        target_vs=vv,
        res_root=output,
        num_every_spe=num_every_spe,
        representative_seqs_pick_method="random",
        extend_bp_num=extend_bp_num,
        deletion_cutoff=0.99,  # parameters for primer design
        mismatch_cutoff=mismatch_cutoff,
        primer_lens_list=primer_lens_list,
        taxo_df=taxo_df,
        SILVA_set_pred_16sDeepSeg=SILVA_set_pred_16sDeepSeg,
        rand_seed=rand_seed,
        rm_tmp_files=rm_tmp_files,
        step_search=step_search,
    )


def primer_design_main(
    input_csv,
    output,
    vs_list=["v3v4"],
    extend_bp_num=50,
    num_every_spe=50,
    primer_lens_list=[20],
    NGS_mode="Single_end",
    input_type="metagenomic",
    taxo_df=None,
    SILVA_set_pred_16sDeepSeg=None,
    abun_cutoff=1e-5,
    thread_num=10,
    rand_seed=None,
    rm_tmp_files=True,
    step_search=1,
):
    ### get the genus profile
    # output dir
    os.makedirs(output, exist_ok=True)
    df_input = pd.read_csv(input_csv)

    df_samples_genus_table = None
    sample_abun_tb = os.path.join(output, f"samples_abundanceTab.csv")
    if input_type == "genera_profiling":
        # genus profile as input directly
        shutil.copy(input_csv, sample_abun_tb)

    if os.path.exists(sample_abun_tb):
        df_samples_genus_table = pd.read_csv(sample_abun_tb) 
    else:
        for srr, fna in zip(
            list(df_input["srr_id"]), list(df_input[f"{NGS_mode}_file"])
        ):
            bla_out_txt = os.path.join(output, f"{srr}_taxa.csv")
            df_now = blast_to_silvaRef(fna, bla_out_txt)
            if df_samples_genus_table is None:
                df_samples_genus_table = df_now
            else:
                df_samples_genus_table = pd.merge(
                    df_samples_genus_table, df_now, on="Genus", how="outer"
                )
        df_samples_genus_table.to_csv(sample_abun_tb, index=False)

    # clear the genus table
    df_samples_genus_table = df_samples_genus_table[
        df_samples_genus_table.iloc[:, 1:].sum(axis=1) > abun_cutoff
    ]
    df_samples_genus_table["Genus"] = df_samples_genus_table["Genus"].apply(
        lambda x: ("_".join(x.split())).lower()
    )
    df_samples_genus_table = df_samples_genus_table[
        df_samples_genus_table["Genus"].isin(list(taxo_df["genus_name"]))
    ]  # genu name should be in the database
    df_samples_genus_table.reset_index(inplace=True, drop=True)
    df_samples_genus_table.to_csv(
        sample_abun_tb.replace("abundanceTab", "abundanceTab_clean"), index=False
    )

    ### get the core microbiota for primer design
    core_micro = list(df_samples_genus_table["Genus"])
    pool = multiprocessing.Pool(processes=thread_num)  # parallel to accelerate
    for vv in tqdm(vs_list):
        pool.apply_async(
            design_pri_vregion,
            (
                vv,
                core_micro,
                output,
                num_every_spe,
                primer_lens_list,
                taxo_df,
                SILVA_set_pred_16sDeepSeg,
                extend_bp_num,
                rand_seed,
                rm_tmp_files,
                step_search,
            ),
        )
    # end of this procedure
    pool.close()
    pool.join()


def parse_args():
    parser = argparse.ArgumentParser(
        description="Design best primers for specific microbial community."
    )
    parser.add_argument("--input", type=str, help="The input tsv file.")
    parser.add_argument("--out_root", type=str, help="The output root.")
    parser.add_argument(
        "--target_vs", type=str, default="v3v4;v4v4", help="The target V-regions."
    )
    parser.add_argument(
        "--NGS_mode",
        type=str,
        default="Single_end",
        help="The metagenomic data type (pair-end or single-end).",
    )
    parser.add_argument(
        "--input_type",
        type=str,
        default="metagenomic",
        help="Input data type of KuafuPrimer (metagenomic data or relevant genera profiling).",
    )

    # parameters for trim corresponding variable regions using 16sDeepSeg.
    parser.add_argument(
        "--num_every_spe",
        type=float,
        default=5,
        help="Number of selected seqs for every genus.",
    )
    parser.add_argument(
        "--extend_bp_num",
        type=int,
        default=20,
        help="Number of extensive bp that will be included to design primer.",
    )
    parser.add_argument(
        "--step_search",
        type=int,
        default=1,
        help="Number of step size for designing primer.",
    )

    # parameters for primer lens range.
    parser.add_argument(
        "--primer_lens_range", type=str, default="20,21", help="Length of primer."
    )

    # 20231107: parameters for reference databases
    parser.add_argument(
        "--SILVA_set_pred_16sDeepSeg",
        type=str,
        default="Model_data/Silva_ref_data/SILVA_set_segmentation_DeepAnno16.csv", # the right ID version
        help="Reference dataset 1.",
    )
    parser.add_argument(
        "--ref_meta_df",
        type=str,
        default="Model_data/Silva_ref_data/silva_145137_EcK12_8f_1492r_TaxaMeta.csv",
        help="Reference dataset 2.",
    )
    parser.add_argument(
        "--ref_id_seqs",
        type=str,
        default="Model_data/Silva_ref_data/silva_145137_EcK12_8f_1492r_IdSeqs.csv",
        help="Reference dataset 3.",
    )

    # other parameters
    parser.add_argument("--thread", type=int, default=10, help="Parallel running.")
    parser.add_argument(
        "--rand_seed",
        type=int,
        default=None,
        help="Random seed to choose silva seqs to do alignment.",
    )
    parser.add_argument(
        "--rm_tmp_files",
        type=bool,
        default=True,
        help="Remove the template files generated by the script.",
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    input_csv = args.input
    output = args.out_root
    target_vregions = args.target_vs.split(";")  # target_vregions
    NGS_mode = args.NGS_mode
    extend_bp_num = args.extend_bp_num
    num_every_spe = args.num_every_spe
    input_type = args.input_type
    step_search = args.step_search

    lens_min, lens_max = args.primer_lens_range.split(",")
    primer_lens_list = [i for i in range(int(lens_min), int(lens_max))]
    thread_num = args.thread  # multi-thread
    rand_seed = args.rand_seed  # default None
    rm_tmp_files = args.rm_tmp_files # remove the tmp files
    
    ## 20231107: set the dataset ##
    # taxonomic informations
    taxo_df = pd.read_csv(args.ref_meta_df)
    taxo_df["genus_name"] = taxo_df["Genus"].apply(
        lambda x: ("_".join(x.split())).lower()
    )
    id_seqs_df = pd.read_csv(args.ref_id_seqs)
    # read 16sDeepSeg demarcation results
    pred_16sDeepSeg_df = pd.read_csv(args.SILVA_set_pred_16sDeepSeg)
    pred_16sDeepSeg_df.rename(
        columns={"id": "silva_id"}, inplace=True
    )  # use correct ID version
    # merge dfs
    SILVA_set_pred_16sDeepSeg = pd.merge(
        id_seqs_df, pred_16sDeepSeg_df, on="silva_id", how="inner"
    )
    print(f"Total dataset for primer design has {SILVA_set_pred_16sDeepSeg.shape[0]} seqs.")

    ## running the main func ##
    start_time = datetime.datetime.now()
    primer_design_main(
        input_csv,
        output,
        target_vregions,
        extend_bp_num=extend_bp_num,
        num_every_spe=num_every_spe,
        primer_lens_list=primer_lens_list,
        NGS_mode=NGS_mode,
        input_type=input_type,
        taxo_df=taxo_df,
        SILVA_set_pred_16sDeepSeg=SILVA_set_pred_16sDeepSeg,
        thread_num=thread_num,
        rand_seed=rand_seed,
        rm_tmp_files=rm_tmp_files,
        step_search=step_search, # candidate primer search step
    )
    os.system("rm ./*mothur*")
    end_time = datetime.datetime.now()

    with open(f"{output}/log_file.txt", "w") as log_file:
        print(
            f"KuafuPrimer designing procedure running time: {end_time - start_time}",
            file=log_file,
        )
