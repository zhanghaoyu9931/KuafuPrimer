import os
import random
import pandas as pd
import numpy as np
from pandas.errors import EmptyDataError

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import re
import argparse

import multiprocessing
import time
import datetime
from common_utils import *

## global vars
flip_bp_add = 20  # 20bp extension at both ends of the amplicon
NGS_platform_error = (
    {  # from: "Sequencing error profiles of Illumina sequencing instruments"
        "MiSeq": [0.473, 0.983],
        "MiniSeq": [0.613, 0.459],
        "NextSeq_500": [0.429, 0.827],
        "NextSeq_550": [0.593, 0.435],
        "HiSeq_2500": [0.112, 0.544],
        "NovaSeq_6000": [0.109, 0.350],
        "HiSeq_X_Ten": [0.087, 0.126],
    }
)


## utils
def parse_blastTxt_getAmplicon(
    blast_txt="/data3/hyzhang/ont/16s_RNA_seg/res/blast_res/query.txt",
    fr_pri_length={"forward_pri": 20, "reverse_pri": 20},
    K=3,
    fr_pri_atgc={},
    silva_ref_atgc=None,
):
    # Read blast res file and get the in-silico PCR amplicons #
    # K: the permitted mismatch number, however need to consider the degenerate base detaily
    with open(blast_txt, "r") as f:
        lines = f.readlines()
    index_line = lines[3]
    index_line = index_line.split(": ")[-1]
    index_line = index_line.split(",")
    index_line = [x.strip(" ") for x in index_line]
    index_line = ["_".join(x.split()) for x in index_line]

    lines = [x.strip("\n").split("\t") for x in lines if not x.startswith("#")]

    blast_df = pd.DataFrame(lines)
    blast_df.columns = index_line
    blast_df = blast_df.apply(pd.to_numeric, errors="ignore")

    # sucessful amplification must meet the following criteria:
    # 1. no gap opens
    blast_df = blast_df[blast_df["gap_opens"] < 1].reset_index(drop=True)
    # blast_df = (
    #     blast_df.groupby("query_acc.ver").head(100).reset_index(drop=True)
    # )  # fast mode for debug
    # 2. alignment length > 5
    blast_df = blast_df[
        blast_df[["query_acc.ver", "alignment_length"]].apply(
            lambda x: x["alignment_length"] > 5, axis=1
        )
    ]  # fr_pri_length[x["query_acc.ver"]] - 5
    # 3. get the primer and template sequence (full-length of primer)
    silva_ref_atgc.index = silva_ref_atgc["silva_id"]
    blast_df["pri_seq"] = blast_df.apply(
        lambda x: fr_pri_atgc[x["query_acc.ver"]], axis=1
    )  # full-legth of primer

    def get_ref_atgc(x):
        full_ref = silva_ref_atgc.loc[x["subject_acc.ver"], "16s_rna"]
        start_, end_ = x["s._start"], x["s._end"]
        start_q, end_q = x["q._start"], x["q._end"]
        if end_ > start_:
            start_ = start_ - (start_q - 1)
            end_ = end_ + (fr_pri_length["forward_pri"] - end_q)
            return full_ref[start_ - 1 : end_]
        else:
            start_ = start_ + (start_q - 1)
            end_ = end_ - (fr_pri_length["reverse_pri"] - end_q)
            ref_seq_ = full_ref[start_ - 1 : end_ - 2 : -1]
            ref_seq_ = "".join(
                [
                    (atgc_to_complement[degenerate_base_table[s][0]])
                    for s in list(ref_seq_)
                ]
            )
            return ref_seq_
    # blast_df = blast_df[
    #     (
    #         (blast_df["query_acc.ver"] == "forward")
    #         & (blast_df["s._end"] > blast_df["s._start"])
    #     )
    #     | (
    #         (blast_df["query_acc.ver"] == "reverse")
    #         & (blast_df["s._end"] < blast_df["s._start"])
    #     )
    # ].reset_index(drop=True) # [maybe open] forward primer needs to be in the forward strand, and reverse primer needs to be in the reverse strand
    blast_df["ref_seq"] = blast_df.apply(lambda x: get_ref_atgc(x), axis=1)
    # 4. get the in-silico PCR amplicons (including length, Tm and binding probability)
    blast_df = blast_df.loc[
        blast_df.groupby(["query_acc.ver", "subject_acc.ver"])["evalue"].idxmin(),
    ]  # de-replicated

    def get_amplicon_info(x):
        x = x.reset_index(drop=True)
        amplicon_info = {
            "silva_id": x.loc[0, "subject_acc.ver"],
            "pcr_start": -1,
            "pcr_end": -1,
            "amplicon_size": -1,
            "f_bind_p": 0,
            "r_bind_p": 0,
            "f_bind_Tm": 0,
            "r_bind_Tm": 0,
        }
        if ("forward_pri" not in list(x["query_acc.ver"])) or (
            "reverse_pri" not in list(x["query_acc.ver"])
        ):
            # successful amplification needs both forward and reverse primers matched
            return pd.DataFrame([amplicon_info])
        x.index = x["query_acc.ver"]
        amplicon_info["pcr_start"], amplicon_info["pcr_end"] = (
            x.loc["forward_pri", "s._start"],
            x.loc["reverse_pri", "s._start"],
        )
        amplicon_info["amplicon_size"] = (
            amplicon_info["pcr_end"] - amplicon_info["pcr_start"]
        )
        # primer_binding
        f_bind_p, r_bind_p = primer_binding_probability(
            x.loc["forward_pri", "pri_seq"], x.loc["forward_pri", "ref_seq"], K=K
        ), primer_binding_probability(
            x.loc["reverse_pri", "pri_seq"], x.loc["reverse_pri", "ref_seq"], K=K
        )
        f_Tm, r_Tm = Tm_GCcal(
            x.loc["forward_pri", "pri_seq"], x.loc["forward_pri", "ref_seq"]
        ), Tm_GCcal(x.loc["reverse_pri", "pri_seq"], x.loc["reverse_pri", "ref_seq"])
        amplicon_info["f_bind_p"] = f_bind_p
        amplicon_info["r_bind_p"] = r_bind_p

        amplicon_info["f_bind_Tm"] = f_Tm
        amplicon_info["r_bind_Tm"] = r_Tm
        return pd.DataFrame([amplicon_info])

    blast_df = (
        blast_df.groupby("subject_acc.ver")
        .apply(get_amplicon_info)
        .reset_index(drop=True)
    )
    # For debug
    # print('AAAAAAAAAAAAAA: ', blast_df.head())
    return blast_df


def parse_pcr_ali_res(fna_aligned):
    # parse the alignment result
    res_df = []
    recs = SeqIO.parse(fna_aligned, "fasta")
    for rec in recs:
        id_ = rec.id
        seq = str(rec.seq).split("|")[0]
        num_pat = r"\d+"
        num_find = re.findall(num_pat, seq)
        num_find = [int(x) for x in num_find]

        res_now = {
            "silva_id": id_,
            "pcr_start": num_find[0],
            "pcr_end": num_find[2],
            "amplicon_size": num_find[2] - num_find[0],
        }
        res_df.append(res_now)

    if len(res_df) == 0:
        # not found any hit
        res_now = {
            "silva_id": "aaa",
            "pcr_start": -1,
            "pcr_end": -1,
            "amplicon_size": -1,
        }
        res_df.append(res_now)
    res_df = pd.DataFrame(res_df)
    return res_df


def simulate_sequencing_errors(sequence, NGS_platform="", PE_lens=300):
    # func to simulate the sequencing errors
    error_mean, error_std = NGS_platform_error[NGS_platform]
    error_rate = (
        np.random.normal(error_mean, error_std) / 100
    )
    error_sequence = list(sequence)
    overlap = [len(sequence) - PE_lens, PE_lens]

    # error_rate_all_pos = np.arr
    for i in range(len(sequence)):
        current_base = sequence[i]
        current_base = degenerate_base_table[current_base][0]
        if (i > overlap[0]) and (i < overlap[1]):
            # error rate in overlap region is squared
            error_rate_pos = error_rate**2
        else:
            error_rate_pos = error_rate
        if random.random() < error_rate_pos:
            possible_bases = ["A", "T", "G", "C"]
            possible_bases.remove(current_base)
            error_sequence[i] = random.choice(possible_bases) # randomly choose a base
        else:
            error_sequence[i] = current_base

    return "".join(error_sequence)


def write_fasta(df_to_fasta, sequencing_error="no"):
    # write the fasta file
    recs = []
    for i in range(df_to_fasta.shape[0]):
        p_start, p_end, seq_full, f_bind_p, r_bind_p = (
            int(df_to_fasta.loc[i, "pcr_start"]),
            int(df_to_fasta.loc[i, "pcr_end"]),
            df_to_fasta.loc[i, "16s_rna"],
            float(df_to_fasta.loc[i, "f_bind_p"]),
            float(df_to_fasta.loc[i, "r_bind_p"]),
        )
        if p_start < 0:
            # not be successfully amplified
            continue
        if random.uniform(0, 1) > (f_bind_p * r_bind_p):
            # sucessful amplification according to the binding probability
            continue

        p_start = max(p_start - flip_bp_add, 0)
        p_end = min(p_end + flip_bp_add, len(seq_full))
        seq_amplicon = seq_full[p_start:p_end]
        if sequencing_error != "no":
            # simulate the sequencing errors
            seq_amplicon = simulate_sequencing_errors(seq_amplicon, sequencing_error)
        rec = SeqRecord(
            seq=Seq(seq_amplicon),
            id=df_to_fasta.loc[i, "silva_id"],
            description="",
        )
        recs.append(rec)
    return recs


def otu_and_class(fa, minsize, randid):
    uni_fa = f"./uniques_{randid}.fa"
    # find unique read sequences and abundances
    cmd1 = f"usearch -fastx_uniques {fa} -sizeout -relabel Uniq -fastaout {uni_fa}"

    # make 97% OTUs and filter chimeras
    otu_fa = f"./otu_{randid}.fa"
    zotu_fa = f"./zotu_{randid}.fa"
    cmd2 = (
        f"usearch -cluster_otus {uni_fa} -otus {otu_fa} -relabel Otu -minsize {minsize}"
    )
    cmd3 = f"usearch -unoise3 {uni_fa} -zotus {zotu_fa}"

    # classify every amplicon
    silva_template = "./Model_data/Silva_ref_data/silva.nr_v138_1.align"
    taxo_ref = (
        "./Model_data/Silva_ref_data/silva.nr_v138_1.tax"
    )

    cmd4 = f'mothur "#classify.seqs(fasta={fa},reference={silva_template},taxonomy={taxo_ref},cutoff=50,processors=4)"'
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)
    os.system(cmd4)

    # cnt otu and zotu number
    otu_recs = SeqIO.parse(otu_fa, "fasta")
    otu_n = len(list(otu_recs))
    zotu_recs = SeqIO.parse(zotu_fa, "fasta")
    zotu_n = len(list(zotu_recs))

    return otu_n, zotu_n


def pcr_match(ref_fa, primer_txt, output, K=3):
    # primer_match_cmd = f"pcr_match -i {ref_fa} -P {primer_txt} -o {output} -M 1000 -r -K {K} -5 7 -a -v"
    primer_match_cmd = (
        f"pcr_match -i {ref_fa} -P {primer_txt} -o {output} -M 2000 -r -K {K} -3 3 -a"
    )
    os.system(primer_match_cmd)


def after_pcr(pcr_df, ref_df_selectedGenus, rand_label, sequencing_error="MiSeq"):
    # merge with ref metainfo
    res_of_this_pri_df = pd.merge(
        ref_df_selectedGenus,  # ref_metadf,
        pcr_df,
        how="left",  # merge to all seqs
        left_on="silva_id",
        right_on="silva_id",
    )
    for pcr_ind in ["pcr_start", "pcr_end", "amplicon_size"]:
        res_of_this_pri_df[pcr_ind].fillna(-1, inplace=True)
    for pcr_ind in ["f_bind_p", "r_bind_p", "f_bind_Tm", "r_bind_Tm"]:
        res_of_this_pri_df[pcr_ind].fillna(0, inplace=True)

    # simulate amplicon, cluster otu and classify
    amplicon_recs = write_fasta(res_of_this_pri_df, sequencing_error=sequencing_error)
    if len(amplicon_recs) < 10:
        print("Amplicon_recs have small size!")
        res_of_this_pri_df.to_csv(f"./res_of_this_pri_df_{rand_label}.csv", index=False)
        return None, 0.0, 0.0

    SeqIO.write(amplicon_recs, f"./amplicon_{rand_label}.fasta", "fasta")
    otu_n, zotu_n = otu_and_class(
        f"./amplicon_{rand_label}.fasta", minsize=1, randid=rand_label
    )

    # drop seqs cols
    res_of_this_pri_df.drop(
        [
            "taxa_info_summary",
            # "16s_rna",  # deplete this col in meta info
        ],
        axis=1,
        inplace=True,
    )

    # parse prediction res
    pred_taxa_df = f"amplicon_{rand_label}.nr_v138_1.wang.taxonomy"
    pred_taxa_df = pd.read_csv(pred_taxa_df, sep="\t", header=None)
    pred_taxa_df.columns = ["silva_id", "taxainfo_pred"]
    for i, tax_ in enumerate(taxa_level):
        pred_taxa_df[tax_ + "_pred"] = pred_taxa_df["taxainfo_pred"].apply(
            lambda x: x.split(";")[i].split("(")[0] if len(x.split(";")) > i else ""
        )

    res_of_this_pri_df = pd.merge(
        res_of_this_pri_df, pred_taxa_df, on="silva_id", how="left"
    )
    return res_of_this_pri_df, otu_n, zotu_n


def primer_classify_acc(pri_insilico_df):
    pri_insilico_df.fillna("", inplace=True)

    # sequence level accuracy
    compare_index_all = {}
    all_n = pri_insilico_df.shape[0]
    pri_amp_eff_all = pri_insilico_df[
        (pri_insilico_df["amplicon_size"] > 10)
        & (pri_insilico_df["f_bind_p"] * pri_insilico_df["r_bind_p"] > 0.5)
    ].shape[
        0
    ]
    pri_amp_eff_all /= all_n
    compare_index_all["pri_amp_eff"] = pri_amp_eff_all
    compare_index_all["pri_amp_eff_Tm"] = (
        pri_insilico_df[
            (pri_insilico_df["amplicon_size"] > 10)
            & (pri_insilico_df["f_bind_Tm"] >= 55)
            & (pri_insilico_df["r_bind_Tm"] >= 55)
        ].shape[0]
        / all_n
    )  # meet the Tm requirement
    # 1.2: the forward and reverse primer binding probability
    compare_index_all["pri_amp_eff_forward"] = pri_insilico_df["f_bind_p"].mean()
    compare_index_all["pri_amp_eff_reverse"] = pri_insilico_df["r_bind_p"].mean()
    for tax in taxa_level[:6]:  # drop species level
        pri_withTruelabel = pri_insilico_df[
            pri_insilico_df[tax].apply(lambda x: len(x) > 1)
        ]
        pri_withTruelabel.reset_index(inplace=True, drop=True)

        # clear the tax names
        pri_withTruelabel[tax] = pri_withTruelabel[tax].apply(
            lambda x: ("_".join(x.split())).lower()
        )
        pri_withTruelabel[tax + "_pred"] = (
            pri_withTruelabel[tax + "_pred"]
            .fillna("-")
            .apply(lambda x: ("_".join(x.split())).lower())
        )

        n_with_trueLabel = pri_withTruelabel.shape[0]
        n_predTrue = pri_withTruelabel[
            pri_withTruelabel[tax] == pri_withTruelabel[tax + "_pred"]
        ].shape[0]

        compare_index_all[tax + "_accuracy"] = n_predTrue / n_with_trueLabel

    return compare_index_all


def check_dimer_hairpin(
    f_atgc, r_atgc, DimerAndHairpin_exam_root="./DimerAndHairpin_exam"
):
    # dimer和hairpin检查
    os.makedirs(DimerAndHairpin_exam_root, exist_ok=True)
    randlabel = random.randint(1, 999999)
    while os.path.exists(f"{DimerAndHairpin_exam_root}/primer_pair{randlabel}.fa"):
        # avoid the conflict
        randlabel = random.randint(1, 999999)
    dimer_flag, hairpin_flag = False, False

    f_pri = SeqRecord(Seq(f_atgc), id="forward_primer", description="forward")
    r_pri = SeqRecord(Seq(r_atgc), id="reverse_primer", description="reverse")
    pri_fa = [f_pri, r_pri]
    SeqIO.write(
        pri_fa, f"{DimerAndHairpin_exam_root}/primer_pair{randlabel}.fa", "fasta"
    )

    os.system(
        f"mfeprimer dimer -i {DimerAndHairpin_exam_root}/primer_pair{randlabel}.fa -o {DimerAndHairpin_exam_root}/dimer{randlabel}.txt -t 35"
    )
    os.system(
        f"mfeprimer hairpin -i {DimerAndHairpin_exam_root}/primer_pair{randlabel}.fa -o {DimerAndHairpin_exam_root}/hairpin{randlabel}.txt -t 35"
    )

    with open(
        f"{DimerAndHairpin_exam_root}/dimer{randlabel}.txt", "r", encoding="utf-8"
    ) as f:
        dimer_res = f.readlines()
    if "No dimer found.\n" not in dimer_res:
        dimer_flag = True

    with open(
        f"{DimerAndHairpin_exam_root}/hairpin{randlabel}.txt", "r", encoding="utf-8"
    ) as f:
        hairpin_res = f.readlines()
    if "No hairpins found.\n" not in hairpin_res:
        hairpin_flag = True

    os.system(f"rm {DimerAndHairpin_exam_root}/*{randlabel}*")
    return dimer_flag, hairpin_flag


### main func: user interface
def main_func(
    ref_fa,
    ref_id_seq,
    pri_pair,
    genus_select,
    K=3,
    very_fast_mode=False,
    amplicon_method="pcr_match",
    sequencing_error="MiSeq",
    rand_seed=99419,
    num_every_spe=50,
    envi_forEva_selectedID=None # select the silva id to be used
):
    rand_label = random.randint(9931, 99419)
    while os.path.exists(f"./testPrimer{rand_label}.txt"):
        rand_label = random.randint(9931, 99419)

    # write the primer txt
    primer_txt = f"./testPrimer{rand_label}.txt"
    with open(primer_txt, "w") as f:
        if amplicon_method == "pcr_match":
            f.write("\n".join(pri_pair))
        elif amplicon_method == "haoyu":
            f.write(">forward_pri\n" + pri_pair[0] + "\n")
            f.write(">reverse_pri\n" + pri_pair[1] + "\n")

    # get the target genus profile
    ref_df_selectedGenus = ref_metadf[
        ref_metadf["genus_name"].isin(genus_select)
    ]
    if envi_forEva_selectedID is not None:
        # select the silva id to be used
        ref_df_selectedGenus = ref_df_selectedGenus[ref_df_selectedGenus['silva_id'].isin(envi_forEva_selectedID)]

    def rand_samp(group, gp_num_fast_mode=num_every_spe):
        # Key: the global parameter couble be modifed
        return group.sample(
            n=min(gp_num_fast_mode, len(group)), random_state=rand_seed
        )  # set the random seed

    if very_fast_mode:
        ref_df_selectedGenus = pd.DataFrame(
            ref_df_selectedGenus.groupby("Genus").apply(rand_samp)
        )
        ref_df_selectedGenus.reset_index(inplace=True, drop=True)
    ref_df_selectedGenus.reset_index(inplace=True, drop=True)

    # accelerate the process
    if amplicon_method == "pcr_match":
        # run the pcr_match
        output = f"temp_res{rand_label}.txt"
        pcr_match(ref_fa=ref_fa, primer_txt=primer_txt, K=K, output=output)

        # after pcr
        align_df = parse_pcr_ali_res(output)
    elif amplicon_method == "haoyu":
        # use blastn to get the amplicon
        print(ref_df_selectedGenus.head())
        selectedGenus_recs = [
            SeqRecord(
                seq=Seq(ref_df_selectedGenus.loc[i, "16s_rna"]),
                id=ref_df_selectedGenus.loc[i, "silva_id"],
                description="",
            )
            for i in range(ref_df_selectedGenus.shape[0])
        ]
        selectedGenus_fa = os.path.join(blast_db_temp_dir, f"{rand_label}.fa")
        SeqIO.write(selectedGenus_recs, selectedGenus_fa, "fasta")
        selectedGenus_fa_ncbidb = selectedGenus_fa.replace(".fa", "_bladb")
        os.system(
            f"makeblastdb -in {selectedGenus_fa} -dbtype nucl -out {selectedGenus_fa_ncbidb}"
        )
        output = f"temp_res_blast{rand_label}.txt"
        blast_cmd(primer_txt, db=selectedGenus_fa_ncbidb, out_path=output)
        
        # read the blast res
        try:
            blast_df = parse_blastTxt_getAmplicon(
                blast_txt=output,
                fr_pri_length={
                    "forward_pri": len(pri_pair[0]),
                    "reverse_pri": len(pri_pair[1]),
                },
                K=K,
                fr_pri_atgc={
                    k: v for k, v in zip(["forward_pri", "reverse_pri"], pri_pair)
                },
                silva_ref_atgc=ref_id_seq,
            )
            align_df = blast_df.copy()
        except:
            os.system(f"rm ./*{rand_label}* ./*mothur* {selectedGenus_fa_ncbidb}* {selectedGenus_fa}")
            return None, 0, 0, None
    if (
        align_df[
            (align_df["amplicon_size"] > 20)
            & (align_df["f_bind_p"] * align_df["r_bind_p"] > 1e-4)
        ].shape[0]
        < 3
    ):
        # too few amplicons
        os.system(f"rm ./*{rand_label}* ./*mothur* {selectedGenus_fa_ncbidb}* {selectedGenus_fa}")
        return None, 0, 0, None

    try:
        res_of_this_pri_df, otu_n, zotu_n = after_pcr(
            align_df,
            ref_df_selectedGenus,
            rand_label=rand_label,
            sequencing_error=sequencing_error,
        )
        acc_index = primer_classify_acc(res_of_this_pri_df)
    except:
        res_of_this_pri_df, otu_n, zotu_n, acc_index = None, 0, 0, None

    # rm files
    os.system(f"rm ./*{rand_label}* ./*mothur* {selectedGenus_fa_ncbidb}* {selectedGenus_fa}")
    return res_of_this_pri_df, otu_n, zotu_n, acc_index


def get_designed_primer(
    designed_pri_dir="./demo_dir", v_region="v3v4", temp_dir="", K=1
):
    # get the desigend primer pairs
    envi_nm = os.path.basename(designed_pri_dir)
    vs_dir = f"{designed_pri_dir}/designed_primers"
    primer_pairs = []
    primer_pairs_nm = []
    if not os.path.exists(f"{vs_dir}/{v_region}"):
        return primer_pairs, primer_pairs_nm

    try:
        forward_df = pd.read_csv(f"{vs_dir}/{v_region}/forward_primer.csv")
        forward_df.sort_values(by=f"cov_num_k_{K}", ascending=False, inplace=True)
        forward_df.reset_index(inplace=True, drop=True)

        reverse_df = pd.read_csv(f"{vs_dir}/{v_region}/reverse_primer.csv")
        reverse_df.sort_values(by=f"cov_num_k_{K}", ascending=False, inplace=True)
        reverse_df.reset_index(inplace=True, drop=True)

        DimerAndHairpin_exam_root = os.path.join(temp_dir, "DimerAndHairpin_exam")
        os.makedirs(DimerAndHairpin_exam_root, exist_ok=True)

        # make primer pairs between top 10 forward and reverse primers (100 candidates at most)
        fp_list = list(forward_df["primer_atgc"])[:10]
        rp_list = list(reverse_df["primer_atgc"])[:10]

        pp_loc_list = []
        # [0, 0] [1, 0] [0, 1] [1, 1] is the order
        for i in range(min(len(fp_list), len(rp_list))):
            for j in range(i + 1):
                for k in range(i + 1):
                    if [j, k] not in pp_loc_list:
                        pp_loc_list.append([j, k])
        if len(fp_list) < len(rp_list):
            for i in range(len(fp_list), len(rp_list)):
                for j in range(len(fp_list)):
                    pp_loc_list.append([j, i])
        elif len(fp_list) > len(rp_list):
            for i in range(len(rp_list), len(fp_list)):
                for j in range(len(rp_list)):
                    pp_loc_list.append([i, j])

        for pp_loc in pp_loc_list:
            f_i = pp_loc[0]
            r_i = pp_loc[1]
            f_atgc = fp_list[f_i]
            r_atgc = rp_list[r_i]
            pri_atgc = [f_atgc, r_atgc]
            pri_nm = f"{envi_nm}_{v_region}_{f_i}f{r_i}r"
            f_Tm, r_Tm = forward_df.loc[f_i, "Tm"], reverse_df.loc[r_i, "Tm"]
            if abs(f_Tm - r_Tm) > 5:
                # Tm difference should be less than 5
                continue

            # check dimer and hairpin
            dimer_flag, hairpin_flag = check_dimer_hairpin(f_atgc, r_atgc)
            if dimer_flag or hairpin_flag:
                continue

            primer_pairs.append(pri_atgc)
            primer_pairs_nm.append(pri_nm)

        # rm temp files
        os.system(f"rm -r {DimerAndHairpin_exam_root}")
        return primer_pairs, primer_pairs_nm

    # except EmptyDataError:
    except:
        return primer_pairs, primer_pairs_nm


def get_additional_primers(
    pris_xlsx="data/primers.xlsx", targetVregion=""
):
    if pris_xlsx.endswith(".csv"):
        pris_df = pd.read_csv(pris_xlsx)
    else:
        pris_df = pd.read_excel(pris_xlsx)
    # add universal primers
    primer_pairs = []
    primer_pairs_nm = []

    # get primers of the target V-region
    pris_df = pris_df[(pris_df["v_region"] == targetVregion)]
    pris_df.reset_index(inplace=True, drop=True)
    for pp_i in range(pris_df.shape[0]):
        pripair_atgc = [
            pris_df.loc[pp_i, "forward_seq"],
            pris_df.loc[pp_i, "reverse_seq"],
        ]
        p_nm = pris_df.loc[pp_i, "pri_nm"]
        if targetVregion not in p_nm:
            p_nm = p_nm + "_" + targetVregion

        primer_pairs.append(pripair_atgc)
        primer_pairs_nm.append(p_nm)

    print(f"Get {len(primer_pairs)} additional primers for {targetVregion} V-region")
    return primer_pairs, primer_pairs_nm


# parse the args
def parse_args():
    parser = argparse.ArgumentParser(
        description="Design best primers for specific environment."
    )
    parser.add_argument(
        "--envi_forEva",
        type=str,
        help="The microbial community in which the primer performance is to be evaluated.",
    )
    parser.add_argument( 
        "--envi_forEva_selectedID",
        type=str,
        default="",
        help="Set silva ids that are allowed to eva.",
    )
    parser.add_argument(
        "--primers_forEva", type=str, help="The designed primers to be evaluated."
    )
    parser.add_argument("--amplicon_method", default="haoyu", type=str, help=".")
    parser.add_argument(
        "--sequencing_error",
        default="no",  # NextSeq_550
        type=str,
        help="Platform on which to simulate the sequencing error.",
    )
    parser.add_argument("--K", type=int, help="The permitted mismatch numbers.")
    parser.add_argument(
        "--target_vs", type=str, default="v3v4;v4v4", help="The target v regions."
    )
    parser.add_argument(
        "--output",
        type=str,
        default="./output/demo_output_InSilicoPCR",
        help="The output dir.",
    )
    parser.add_argument(
        "--additional_primers",
        type=str,
        default="Model_data/Universal_primers/primers.xlsx",
        help="User can add their own primers for the comparison.",
    )
    parser.add_argument(
        "--num_condidate_pris",
        type=int,
        default=9,
        help="How many candidate primer pairs from a V-region to compare.",
    )
    parser.add_argument("--very_fast", action="store_true", default=False, help=".")
    parser.add_argument(
        "--rand_seed", type=int, default=99419, help="rand_seed for choosing seqs."
    )
    parser.add_argument(
        "--num_every_spe", type=int, default=20, help="num for choosing seqs."
    )
    # parameters to set the reference database
    parser.add_argument(
        "--ref_fa",
        type=str,
        default="Model_data/Silva_ref_data/silva_seqs_145137_EcK12_8f_1492r.fasta",
        help=".",
    )
    parser.add_argument(
        "--ref_meta_df",
        type=str,
        default="Model_data/Silva_ref_data/silva_145137_EcK12_8f_1492r_TaxaMeta.csv",
        help=".",
    )
    parser.add_argument(
        "--ref_id_seqs",
        type=str,
        default="Model_data/Silva_ref_data/silva_145137_EcK12_8f_1492r_IdSeqs.csv",
        help=".",
    )

    parser.add_argument("--thread", type=int, default=10, help="Parallel running.")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    primers_forEva = [f"{x}" for x in args.primers_forEva.split(";")]  # multiple primer sets to be evaluated
    envi_forEva = args.envi_forEva
    envi_forEva_selectedID = args.envi_forEva_selectedID
    if envi_forEva_selectedID.endswith('.txt'):
        with open(envi_forEva_selectedID, 'r') as f:
            envi_forEva_selectedID = list(f.readlines())
            envi_forEva_selectedID = [x.strip('\n') for x in envi_forEva_selectedID]
    else:
        envi_forEva_selectedID = None # default None
     
    targetVregions = args.target_vs.split(
        ";"
    )
    amplicon_method = args.amplicon_method
    sequencing_error = args.sequencing_error
    K = args.K
    num_condidate_pris = int(args.num_condidate_pris) 
    additional_pri_excel = (
        args.additional_primers
    )
    thread_num = args.thread
    # reference database parameters
    ref_fa = args.ref_fa
    ref_metadf = pd.read_csv(args.ref_meta_df)
    ref_idseqs = pd.read_csv(args.ref_id_seqs)
    ref_idseqs = pd.DataFrame(
        ref_idseqs,
        columns=[
            "silva_id",
            "16s_rna",
        ],
    )
    ref_metadf = pd.merge(ref_metadf, ref_idseqs, on="silva_id", how="inner")
    ref_metadf["genus_name"] = ref_metadf["Genus"].apply(
        lambda x: ("_".join(x.split())).lower()
    )
    print(f"Reference database Seqs num: {ref_metadf.shape}")
    # very fast runing mode
    very_fast_mode = args.very_fast
    rand_seed = args.rand_seed  # rand_seed
    num_every_spe = args.num_every_spe
    if very_fast_mode:
        print(f"Running in fast mode at the expense of some accuracy.")

    # directory and func for in-silico PCR
    root_dir = f"{args.output}_K{K}"
    os.makedirs(root_dir, exist_ok=True)

    def single_pri_cov_eva(
        pri_pair,
        pri_nm,
        ref_fa,
        genus_select,
        K,
        very_fast_mode,
        temp_dir,
        sequencing_error="",
        rand_seed=99419,
        num_every_spe=50,
        envi_forEva_selectedID=None,
    ):
        primer_pos = get_PP_position_Ecoli_K12(pri_pair[0], pri_pair[1])

        if os.path.exists(f"{temp_dir}/{pri_nm}_res.csv"):
            # avoid the repeat running
            acc_index = {"pri_amp_eff": 0.01}
            for tax_l in taxa_level:
                acc_index[f"{tax_l}_accuracy"] = 0.01
            acc_index["otu_n"] = 0
            acc_index["zotu_n"] = 0
            acc_index["pri_nm"] = pri_nm
            acc_index["pri_pair"] = pri_pair
            acc_index.update(primer_pos)
            return acc_index

        res_of_this_pri_df, otu_n, zotu_n, acc_index = main_func(
            ref_fa=ref_fa,
            ref_id_seq=ref_idseqs,
            pri_pair=pri_pair,
            genus_select=genus_select,
            K=K,
            very_fast_mode=very_fast_mode,
            amplicon_method=amplicon_method,
            sequencing_error=sequencing_error,
            rand_seed=rand_seed,
            num_every_spe=num_every_spe,
            envi_forEva_selectedID=envi_forEva_selectedID,
        )
        if res_of_this_pri_df is None:
            # default index
            acc_index = {"pri_amp_eff": 0.01}
            for tax_l in taxa_level:
                acc_index[f"{tax_l}_accuracy"] = 0.01
            acc_index["otu_n"] = 0
            acc_index["zotu_n"] = 0
            acc_index["pri_nm"] = pri_nm
            acc_index["pri_pair"] = pri_pair 
            acc_index.update(primer_pos)
            return acc_index

        res_of_this_pri_df.to_csv(f"{temp_dir}/{pri_nm}_res.csv", index=False)
        acc_index["otu_n"] = otu_n
        acc_index["zotu_n"] = zotu_n
        acc_index["pri_nm"] = pri_nm
        acc_index["pri_pair"] = pri_pair
        acc_index.update(primer_pos)
        return acc_index

    def pri_cov_eva(
        targetVregion="v3v4",
    ):
        # evaluate the primer performance of a specific region
        temp_dir = os.path.join(root_dir, targetVregion)
        # if os.path.exists(f"{temp_dir}/resInfo.csv"):
        #     return
        os.makedirs(temp_dir, exist_ok=True)
        print(f"Start to evaluate the primer performance of {targetVregion} region.")

        ## primer list to be evaluated
        primer_pairs = []
        primer_pairs_nm = []

        for primer_forEva in primers_forEva:
            primer_pairs_t, primer_pairs_nm_t = get_designed_primer(
                primer_forEva, targetVregion, temp_dir, K=K 
            )

            primer_pairs += primer_pairs_t[
                :num_condidate_pris
            ] 
            primer_pairs_nm += primer_pairs_nm_t[:num_condidate_pris]

        primer_pairs_t, primer_pairs_nm_t = get_additional_primers(
            pris_xlsx=additional_pri_excel, targetVregion=targetVregion
        )
        primer_pairs += primer_pairs_t
        primer_pairs_nm += primer_pairs_nm_t

        if len(primer_pairs) == 0:
            return None

        ## get the target genus list
        genus_select = pd.read_csv(envi_forEva)  # input is a csv file
        genus_select["Genus"] = genus_select["Genus"].apply(
            lambda x: ("_".join(x.split())).lower()
        )
        genus_select.to_csv(
            f"{root_dir}/genus_profiling_for_evaluation.csv", index=False
        )
        genus_select = list(genus_select["Genus"])

        ## 开始验证
        # avoid the repeat running
        primer_nm_atgc_no_overwrite = {pri_nm: pri_pair for pri_pair, pri_nm in zip(primer_pairs, primer_pairs_nm) if (not os.path.exists(f"{temp_dir}/{pri_nm}_res.csv"))}
        primer_pairs, primer_pairs_nm = list(primer_nm_atgc_no_overwrite.values()), list(primer_nm_atgc_no_overwrite.keys())
        
        pool = multiprocessing.Pool(processes=thread_num)  # parallel to accelerate
        pris_df = []
        for pri_pair, pri_nm in zip(primer_pairs, primer_pairs_nm):
            pris_df.append(
                pool.apply_async(
                    single_pri_cov_eva,
                    (
                        pri_pair,
                        pri_nm,
                        ref_fa,
                        genus_select,
                        K,
                        very_fast_mode,
                        temp_dir,
                        sequencing_error,
                        rand_seed,
                        num_every_spe,
                        envi_forEva_selectedID,
                    ),
                )
            )
        # end of this procedure
        pool.close()
        pool.join()

        pris_df = [x.get() for x in pris_df]
        if len(pris_df) == 0:
            return
        pris_df = pd.DataFrame(pris_df)
        if os.path.exists(f"{temp_dir}/resInfo.csv"):
            # merge the results
            pris_df_ori = pd.read_csv(f"{temp_dir}/resInfo.csv") 
            pris_df = pd.concat([pris_df_ori, pris_df], axis=0)
            
        pris_df.to_csv(f"{temp_dir}/resInfo.csv", index=False)

    ## running
    start_time = datetime.datetime.now()
    for targetVregion in targetVregions:
        pri_cov_eva(targetVregion)
        # pool.apply_async(pri_cov_eva, (targetVregion,))
    end_time = datetime.datetime.now()

    with open(f"{root_dir}/log_file.txt", "w") as log_file:
        print(
            f"KuafuPrimer In-silico PCR procedure running time: {end_time - start_time}",
            file=log_file,
        )