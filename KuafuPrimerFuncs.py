import random
import os

import numpy as np
import pandas as pd
from sklearn.utils import shuffle
from tqdm import tqdm

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import Counter

from common_utils import *

#### funcs
super_conserve_pos_cutoff = 0.995
possible_conserve_pos_cutoff = 0.98  # TODO: set small for v1 and v9 region
degebase_cutoff = 0.001  # when the base frequency is higher than this, it will be considered as degebase
cov_k_by = 0

# atgc order table
atgc_order = ["A", "T", "G", "C"]


# GC ratio and Tm calculation
def TmGc(primer_atgc="AA"):
    for dege_base, subs_nor_base in degenerate_base_table.items():
        subs_nor_base = subs_nor_base[0]
        primer_atgc = primer_atgc.replace(dege_base, subs_nor_base)
    # 计算GC含量 Tm
    atgc_cnt = dict(Counter(primer_atgc))
    for base in atgc_order:
        if base not in atgc_cnt:
            atgc_cnt[base] = 0
    GC_ratio = (atgc_cnt["G"] + atgc_cnt["C"]) / len(primer_atgc)
    Tm = 4 * (atgc_cnt["G"] + atgc_cnt["C"]) + 2 * (atgc_cnt["A"] + atgc_cnt["T"])
    return GC_ratio, Tm


# using blast to do this
def parse_pcr_ali_res_singlePri(fna_aligned):
    # parse the alignment result
    res_df = []
    recs = SeqIO.parse(fna_aligned, "fasta")
    for rec in recs:
        id_ = rec.id
        res_now = {"silva_id": id_, "pcr_start": -1, "pcr_end": -1}
        res_df.append(res_now)

    if len(res_df) == 0:
        # not found any hit
        res_now = {
            "silva_id": "aaa",
            "pcr_start": -1,
            "pcr_end": -1,
        }
        res_df.append(res_now)
    res_df = pd.DataFrame(res_df) 
    return res_df


def parse_pcr_ali_res_singlePri_haoyu(
    fna_aligned, K=3, pri_atgc="GTGTAGCGGTGAAATGCKTA", ref_fa=None, forward_reverse=""
):
    # read the blast result and get the in-silico amplicon
    recs = SeqIO.parse(ref_fa, "fasta")
    recs = {str(x.id): str(x.seq) for x in recs}

    with open(fna_aligned, "r") as f:
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

    # requirements for successful amplicon
    # 1.no gap opens
    blast_df = blast_df[blast_df["gap_opens"] < 1].reset_index(drop=True)
    # 2.alignment length is enough
    blast_df = blast_df[
        blast_df["alignment_length"].apply(lambda x: x >= len(pri_atgc) - 5)
    ]
    # 3.get the ref and primer sequence
    blast_df["pri_seq"] = blast_df.apply(lambda x: pri_atgc, axis=1)

    def get_ref_atgc(x):
        # TODO: some pri-ref maybe align in reverse direction (need to check)
        if x["subject_acc.ver"] not in recs:
            print(x["subject_acc.ver"])
            return "".join(["A"] * len(pri_atgc))

        full_ref = recs[x["subject_acc.ver"]]
        start_, end_ = x["s._start"], x["s._end"]
        start_q, end_q = x["q._start"], x["q._end"]
        try:
            if end_ > start_:
                start_ = start_ - (start_q - 1)
                end_ = end_ + (len(pri_atgc) - end_q)
                return_refseq = full_ref[start_ - 1 : end_]
                return return_refseq
            else:
                start_ = start_ + (start_q - 1)
                end_ = end_ - (len(pri_atgc) - end_q)
                ref_seq_ = full_ref[start_ - 1 : end_ - 2 : -1]
                ref_seq_ = "".join(
                    [
                        (atgc_to_complement[degenerate_base_table[s][0]])
                        for s in list(ref_seq_)
                    ]
                )
                return ref_seq_
        except:
            # print("bbbbbbb") # for debug
            return ""

    blast_df["ref_seq"] = blast_df.apply(lambda x: get_ref_atgc(x), axis=1)
    # 4.possibility of sucessful amplicon
    blast_df = blast_df.loc[
        blast_df.groupby(["query_acc.ver", "subject_acc.ver"])["evalue"].idxmin(),
    ]  # de-duplicate
    blast_df.reset_index(inplace=True, drop=True)

    def get_amplicon_info(x):
        x = x.reset_index(drop=True)
        amplicon_info = {
            "silva_id": x.loc[0, "subject_acc.ver"],
            "pcr_start": -1,
            "pcr_end": -1,
            "bind_p": 0,
            "Tm_gc_cal": 0,
        }
        # primer_binding
        bind_p = primer_binding_probability(
            x.loc[0, "pri_seq"], x.loc[0, "ref_seq"], K=K
        )
        # Tm calculation
        Tm_gc = Tm_GCcal(x.loc[0, "pri_seq"], x.loc[0, "ref_seq"])
        amplicon_info["bind_p"] = bind_p
        amplicon_info["Tm_gc_cal"] = Tm_gc
        return pd.DataFrame([amplicon_info])

    blast_df = (
        blast_df.groupby("subject_acc.ver")
        .apply(get_amplicon_info)
        .reset_index(drop=True)
    )
    return blast_df


def get_single_pri_match_num(
    ref_fa,
    output,
    pri_atgc="GTGTAGCGGTGAAATGCKTA",
    Ks_list=[0, 1],
    pri_bind_method="primer_match",
    forward_reverse="forward",
):
    return_info = {}
    for K in Ks_list:
        os.system(f"rm {output}")
        if pri_bind_method == "primer_match":
            primer_match_cmd = (
                f"primer_match -i {ref_fa} -p {pri_atgc} -o {output} -r -K {K} -3 3"
            )
            os.system(primer_match_cmd)
            
            res_df = parse_pcr_ali_res_singlePri(output)
            return_info[f"cov_num_k_{K}"] = res_df.shape[0]
        elif pri_bind_method == "haoyu":
            # using blast to do in-silico PCR (by default)
            selectedGenus_fa_ncbidb = ref_fa.replace(".fasta", "_bladb")
            if not os.path.exists(selectedGenus_fa_ncbidb + ".nin"):
                os.system(
                    f"makeblastdb -in {ref_fa} -dbtype nucl -out {selectedGenus_fa_ncbidb}"
                )
            rand_label = random.randint(9931, 99419)
            while os.path.exists(f"./testPrimer{rand_label}.txt"):
                # avoid the same label
                rand_label = random.randint(9931, 99419)
            primer_txt = f"./testPrimer{rand_label}.txt"
            with open(primer_txt, "w") as f:
                f.write(f">{forward_reverse}_pri\n" + pri_atgc + "\n")
            blast_cmd(primer_txt, db=selectedGenus_fa_ncbidb, out_path=output)
            os.system(f"rm {primer_txt}")
            res_df = parse_pcr_ali_res_singlePri_haoyu(
                output,
                K=K,
                pri_atgc=pri_atgc,
                ref_fa=ref_fa,
                forward_reverse=forward_reverse,
            )
            res_df.to_csv(ref_fa.replace(".fasta", "_temp_pri.csv"), index=False)
            try:
                return_info[f"cov_num_k_{K}"] = res_df["bind_p"].sum()
            except:
                return_info[f"cov_num_k_{K}"] = 0
            try:
                res_df["Tm_gc_cal"].fillna(0, inplace=True)
                return_info[f"Tm_cut"] = res_df[res_df["Tm_gc_cal"] >= 55].shape[
                    0
                ]
            except:
                return_info[f"Tm_cut"] = 0

    return return_info


# func to design primer from MSA file
def primer_design_from_MAS(
    mas_fna="designTool/gut/v3v4/forward_conserved_afterMuscle.fasta",
    degebase_cutoff=0.01,  
    deletion_cutoff=0.99,  
    mismatch_cutoff=0.99,  
    primer_lens_list=[20],  # primer length
    forward_reverse="forward",
    step_search=1,  # search step
):
    recs = SeqIO.parse(mas_fna, "fasta")
    seqs_np = []
    for rec in recs:
        seq_now = str(rec.seq)
        for dege_base, subs_nor_base in degenerate_base_table.items():
            subs_nor_base = subs_nor_base[
                0
            ]  # simply replace the degebase to the first candidate base
            seq_now = seq_now.replace(dege_base, subs_nor_base)
        seqs_np.append(list(seq_now))
    seqs_np = np.array(seqs_np)

    # distribution of ATGC at each position
    total_seq_num = seqs_np.shape[0]  # total number of sequences
    pos_conserved = [{"A": -1, "T": -1, "G": -1, "C": -1, "-": -1}]
    for i in range(seqs_np.shape[1]):
        seqs_this_pos = list(seqs_np[:, i])
        if forward_reverse == "reverse":
            # supplementary for reverse primer
            seqs_this_pos = [atgc_to_complement[x] for x in seqs_this_pos]

        atgc_cnt = Counter(seqs_this_pos)
        pos_conserved.append(dict(atgc_cnt))
    pos_conserved = pd.DataFrame(pos_conserved)
    pos_conserved.fillna(0.0, inplace=True)
    pos_conserved = pos_conserved.iloc[1:, :]
    pos_conserved.reset_index(inplace=True, drop=True)

    # get the label of each position (candidate conserved or degebase)
    def pos_conserve_label(atgc_freq=[], degebase_cutoff=0.1, deletion_cutoff=0.98):
        atgc_poss = [x / sum(atgc_freq) for x in atgc_freq]
        # if conserved or candidate degebase
        if atgc_poss[-1] > 0.05:  # candidate (0.01 - 0.001)
            if atgc_poss[-1] > deletion_cutoff:
                return "deletion_couldDrop", "-", max(atgc_poss), "-", max(atgc_poss)
            else:
                return "deletion_soMuch", "-", max(atgc_poss), "-", max(atgc_poss)
        if max(atgc_poss) > super_conserve_pos_cutoff:
            # super conserved
            return (
                "super_conserve_pos",
                atgc_order[atgc_poss.index(max(atgc_poss))],
                max(atgc_poss),
                atgc_order[atgc_poss.index(max(atgc_poss))],
                max(atgc_poss),
            )

        # candidate degenerate base
        obvious_base_type = []
        dege_bases_freqSum = 0.0
        for i, atgc in enumerate(atgc_order):
            if atgc_poss[i] > degebase_cutoff:
                obvious_base_type.append(atgc)
                dege_bases_freqSum += atgc_poss[i]
        obvious_base_type.sort()

        possible_degeBase = ""
        for degeB, degeList in degenerate_base_table.items():
            degeList.sort()
            if obvious_base_type == degeList:
                possible_degeBase = degeB
                break
        if len(possible_degeBase) > 0:
            return (
                "possible_degebase_pos",
                possible_degeBase,
                dege_bases_freqSum,
                atgc_order[
                    atgc_poss.index(max(atgc_poss))
                ],  # the most frequent base except degebase
                max(atgc_poss),
            )
        else:
            return (
                "normal_pos",
                atgc_order[atgc_poss.index(max(atgc_poss))],
                max(atgc_poss),
                atgc_order[atgc_poss.index(max(atgc_poss))],
                max(atgc_poss),
            )

    degeprimer_df = []
    for pos_i in range(pos_conserved.shape[0]):
        atgc_freq = list(pos_conserved.iloc[pos_i, :])
        res = pos_conserve_label(
            atgc_freq, degebase_cutoff=degebase_cutoff, deletion_cutoff=deletion_cutoff
        )
        degeprimer_df.append(
            {
                "pos_type": res[0],
                "pos_base": res[1],
                "pos_base_freq_max": res[2],
                "pos_base_nodege": res[3],
                "pos_base_freq_max_nodege": res[4],
            }
        )

    degeprimer_df = pd.DataFrame(degeprimer_df)
    pos_conserved = pd.concat([pos_conserved, degeprimer_df], axis=1)
    pos_conserved.to_csv(mas_fna.replace("_afterMuscle.fasta", ".csv"), index=False)
    print(f"All position after alignment: {pos_conserved.shape}")

    # find the potential primers
    # drop the deletion_couldDrop bases
    pos_conserved = pos_conserved[
        pos_conserved["pos_type"] != "deletion_couldDrop"
    ]  # drop base of deletion
    print(f"After drop - positon: {pos_conserved.shape}")
    pos_conserved.reset_index(inplace=True, drop=True)

    potential_primers = []
    for primer_lens in primer_lens_list:
        for start_pos in tqdm(range(
            0, pos_conserved.shape[0] + 1 - primer_lens, step_search
        )):  # step_search to speed up
            end_pos = start_pos + primer_lens
            primer_now = pos_conserved.iloc[start_pos:end_pos, :].copy()
            primer_now.reset_index(inplace=True, drop=True)
            ## constraints for appropriate primer 
            # 1. no -
            if "deletion_soMuch" in list(primer_now["pos_type"]):
                continue
            # 2. <= 3 degebase, and no degebase at both ends
            for ii in range(len(primer_now) - 3, len(primer_now)):
                # no degebase at the end
                if primer_now.loc[ii, "pos_type"] == "possible_degebase_pos":
                    # if (
                    #     primer_now.loc[ii, "pos_base_freq_max_nodege"]
                    #     < possible_conserve_pos_cutoff
                    # ):
                    #     continue

                    # switch to normal base
                    primer_now.loc[ii, "pos_type"] = "possible_conserved_pos"
                    primer_now.loc[ii, "pos_base"] = primer_now.loc[
                        ii, "pos_base_nodege"
                    ]
                    primer_now.loc[ii, "pos_base_freq_max"] = primer_now.loc[
                        ii, "pos_base_freq_max_nodege"
                    ]
                    
            for ii in range(0, 3):
                # no degebase at the start
                if primer_now.loc[ii, "pos_type"] == "possible_degebase_pos":
                    # if (
                    #     primer_now.loc[ii, "pos_base_freq_max_nodege"]
                    #     < possible_conserve_pos_cutoff
                    # ):
                    #     continue

                    # switch to normal base
                    primer_now.loc[ii, "pos_type"] = "possible_conserved_pos"
                    primer_now.loc[ii, "pos_base"] = primer_now.loc[
                        ii, "pos_base_nodege"
                    ]
                    primer_now.loc[ii, "pos_base_freq_max"] = primer_now.loc[
                        ii, "pos_base_freq_max_nodege"
                    ]

            pos_type_summary = dict(Counter(list(primer_now["pos_type"])))
            if "possible_degebase_pos" not in pos_type_summary:
                pos_type_summary["possible_degebase_pos"] = 0
            if pos_type_summary["possible_degebase_pos"] > 3:
                dege_pos_rank = list(
                    primer_now[primer_now["pos_type"] == "possible_degebase_pos"]
                    .sort_values(by="pos_base_freq_max_nodege", ascending=False)
                    .index
                )
                dege_pos_to_fix = dege_pos_rank[
                    : pos_type_summary["possible_degebase_pos"] - 3
                ]
                for dege_i in dege_pos_to_fix:
                    primer_now.loc[dege_i, "pos_type"] = "possible_conserved_pos"
                    primer_now.loc[dege_i, "pos_base"] = primer_now.loc[
                        dege_i, "pos_base_nodege"
                    ]
                    primer_now.loc[dege_i, "pos_base_freq_max"] = primer_now.loc[
                        dege_i, "pos_base_freq_max_nodege"
                    ]
                # continue
            # 3. no more than 3 mismatch
            possibility_ls = list(primer_now["pos_base_freq_max"])
            mismatch = [x for x in possibility_ls if x < mismatch_cutoff]
            if len(mismatch) > 3:
                continue
            # 4. Tm and GC ratio constraints
            if forward_reverse == "reverse":
                # reverse primer
                primer_atgc = "".join(list(primer_now["pos_base"])[::-1])
            else:
                primer_atgc = "".join(list(primer_now["pos_base"]))

            GC_ratio, Tm = TmGc(primer_atgc)
            if GC_ratio < 0.4 or GC_ratio > 0.65:
                continue
            if Tm < 55 or Tm > 66:  #  65
                continue
            # 5. get the in-silico PCR of the primer
            pri_cov_num_info = get_single_pri_match_num(
                ref_fa=mas_fna.replace("_afterMuscle.fasta", ".fasta"),
                pri_atgc=primer_atgc,
                output=mas_fna.replace(
                    "_conserved_afterMuscle.fasta", "_pri_cov_temp.txt"
                ),
                Ks_list=[0, 1, 2, 3],
                pri_bind_method="haoyu",
                forward_reverse=forward_reverse,
            )
            for ky in pri_cov_num_info.keys():
                pri_cov_num_info[ky] = pri_cov_num_info[ky] / total_seq_num
            if pri_cov_num_info[f"cov_num_k_{cov_k_by}"] < 0.85:
                # coverage is too low
                continue

            # get the start position of the primer on the reference
            primer_degenum = pos_type_summary["possible_degebase_pos"]
            primer_base_coverage_ave = np.mean(possibility_ls)
            start_Ecoli = get_PP_position_Ecoli_K12(primer_atgc, primer_atgc)

            possible_primer_info_now = {
                "primer_atgc": primer_atgc,
                "mismatch_num": len(mismatch),
                "primer_degenum": primer_degenum,
                "primer_base_coverage_ave": primer_base_coverage_ave,
                "GC_ratio": GC_ratio,
                "Tm": Tm,
                "start_pos_related": start_pos,
                "start_pos_Ecoli": start_Ecoli[f"{forward_reverse}_start"],
            }
            possible_primer_info_now.update(pri_cov_num_info)
            potential_primers.append(possible_primer_info_now)

    potential_primers = pd.DataFrame(potential_primers)
    return potential_primers


# func to create fna file for primer design
def creat_primer_fna(rna_seq, pos, id, des, fna_file):
    ## create a fna to store the sub-regions
    records = []
    lens = []
    # merge
    for i in range(min(len(rna_seq), 50000)):
        spe_id = id[i]

        # create records
        primer_seq = rna_seq[i]
        _from, _to = pos[i]  # [a, b, c, d]

        # Get the variable region
        if _from[1] == -1 or _to[0] == -1:
            print("Have -1.", _from[1], _to[0])
            continue
        primer_seq = primer_seq[_from[1] : _to[0]]

        lens.append(len(primer_seq))
        primer_seq = Seq(primer_seq)

        rec = SeqRecord(primer_seq)
        rec.id = spe_id + f"_{i}"
        rec.description = des[i]

        records.append(rec)

    # remove > 3 std
    len_m, len_std = np.mean(lens), np.std(lens)
    records = [
        x
        for x, l in zip(records, lens)
        if (l >= len_m - 3 * len_std) and (l <= len_m + 3 * len_std)
    ]

    # write to fna file
    print(f"Get {len(records)} records.")

    SeqIO.write(records, fna_file + ".fasta", "fasta")
    lens = [len(r) for r in records]
    return lens


# parse the annotation res
def str_list(str='"[1204, 1224]"'):
    if "notfound" in str or "wrongorder" in str:
        return [0, 0]

    if "[" not in str:
        # from vxtractor
        str = str.split()[0]
        str = str.strip("'")
        str = str.split("-")
    else:
        str = str.strip('"[').strip(']"')
        str = str.split(", ")

    str = [int(float(x)) for x in str]
    return str


# func to remove _idx from fasta file after multi-alignment
def fasta_description_change(fasta_input="v3v3_Aquaculture_Vregion_afterMuscle.fasta"):
    records = SeqIO.parse(fasta_input, "fasta")
    records_after = []
    for i, rec in enumerate(records):
        idx = rec.id
        idx = idx.split("_")[:2]
        idx = "_".join(idx)

        rec.id = idx
        rec.description = ""
        records_after.append(rec)
    SeqIO.write(records_after, fasta_input, "fasta")


#### main funcs
def design_primer(
    microbiota_target="gut",
    core_microbiota=[],
    target_vs="v3v3",
    res_root="/data1/hyzhang/Projects/16sDeepSeg_summary/Evi_specific_primers_database/Results",
    num_every_spe=5000,
    representative_seqs_pick_method="random",
    extend_bp_num=50,
    deletion_cutoff=0.99,  
    mismatch_cutoff=0.99,  # usually set to 0.9
    primer_lens_list=[20],
    taxo_df=None,
    SILVA_set_pred_16sDeepSeg=None,
    rand_seed=None,
    rm_tmp_files=True, # if remove the temp files
    step_search=1, # search step
):
    # num_every_spe could be a ratio
    if num_every_spe < 1:
        ratio_every_spe = num_every_spe
    else:
        ratio_every_spe = None
        num_every_spe = int(num_every_spe)
    
    # get core species
    core_microbiota = [x.strip().lower() for x in core_microbiota]
    print(f"The target microbial community has {len(core_microbiota)} genus.")

    # extract corresponding v-regions and flanking conserved regions
    print("Target " + target_vs)
    target_id = target_vs + "_target"
    target_v_root = os.path.join(f"{res_root}/{microbiota_target}", target_vs)
    os.makedirs(target_v_root, exist_ok=True)
    if os.path.exists(os.path.join(target_v_root, 'reverse_primer.csv')):
        # already have the primer designed, skip to save time
        return 1

    target_vs = target_vs.split("v")[1:]
    target_vs = [int(x) for x in target_vs]

    # func to design forward or reverse primer
    def design_forwardOrReverse(target_v=3, fr="reverse", num_every_spe = 100, ratio_every_spe = None):
        if fr == "forward":
            start_primer = f"v{target_v - 1}"
            end_primer = f"v{target_v}"
        else:
            start_primer = f"v{target_v}"
            end_primer = f"v{target_v + 1}"
        ## get the conserved region for MSA and primer design
        # scale = 1
        # extract regions
        start_pos, end_pos, ids, rna_seqs = [], [], [], []

        for micro_genus in core_microbiota:
            if len(micro_genus.split()) > 1:
                print(f"Maybe unusual genus name: {micro_genus}.")

            # get representative seqs for target genera
            taxo_df_micro_now = taxo_df[taxo_df["genus_name"] == micro_genus]
            if ratio_every_spe is not None:
                num_every_spe = int(ratio_every_spe * taxo_df_micro_now.shape[0])
                
            if representative_seqs_pick_method == "lens":
                # select the longest seqs as representative
                taxo_df_micro_now.sort_values(by="lens", ascending=False, inplace=True)
            else:
                # randomly select
                if rand_seed is not None:
                    taxo_df_micro_now = taxo_df_micro_now.sample(
                        frac=1, random_state=rand_seed
                    )
                else:
                    taxo_df_micro_now = taxo_df_micro_now.sample(frac=1)

            taxo_df_micro_now.reset_index(inplace=True, drop=True)
            taxo_df_micro_now = taxo_df_micro_now.iloc[:num_every_spe, :]  # 本来有+ 1删掉了

            micro_genus_id = list(
                taxo_df_micro_now["silva_id"] # ["silva_id_wrong"]
            )
            
            SILVA_set_pred_16sDeepSeg_t = SILVA_set_pred_16sDeepSeg[
                SILVA_set_pred_16sDeepSeg["silva_id"].isin(micro_genus_id) # ["silva_id_wrong"]
            ]
            SILVA_set_pred_16sDeepSeg_t = shuffle(SILVA_set_pred_16sDeepSeg_t)
            SILVA_set_pred_16sDeepSeg_t.reset_index(inplace=True, drop=True)
            if SILVA_set_pred_16sDeepSeg_t.shape[0] < num_every_spe - 5:
                # randomly sample with replacement
                SILVA_set_pred_16sDeepSeg_t = SILVA_set_pred_16sDeepSeg_t.sample(
                    n=num_every_spe, replace=True
                )

            # not any representative seqs, continue
            if SILVA_set_pred_16sDeepSeg_t.shape[0] < 1:
                continue

            # target region length
            if start_primer not in SILVA_set_pred_16sDeepSeg_t.columns:
                start_pos += [
                    [0, 0]
                    for x in list(SILVA_set_pred_16sDeepSeg_t["v1"])[:num_every_spe]
                ]
            else:
                start_pos += [
                    str_list(x)
                    for x in list(SILVA_set_pred_16sDeepSeg_t[start_primer])[
                        :num_every_spe
                    ]
                ]

            if end_primer not in SILVA_set_pred_16sDeepSeg_t.columns:
                end_pos += [
                    [-2, -2]
                    for x in list(SILVA_set_pred_16sDeepSeg_t["v1"])[:num_every_spe]
                ]
            else:
                end_pos += [
                    str_list(x)
                    for x in list(SILVA_set_pred_16sDeepSeg_t[end_primer])[
                        :num_every_spe
                    ]
                ]

            ids += list(SILVA_set_pred_16sDeepSeg_t["silva_id"])[:num_every_spe] # ["silva_id_wrong"]
            rna_seqs += list(SILVA_set_pred_16sDeepSeg_t["16s_rna"])[:num_every_spe]

        # print(num_every_spe, 'AAAAA')
        with open(os.path.join(target_v_root, "silva_id_used.txt"), "w") as f:
            # save the used silva_id
            f.write("\n".join(ids))

        # extend the conserved region
        extend_bp_num_fr = extend_bp_num
        if start_primer == "v0" or end_primer == "v10":
            # for v1 and v9 region
            extend_bp_num_fr = 2 * extend_bp_num_fr
        pos_info = [
            [
                (st[0], max(st[1] - extend_bp_num_fr, 1)),
                (en[0], en[1]) if (en[0] < 0) else (en[0] + extend_bp_num_fr, en[1]),
            ]
            for st, en in zip(start_pos, end_pos)
        ]
        describe = [""] * len(rna_seqs)
        ids = [target_id] * len(rna_seqs)
        # create fna file
        fna_name = f"{fr}_conserved"
        lens = creat_primer_fna(
            rna_seqs,
            pos=pos_info,
            id=ids,
            des=describe,
            fna_file=os.path.join(target_v_root, fna_name),
        )

        # run muscle to multiple-alignment
        from Bio.Align.Applications import MuscleCommandline

        fna_name_muscle = fna_name
        cline = MuscleCommandline(
            input=os.path.join(target_v_root, fna_name_muscle + ".fasta"),
            out=os.path.join(target_v_root, fna_name_muscle + "_afterMuscle.fasta"),
        )
        with open(f"{target_v_root}/run_muscle.bash", "w") as f:
            # add otu cluster to reduce the time cost [deprecated]
            f.write(
                "## Commandline for clustering seqs to OTUs with usearch tool.\n"
            )  # add comments

            file_in = os.path.join(target_v_root, fna_name_muscle + ".fasta")
            uniq_file = os.path.join(target_v_root, "uniques.fa")
            usearch_cmd = (
                f"usearch -fastx_uniques {file_in} -fastaout {uniq_file} -sizeout -relabel Uniq\n"
                + f'usearch -cluster_otus {uniq_file} -otus {file_in} -relabel {fna_name_muscle.split("_")[0]}_target_\n'
            )
            # f.write(usearch_cmd) # 

            # Do multi-alignment using muscle.
            f.write(
                "\n## Commandline for multi-alignment seqs with muscle tool.\n"
            )  # add comments
            cline = str(cline).replace("-in", "-align").replace("-out", "-output")
            cline = cline.replace("-align", "-super5") + " -threads 32"
            # cline += ' -replicates 5000' # allow N replicates
            f.write(str(cline))

        if len(lens) < 5:
            # too few seqs, so the primer design process failed.
            return None

        # then run the bash file in corresponding dir
        os.system("bash " + f"{target_v_root}/run_muscle.bash")
        # post-process
        post_muscle_file = os.path.join(
            target_v_root, fna_name_muscle + "_afterMuscle.fasta"
        )
        recs = SeqIO.parse(post_muscle_file, "fasta")

        post_recs = []
        for rec in recs:
            rec.description = rec.description.split()[0]
            post_recs.append(rec)
        with open(post_muscle_file, "w") as f:
            SeqIO.write(post_recs, f, "fasta")

        # # remove the idx from final multi-alignment file
        # fasta_description_change(post_muscle_file)
        fr_primer_df = primer_design_from_MAS(
            post_muscle_file,
            degebase_cutoff=degebase_cutoff,
            deletion_cutoff=deletion_cutoff,
            mismatch_cutoff=mismatch_cutoff,
            primer_lens_list=primer_lens_list,
            forward_reverse=fr,
            step_search=step_search
        )

        return fr_primer_df
    
    fg = design_forwardOrReverse(target_vs[0], "forward", num_every_spe = num_every_spe, ratio_every_spe = ratio_every_spe)
    rg = design_forwardOrReverse(target_vs[1], "reverse", num_every_spe = num_every_spe, ratio_every_spe = ratio_every_spe)
    # rm tmp files
    if rm_tmp_files:
        for rm_file_nm in ["conserved", "pri_cov_temp", "run_muscle.bash"]:
            os.system(f"rm {target_v_root}/*{rm_file_nm}*")
    if (fg is None) or (rg is None):
        return -1  # unable to design primer
    else:
        fg.to_csv(os.path.join(target_v_root, "forward_primer.csv"), index=False)
        rg.to_csv(os.path.join(target_v_root, "reverse_primer.csv"), index=False)

if __name__ == "__main__":
    design_primer(
        microbiota_target="gut",
        core_microbiota=["Pseudomonas", "Dickeya", "Prevotella"],
        target_vs="v3v4",
        res_root="designTool",
        num_every_spe=5,
        extend_bp_num=50,
    )
