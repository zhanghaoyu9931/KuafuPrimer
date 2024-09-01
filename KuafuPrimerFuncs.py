## 0409: 全部使用自己的方法来做，不使用Decipher; 针对多序列比对结果，设计forward、backward的引物

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

#### 用到的一些funcs
# 1219：将定义超保守的这些参数提取到此定义
super_conserve_pos_cutoff = 0.995
possible_conserve_pos_cutoff = 0.98  #
degebase_cutoff = 0.001  # 大于多少频率的碱基纳入degebase的考虑范围
cov_k_by = 0

# 一些依据的table
atgc_order = ["A", "T", "G", "C"]


# 计算GC含量、Tm等
def TmGc(primer_atgc="AA"):
    # 把其他dege的base换成一个普通base
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


# 使用primer_match来计算primer能够match上的序列数目
def parse_pcr_ali_res_singlePri(fna_aligned):
    # 这个函数解析match后的结果
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
    res_df = pd.DataFrame(res_df)  # 保留id信息是为了后续需要做genus的匹配做准备
    return res_df


def parse_pcr_ali_res_singlePri_haoyu(
    fna_aligned, K=3, pri_atgc="GTGTAGCGGTGAAATGCKTA", ref_fa=None, forward_reverse=""
):
    # 这个函数解析blast之后的结果
    # 读取blast结果，同时确认是否能够amplicon
    # K: 允许的错配数量，但是不能简单用blastn结果，因为有兼并碱基
    # id没有具体意义，不是silva-id，只是1-n编号
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

    # 进行amplicon与否的几个条件筛选
    # 1.不允许gap
    blast_df = blast_df[blast_df["gap_opens"] < 1].reset_index(drop=True)
    # 2.覆盖引物全长: 1226放松一些，最多有5个少的alignment
    blast_df = blast_df[
        blast_df["alignment_length"].apply(lambda x: x >= len(pri_atgc) - 5)
    ]
    # 3.拿出引物序列和参考序列
    blast_df["pri_seq"] = blast_df.apply(lambda x: pri_atgc, axis=1)  # 无论如何都是引物全长

    def get_ref_atgc(x):
        # 注意现在有可能q的start-end覆盖不住
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
            print("bbbbbbb")
            return ""

    blast_df["ref_seq"] = blast_df.apply(lambda x: get_ref_atgc(x), axis=1)
    # 4.计算amplicon成功的概率
    blast_df = blast_df.loc[
        blast_df.groupby(["query_acc.ver", "subject_acc.ver"])["evalue"].idxmin(),
    ]  # 去除重复
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
        # 0102：加入Tm计算的code
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
            # 读取结果
            res_df = parse_pcr_ali_res_singlePri(output)
            return_info[f"cov_num_k_{K}"] = res_df.shape[0]
        elif pri_bind_method == "haoyu":
            # todo: 修改为使用blast然后加入primer bind的准确率
            selectedGenus_fa_ncbidb = ref_fa.replace(".fasta", "_bladb")
            if not os.path.exists(selectedGenus_fa_ncbidb + ".nin"):
                os.system(
                    f"makeblastdb -in {ref_fa} -dbtype nucl -out {selectedGenus_fa_ncbidb}"
                )
            rand_label = random.randint(9931, 99419)
            while os.path.exists(f"./testPrimer{rand_label}.txt"):
                # 避免出现和其他文件冲突的情况
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
                ]  # 0102: 加入根据Tm计算的方式
            except:
                return_info[f"Tm_cut"] = 0

    return return_info


# 主要设计primer的函数
def primer_design_from_MAS(
    mas_fna="designTool/gut/v3v4/forward_conserved_afterMuscle.fasta",
    degebase_cutoff=0.01,  # 超过多少频率的base会被考虑到degebase里面；太低频率的认为是噪声
    deletion_cutoff=0.99,  #
    mismatch_cutoff=0.99,  # 小于这个阈值认为是mismatch
    primer_lens_list=[20],  # primer的长度允许区间
    forward_reverse="forward",
    step_search=1,  # 1230: 加速搜索，扩大范围
):
    recs = SeqIO.parse(mas_fna, "fasta")
    seqs_np = []
    for rec in recs:
        seq_now = str(rec.seq)
        # 把其他dege的base换成一个普通base
        for dege_base, subs_nor_base in degenerate_base_table.items():
            subs_nor_base = subs_nor_base[
                0
            ]  # 1219：这部分因为使用的是SILVA中的参考序列，这种dege base比例本身不高，参考师兄们的说法就处理成第一个base就行
            seq_now = seq_now.replace(dege_base, subs_nor_base)
        seqs_np.append(list(seq_now))
    seqs_np = np.array(seqs_np)

    # 统计每个位点的atgc分布情况
    total_seq_num = seqs_np.shape[0]  # 总共有多少条序列参与了引物设计
    pos_conserved = [{"A": -1, "T": -1, "G": -1, "C": -1, "-": -1}]
    for i in range(seqs_np.shape[1]):
        seqs_this_pos = list(seqs_np[:, i])
        if forward_reverse == "reverse":
            # 找到其互补链
            seqs_this_pos = [atgc_to_complement[x] for x in seqs_this_pos]

        atgc_cnt = Counter(seqs_this_pos)
        pos_conserved.append(dict(atgc_cnt))
    pos_conserved = pd.DataFrame(pos_conserved)
    pos_conserved.fillna(0.0, inplace=True)
    pos_conserved = pos_conserved.iloc[1:, :]
    pos_conserved.reset_index(inplace=True, drop=True)

    # 给每一个位点一个保守还是潜在degenerate的label
    # 1231：进行修改以更好地符合
    def pos_conserve_label(atgc_freq=[], degebase_cutoff=0.1, deletion_cutoff=0.98):
        atgc_poss = [x / sum(atgc_freq) for x in atgc_freq]
        # 第一个判断条件，是否可能是保守位点或者degenerate位点
        if atgc_poss[-1] > 0.05:  # 考虑改成--0.001 or 0.01（有必要的）
            if atgc_poss[-1] > deletion_cutoff:
                return "deletion_couldDrop", "-", max(atgc_poss), "-", max(atgc_poss)
            else:
                return "deletion_soMuch", "-", max(atgc_poss), "-", max(atgc_poss)
        if max(atgc_poss) > super_conserve_pos_cutoff:
            # 这边是超级保守的-0.995以上
            return (
                "super_conserve_pos",
                atgc_order[atgc_poss.index(max(atgc_poss))],
                max(atgc_poss),
                atgc_order[atgc_poss.index(max(atgc_poss))],
                max(atgc_poss),
            )

        # 创建degenerate base
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
                ],  # 主要是为了degenerate太多的时候选择性地不dege
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

    # 找到可能的primer序列
    # 首先清除超过99%都是-的position
    pos_conserved = pos_conserved[
        pos_conserved["pos_type"] != "deletion_couldDrop"
    ]  # drop 掉很多deletion的序列
    print(f"After drop - positon: {pos_conserved.shape}")
    pos_conserved.reset_index(inplace=True, drop=True)

    potential_primers = []
    for primer_lens in primer_lens_list:
        for start_pos in tqdm(range(
            0, pos_conserved.shape[0] + 1 - primer_lens, step_search
        )):  # 1230: step_search加速搜索
            end_pos = start_pos + primer_lens
            primer_now = pos_conserved.iloc[start_pos:end_pos, :].copy()
            primer_now.reset_index(inplace=True, drop=True)
            ## 接下来进行一些筛选条件，只有全部通过的才能作为可选primer
            # 1. 不能包含过多的-
            if "deletion_soMuch" in list(primer_now["pos_type"]):
                continue
            # 2. 包含少于等于3个degebase: 20240123修改--头4个和最后4个位点都不允许有degenerate base
            for ii in range(len(primer_now) - 3, len(primer_now)):
                # 末尾3个位点不能有dege
                if primer_now.loc[ii, "pos_type"] == "possible_degebase_pos":
                    # if (
                    #     primer_now.loc[ii, "pos_base_freq_max_nodege"]
                    #     < possible_conserve_pos_cutoff
                    # ):
                    #     # 太小了不允许退化
                    #     continue

                    # 退化为普通位点
                    primer_now.loc[ii, "pos_type"] = "possible_conserved_pos"
                    primer_now.loc[ii, "pos_base"] = primer_now.loc[
                        ii, "pos_base_nodege"
                    ]
                    primer_now.loc[ii, "pos_base_freq_max"] = primer_now.loc[
                        ii, "pos_base_freq_max_nodege"
                    ]
                    
            for ii in range(0, 3):
                # 头3个位点不能有dege
                if primer_now.loc[ii, "pos_type"] == "possible_degebase_pos":
                    # if (
                    #     primer_now.loc[ii, "pos_base_freq_max_nodege"]
                    #     < possible_conserve_pos_cutoff
                    # ):
                    #     # 太小了不允许退化
                    #     continue

                    # 退化为普通位点
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
                # 进行一些修补
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
            # 3. 包含最多1个mismatch(除开degebase之外)
            possibility_ls = list(primer_now["pos_base_freq_max"])
            mismatch = [x for x in possibility_ls if x < mismatch_cutoff]
            if len(mismatch) > 3: # 0104: 1 -> 3，保留更多可能性
                continue
            # 4. 按照GC含量、Tm来筛选，按照文章给定了筛选标准
            if forward_reverse == "reverse":
                # reverse需要反向
                primer_atgc = "".join(list(primer_now["pos_base"])[::-1])
            else:
                primer_atgc = "".join(list(primer_now["pos_base"]))

            GC_ratio, Tm = TmGc(primer_atgc)
            if GC_ratio < 0.4 or GC_ratio > 0.65:
                continue
            if Tm < 55 or Tm > 66:  #  65
                continue
            # 5. 直接判断seq与模板之间的cover情况，去掉3‘端附近3bp有mismatch的那些
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
            if pri_cov_num_info[f"cov_num_k_{cov_k_by}"] < 0.85:  # 1219: 修改为一个比较大的值了
                # coverage 太小筛出
                continue

            # 通过筛选之后的primer，获取结果
            primer_degenum = pos_type_summary["possible_degebase_pos"]
            primer_base_coverage_ave = np.mean(possibility_ls)
            # 引物bind在Ecoli上的位点
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


# 从seg信息生成
def creat_primer_fna(rna_seq, pos, id, des, fna_file):
    ## create a fna to store the sub-regions
    records = []
    lens = []
    # merge
    for i in range(min(len(rna_seq), 50000)):
        spe_id = id[i]

        # Create records
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

    ## 生成fna
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
    num_every_spe=5000, # 考虑加入num_every_spe = 0.1这种
    representative_seqs_pick_method="random",
    extend_bp_num=50,
    deletion_cutoff=0.99,  # 下面是一些设计引物时候的参数
    mismatch_cutoff=0.99,  # 1219：后面用的时候正常是0.85
    primer_lens_list=[20],
    taxo_df=None,
    SILVA_set_pred_16sDeepSeg=None,
    rand_seed=None,  # 为了保证各个region设计的用的seqs一致（in-silico 1 genus验证）
    rm_tmp_files=True, # 是否删除中间文件
    step_search=1, # 默认每隔一个位置搜索一次，这样更细得到的候选引物更多
):
    # 20240103：如果输入num_every_spe < 1，按照每个genus取一定比例
    if num_every_spe < 1:
        ratio_every_spe = num_every_spe
    else:
        ratio_every_spe = None
        num_every_spe = int(num_every_spe)
    
    # get core species
    core_microbiota = [x.strip().lower() for x in core_microbiota]
    print(f"Environment has {len(core_microbiota)} genus.")

    # extract corresponding v-regions and flanking conserved regions
    print("Target " + target_vs)
    target_id = target_vs + "_target"
    target_v_root = os.path.join(f"{res_root}/{microbiota_target}", target_vs)
    os.makedirs(target_v_root, exist_ok=True)

    target_vs = target_vs.split("v")[1:]
    target_vs = [int(x) for x in target_vs]

    # 用于设计forward and reverse primer的函数
    def design_forwardOrReverse(target_v=3, fr="reverse", num_every_spe = 100, ratio_every_spe = None):  # 以设计v3v4为例子
        if fr == "forward":
            start_primer = f"v{target_v - 1}"
            end_primer = f"v{target_v}"
        else:
            start_primer = f"v{target_v}"
            end_primer = f"v{target_v + 1}"
        ## 思路：提取出v3v4区域左右侧的保守区域，做MSA，之后去掉>3SD的那些序列
        # scale = 1
        # extract regions
        start_pos, end_pos, ids, rna_seqs = [], [], [], []

        for micro_genus in core_microbiota:
            if len(micro_genus.split()) > 1:
                print(f"Maybe unusual genus name: {micro_genus}.")

            # 找寻当前core genus的所属seqs
            taxo_df_micro_now = taxo_df[taxo_df["genus_name"] == micro_genus]
            if ratio_every_spe is not None:
                # 0103：按照比例取一定的seqs
                num_every_spe = int(ratio_every_spe * taxo_df_micro_now.shape[0])
                
            if representative_seqs_pick_method == "lens":
                # 按照长度排序取前几位
                taxo_df_micro_now.sort_values(by="lens", ascending=False, inplace=True)
            else:
                # 随机重新排列
                if rand_seed is not None:
                    taxo_df_micro_now = taxo_df_micro_now.sample(
                        frac=1, random_state=rand_seed
                    )
                else:
                    taxo_df_micro_now = taxo_df_micro_now.sample(frac=1)

            taxo_df_micro_now.reset_index(inplace=True, drop=True)
            taxo_df_micro_now = taxo_df_micro_now.iloc[:num_every_spe, :]  # 本来有+ 1删掉了

            micro_genus_id = list(
                taxo_df_micro_now["silva_id"] # ["silva_id_wrong"] # 20240102: 修改为silva-id呢
            )  # 使用的是wrong的id U -> T 
            
            SILVA_set_pred_16sDeepSeg_t = SILVA_set_pred_16sDeepSeg[
                SILVA_set_pred_16sDeepSeg["silva_id"].isin(micro_genus_id) # ["silva_id_wrong"] # 20240102: 修改为silva-id呢
            ]
            SILVA_set_pred_16sDeepSeg_t = shuffle(SILVA_set_pred_16sDeepSeg_t)
            SILVA_set_pred_16sDeepSeg_t.reset_index(inplace=True, drop=True)
            if SILVA_set_pred_16sDeepSeg_t.shape[0] < num_every_spe - 5:
                # 1218：如果不够num_every_spe个，则取重复补全
                SILVA_set_pred_16sDeepSeg_t = SILVA_set_pred_16sDeepSeg_t.sample(
                    n=num_every_spe, replace=True
                )

            # 11.24: 没找到相应genus直接跳过
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

            ids += list(SILVA_set_pred_16sDeepSeg_t["silva_id"])[:num_every_spe] # ["silva_id_wrong"] # 20240102: 修改为silva-id呢
            rna_seqs += list(SILVA_set_pred_16sDeepSeg_t["16s_rna"])[:num_every_spe]

        print(num_every_spe, 'AAAAA')
        with open(os.path.join(target_v_root, "silva_id_used.txt"), "w") as f:
            # 0102: 保存一下是用了哪些silva id的序列来操作
            f.write("\n".join(ids))

        # 11.5: 添加往外扩增一定bp的操作
        extend_bp_num_fr = extend_bp_num
        if start_primer == "v0" or end_primer == "v10":
            # 如果是forward的v0，则右侧extend两倍的bp；若是reverse的v9，亦然
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
            # TODO: 2023.0305: add otu cluster to reduce the time cost
            f.write(
                "## Commandline for clustering seqs to OTUs with usearch tool.\n"
            )  # add comments

            file_in = os.path.join(target_v_root, fna_name_muscle + ".fasta")
            uniq_file = os.path.join(target_v_root, "uniques.fa")
            usearch_cmd = (
                f"usearch -fastx_uniques {file_in} -fastaout {uniq_file} -sizeout -relabel Uniq\n"
                + f'usearch -cluster_otus {uniq_file} -otus {file_in} -relabel {fna_name_muscle.split("_")[0]}_target_\n'
            )
            # f.write(usearch_cmd) # 这里是去掉了么？

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
            step_search=step_search # 引物设计的时候每隔多少bp选择一个candidate primer
        )

        return fr_primer_df
    
    fg = design_forwardOrReverse(target_vs[0], "forward", num_every_spe = num_every_spe, ratio_every_spe = ratio_every_spe)
    rg = design_forwardOrReverse(target_vs[1], "reverse", num_every_spe = num_every_spe, ratio_every_spe = ratio_every_spe)
    # 删除中间文件
    if rm_tmp_files:
        for rm_file_nm in ["conserved", "pri_cov_temp", "run_muscle.bash"]:
            os.system(f"rm {target_v_root}/*{rm_file_nm}*")
    if (fg is None) or (rg is None):
        return -1  # 没设计成功
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
