import os
import pandas as pd
import numpy as np
import random
from tqdm import tqdm
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from common_utils import *


# dimer和hairpin检查
def check_dimer_hairpin(
    f_atgc, r_atgc, DimerAndHairpin_exam_root="./DimerAndHairpin_exam"
):
    os.makedirs(DimerAndHairpin_exam_root, exist_ok=True)
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


## 主要函数
def screen_PPs_main(pcr_dir="", rk_by="Genus_accuracy", ge_seq_num_cutoff = 10, ids_select_list=None, ids_neglect_list=None):
    # 先跑这部分为每个环境跑一个总体结果，以及跑一个大的结果出来
    bacdive_designed_res = []
    # for rk_by in ['Pri_amp_eff', 'Genus_accuracy']:

    df_this_evi = {}
    df_this_evi_pri_info = []
    vvs = os.listdir(pcr_dir)
    vvs = [x for x in vvs if os.path.isdir(os.path.join(pcr_dir, x))]

    for vv in vvs:
        for pri_nm in os.listdir(os.path.join(pcr_dir, vv)):
            if "resInfo" in pri_nm:
                df_resInfo_vv = pd.read_csv(os.path.join(pcr_dir, vv, pri_nm))
                df_resInfo_vv = pd.DataFrame(
                    df_resInfo_vv, columns=["pri_nm", "pri_pair"]
                )
                df_resInfo_vv["pri_vregion"] = [vv] * df_resInfo_vv.shape[0]
                # 1106: 获取primer pair的target position信息
                pos_df = {
                    "forward_start": [],
                    "reverse_start": [],
                    "dimer_flag": [],
                    "hairpin_flag": [],
                }
                for pri_i in range(df_resInfo_vv.shape[0]):
                    f_atgc, r_atgc = eval(df_resInfo_vv.loc[pri_i, "pri_pair"])
                    f_atgc, r_atgc = (
                        f_atgc.replace("5'-", "").replace("-3'", "").strip(),
                        r_atgc.replace("5'-", "").replace("-3'", "").strip(),
                    )  # 去除通用引物两侧的5，3标记
                    pos_info = get_PP_position_Ecoli_K12(f_atgc, r_atgc)
                    pos_df["forward_start"].append(pos_info["forward_start"])
                    pos_df["reverse_start"].append(pos_info["reverse_start"])
                    dimer_flag, hairpin_flag = check_dimer_hairpin(f_atgc, r_atgc)
                    pos_df["dimer_flag"].append(dimer_flag)
                    pos_df["hairpin_flag"].append(hairpin_flag)
                for k, v in pos_df.items():
                    df_resInfo_vv[k] = pos_df[k]
                    # df_resInfo_vv['reverse_start'] = pos_df['reverse_start']

                df_this_evi_pri_info.append(df_resInfo_vv)
                continue
            try:
                df_t = pd.read_csv(os.path.join(pcr_dir, vv, pri_nm))
                pri_eff_res = []
            except:
                continue

            # 1123: 规范化Genus名称
            df_t["Genus"] = df_t["Genus"].apply(lambda x: ("_".join(x.split())).lower())
            unGenus_nms = ['unknown', 'uncultured', 'uncultured_bacterium'] 
            df_t['Genus'] = df_t["Genus"].apply(lambda x: 'unknown' if (x in unGenus_nms) else x) # 把未知物种的几个条目整合一下都叫unknown
            
            df_t["Genus_pred"] = (
                df_t["Genus_pred"]
                .fillna("-")
                .apply(lambda x: ("_".join(x.split())).lower())
            )

            taxa_all = list(set(list(df_t["Genus"])))
            tax_num = [""] * len(taxa_all)
            for ge_i, tax in enumerate(taxa_all):
                df_this_tax = df_t[df_t["Genus"] == tax]
                # 20240103: 加入筛选，在或者不在某个list中的ids
                if ids_select_list is not None:
                    df_this_tax = df_this_tax[df_this_tax['silva_id'].isin(ids_select_list)]
                if ids_neglect_list is not None:                  
                    df_this_tax = df_this_tax[~ df_this_tax['silva_id'].isin(ids_neglect_list)]
                    
                df_this_tax.reset_index(inplace=True, drop=True)
                reads_all_n = df_this_tax.shape[0]
                tax_num[ge_i] = reads_all_n
                try:
                    if rk_by == "Pri_amp_eff":
                        # 只计算amp efficient
                        eff = (
                            df_this_tax[df_this_tax["amplicon_size"] > 10].shape[0]
                            / reads_all_n
                        )
                    else:
                        eff = (
                            df_this_tax[
                                df_this_tax["Genus"] == df_this_tax["Genus_pred"]
                            ].shape[0]
                            / reads_all_n
                        )
                        # eff_amplicon = (
                        #     df_this_tax[df_this_tax["amplicon_size"] > 10].shape[0]
                        #     / reads_all_n
                        # )  # 也一并保存了
                except:
                    eff = -1.0
                pri_eff_res.append(eff)
            df_this_evi[pri_nm[:-8]] = pri_eff_res  # 删除_res.csv

    df_this_evi["tax_num"] = tax_num
    df_this_evi["tax_name"] = taxa_all
    df_this_evi = pd.DataFrame(df_this_evi)
    if ge_seq_num_cutoff != 10:
        df_this_evi.to_csv(os.path.join(pcr_dir, f"detail_{rk_by}_taxnumCut{ge_seq_num_cutoff}.csv"), index=False)
    else:
        df_this_evi.to_csv(os.path.join(pcr_dir, f"detail_{rk_by}.csv"), index=False)

    df_this_evi_pri_info = pd.concat(df_this_evi_pri_info, axis=0)
    # 输出最优的引物组合
    df_this_evi = df_this_evi[df_this_evi["tax_num"] > ge_seq_num_cutoff].reset_index(
        drop=True
    )  # 20240103: 筛去那些非常少genus seqs的，改为使用超参数
    df_this_evi_pri_means = df_this_evi.iloc[:, :-2].mean(
        axis=0
    )  # 最后两列是genus和genus number
    best_pri_index = df_this_evi_pri_means.argmax()
    best_pri_nm, best_pri_pcr_res = (
        df_this_evi_pri_means.index[best_pri_index],
        df_this_evi_pri_means[best_pri_index],
    )
    (
        best_pri_pair,
        best_pri_vv,
        best_pri_st,
        best_pri_ed,
        dimer_flag,
        hairpin_flag,
    ) = df_this_evi_pri_info[df_this_evi_pri_info["pri_nm"] == best_pri_nm].values[
        0, 1:
    ]
    print(
        f"Designed ecosystem specific primer pair for {pcr_dir[:-3]} is {best_pri_pair}, targeting {best_pri_vv} V-region, forward and reverse position in E.coli are ({best_pri_st},{ best_pri_ed}), with {best_pri_pcr_res * 100:.2f}% in-silico PCR amplicon accuracy. The dimer_flag and hairpin_flag are ({dimer_flag}, {hairpin_flag})"
    )

    # 计算一下每个环境的平均值,pri的std，用作评判扩增谱系是否均匀的标准
    mean_ge_val = dict(df_this_evi.iloc[:, :-2].mean(axis=0))
    for pri_name, pri_val in mean_ge_val.items():
        pri_std = df_this_evi[pri_name].std()
        bacdive_designed_res.append(
            {
                "pri_nm": pri_name,
                rk_by: pri_val,
                rk_by + "_Std": pri_std,
            }
        )

    # 保存metainfo table，把accuracy也加上并排序
    bacdive_designed_res = pd.DataFrame(bacdive_designed_res)
    df_this_evi_pri_info = pd.merge(
        df_this_evi_pri_info, bacdive_designed_res, on="pri_nm", how="left"
    )
    df_this_evi_pri_info.sort_values(by=rk_by, inplace=True, ascending=False)
    df_this_evi_pri_info.reset_index(inplace=True, drop=True)
    if ge_seq_num_cutoff != 10:
        df_this_evi_pri_info.to_csv(
            os.path.join(pcr_dir, f"pri_metainfo_{rk_by}_taxnumCut{ge_seq_num_cutoff}.csv"), index=False
        )
    else:
        df_this_evi_pri_info.to_csv(
            os.path.join(pcr_dir, f"pri_metainfo_{rk_by}.csv"), index=False
        )

    # 1124：组合并生成一个和uni primer拼接在一起的列表，只是为了后面的比较方便
    top_k = 5
    uni_pris_compare = pd.read_excel(
        "/data1/hyzhang/Projects/EcoPrimer_git/DeepEcoPrimer_v2/data/primers_sequencing_bacdiveComp.xlsx"
    )
    best_pri_pair_info = []
    df_this_evi_pri_info_designed = df_this_evi_pri_info[
        df_this_evi_pri_info["pri_nm"].str.contains("design")
    ].reset_index(drop=True)
    for top_i in range(top_k):
        if top_i >= df_this_evi_pri_info_designed.shape[0]:
            break
        
        (
            best_pri_pair,
            best_pri_vv,
            best_pri_st,
            best_pri_ed,
            dimer_flag,
            hairpin_flag,
        ) = df_this_evi_pri_info_designed.values[top_i, 1:-2]
        best_pri_pair_info.append(
            {
                "seq_type": f"DEcoPrimer_designed",
                "pri_nm": f"DEcoPrimer_designed_rank{top_i+1}",
                "forward_seq": eval(best_pri_pair)[0],
                "reverse_seq": eval(best_pri_pair)[1],
                "lens": f"{int(best_pri_ed) - int(best_pri_st)}bp",
                "platform": "PE",
                "v_region": best_pri_vv,
            }
        )

    best_pri_pair_info = pd.DataFrame(best_pri_pair_info)
    uni_pris_compare = pd.concat(
        [uni_pris_compare, best_pri_pair_info], axis=0
    ).reset_index(drop=True)
    uni_pris_compare.to_csv(
        os.path.join(pcr_dir, "uni_designed_pps_for_compare.csv"), index=False
    )

    return best_pri_pair, best_pri_vv, best_pri_pcr_res


# parse the args
def parse_args():
    parser = argparse.ArgumentParser(
        description="Design best primers for specific environment."
    )
    parser.add_argument(
        "--pcr_dir", type=str, help="Directory to the in-silico PCR res."
    )
    parser.add_argument("--rk_by", type=str, help=".")
    parser.add_argument("--ge_seq_num_cutoff", type=int, default=10, help=".")
    parser.add_argument("--ids_select_neglect", type=str, default="no;no", help=".") # 只包含或者不包含某些silva-ids
    
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    pcr_dir = args.pcr_dir
    rk_by = args.rk_by
    ge_seq_num_cutoff = args.ge_seq_num_cutoff
    ids_select_neglect = args.ids_select_neglect
    ids_select_neglect_list = []
    for ids_s_n in ids_select_neglect.split(';'):
        if ids_s_n.endswith('.txt'):
            with open(ids_s_n, 'r') as f:
                ids_s_n_list = list(f.readlines())
                ids_s_n_list = [x.strip('\n') for x in ids_s_n_list]
        else:
            ids_s_n_list = None # 默认情况不设定该参数
        ids_select_neglect_list.append(ids_s_n_list)
    ids_select_list, ids_neglect_list = ids_select_neglect_list[:2]
    print(ids_select_list, ids_neglect_list)
    
    screen_PPs_main(pcr_dir=pcr_dir, rk_by=rk_by, ge_seq_num_cutoff=ge_seq_num_cutoff, ids_select_list=ids_select_list, ids_neglect_list=ids_neglect_list)
