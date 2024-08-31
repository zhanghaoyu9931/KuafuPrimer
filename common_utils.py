import os
import random
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

## 全局变量
taxa_level = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
atgc_to_complement = {"A": "T", "T": "A", "G": "C", "C": "G", "-": "-"}
degenerate_base_table = {
    "A": ["A"],
    "T": ["T"],
    "G": ["G"],
    "C": ["C"],
    "Y": ["C", "T"],
    "R": ["A", "G"],
    "W": ["A", "T"],
    "S": ["G", "C"],
    "K": ["T", "G"],
    "M": ["C", "A"],
    "N": ["A", "T", "G", "C"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"],
    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
}
blast_db_temp_dir = (
    "/data1/hyzhang/Projects/EcoPrimer_git/DeepEcoPrimer_v2/tools/temp_blastdb"
)


## 有用的funcs
def blast_cmd(query_seq, db="", out_path=""):
    cmd = f"blastn -db {db} -query "
    cmd += query_seq
    cmd += " -outfmt 7 "
    cmd += "-out "
    cmd += out_path
    cmd += " -task blastn-short -word_size 4 -evalue 1000 -max_target_seqs 1000000"
    os.system(cmd)


def Tm_GCcal(primer_atgc="AY", sup_seq="AA"):
    if len(primer_atgc) != len(sup_seq):
        # print('Primer and template seq have different length!')
        return 0

    Tm_C = {"A": 2, "T": 2, "G": 4, "C": 4}
    Tm = 0
    for i in range(len(primer_atgc)):
        p_atgc, s_atgc = primer_atgc[i], sup_seq[i]
        p_atgc_list = degenerate_base_table[p_atgc]
        if s_atgc not in p_atgc_list:
            continue
        Tm += Tm_C[s_atgc] / len(p_atgc_list)  # degenerate base的结合能力要打个折扣
    return Tm


def primer_binding_probability(primer_sequence, template_sequence, K=3):
    # 引物和模板链匹配概率的计算函数
    if len(primer_sequence) != len(template_sequence):
        # print('Primer and template seq have different length!')
        return 0
    length = len(primer_sequence)
    match_list = []

    ## 获取match信息
    for i in range(length):
        primer_base = primer_sequence[i]
        template_base = template_sequence[i]

        if primer_base == template_base:
            match_list.append("match")
        elif (
            primer_base in degenerate_base_table
            and template_base in degenerate_base_table[primer_base]
        ):
            match_list.append("Degenerate match")
        else:
            mismatch = primer_base + "->" + template_base + " mismatch"
            match_list.append(mismatch)

    ## 计算bind概率
    # 1.最后3个位点不能mismatch & 最后三个位点不能degenerate
    match_list_3end = [x for x in match_list[-3:]]
    if len([x for x in match_list_3end if "mismatch" in x]) > 0:
        return 0
    if (
        len([x for x in match_list_3end if "Degenerate match" in x]) > 1
    ):  # 0104: 不允许3‘端存在多余1个degenerate
        return 0
    # 2.不能超过K个错配
    mismatch_num = len([x for x in match_list if "mismatch" in x])
    if mismatch_num > K:
        return 0

    # 3.能够amplify的情况，计算amp概率，最简单都为1
    return 1


def parse_blastTxt(blast_txt="/data3/hyzhang/ont/16s_RNA_seg/res/blast_res/query.txt"):
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

    # os.remove(blast_txt)
    return blast_df


# 获取Ecoli相对位置的函数
def get_PP_position_Ecoli_K12(
    f_pri,
    r_pri,
    db="/data1/hyzhang/Projects/EcoPrimer_git/DeepEcoPrimer_v2/data/ecoli_db_K12/Ecoli_K12",
):
    rand_int = random.randint(9931, 99419)
    temp_dir = f"temp_align_{rand_int}"
    os.makedirs(temp_dir, exist_ok=True)
    records = [
        SeqRecord(seq=Seq(f_pri), id="forward", description=""),
        SeqRecord(seq=Seq(r_pri), id="reverse", description=""),
    ]

    with open(f"./{temp_dir}/temp_al.fna", "w") as f:
        SeqIO.write(records, f, "fasta")

    blast_cmd(
        f"./{temp_dir}/temp_al.fna", db=db, out_path=f"./{temp_dir}/temp_blast.txt"
    )
    blast_df = parse_blastTxt(blast_txt=f"./{temp_dir}/temp_blast.txt")

    info_t = {}
    try:
        for pri_ty in ["forward", "reverse"]:
            df_ = blast_df[blast_df["query_acc.ver"] == pri_ty].reset_index(drop=True)
            ## TODO: 可能还要加入，使得blast的序列长度至少15
            df_.sort_values(by="evalue", inplace=True, ascending=True)
            if pri_ty == "forward":
                info_t[pri_ty + "_start"] = int(df_.loc[0, "s._start"]) - (
                    int(df_.loc[0, "q._start"]) - 1
                )
            else:
                info_t[pri_ty + "_start"] = int(df_.loc[0, "s._start"]) + (
                    int(df_.loc[0, "q._start"]) - 1
                )
    except:
        info_t = {"forward_start": -1, "reverse_start": -1}

    os.system(f"rm -r ./{temp_dir}")
    return info_t
