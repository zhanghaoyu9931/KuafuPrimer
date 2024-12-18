import os
import random
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

## global vars
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
blast_db_temp_dir = "output/temp"  # temperate dir for blast db
blastn_columns = [
    "query_acc.ver",
    "subject_acc.ver",
    "%_identity",
    "alignment_length",
    "mismatches",
    "gap_opens",
    "q._start",
    "q._end",
    "s._start",
    "s._end",
    "evalue",
    "bit_score",
]


## Useful functions
def blast_cmd(query_seq, db="", out_path="", evalue=1000):  # original evalue=1000
    cmd = f"blastn -db {db} -query "
    cmd += query_seq
    cmd += " -outfmt 7 "
    cmd += "-out "
    cmd += out_path
    cmd += f" -task blastn-short -word_size 4 -evalue {evalue} -max_target_seqs 1000000 -num_threads 16"
    os.system(cmd)


def Tm_GCcal(primer_atgc="AY", sup_seq="AA"):
    # calculate Tm of primer
    if len(primer_atgc) != len(sup_seq):
        return 0

    Tm_C = {"A": 2, "T": 2, "G": 4, "C": 4}
    Tm = 0
    for i in range(len(primer_atgc)):
        p_atgc, s_atgc = primer_atgc[i], sup_seq[i]
        p_atgc_list = degenerate_base_table[p_atgc]
        if s_atgc not in p_atgc_list:
            continue
        Tm += Tm_C[s_atgc] / len(p_atgc_list)  # Weighted degenerate base Tm
    return Tm


def primer_binding_probability(primer_sequence, template_sequence, K=3):
    # func to calculate the successful amplification probability between primer and template
    if len(primer_sequence) != len(template_sequence):
        return 0
    length = len(primer_sequence)
    match_list = []

    ## get the match types in each base
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

    ## calculate the sucessful amplification probability
    # 1. there is no mismatch or degenerate match in the last 3 bases
    match_list_3end = [x for x in match_list[-3:]]
    if len([x for x in match_list_3end if "mismatch" in x]) > 0:
        return 0
    if (
        len([x for x in match_list_3end if "Degenerate match" in x]) > 1
    ):  # at most 1 degenrate match in the last 3 bases
        return 0
    # 2. mismatch number <= K
    mismatch_num = len([x for x in match_list if "mismatch" in x])
    if mismatch_num > K:
        return 0
    # 3. sucessful amplification
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
    if len(lines) == 0:
        blast_df = pd.DataFrame(columns=blastn_columns)
    else:
        blast_df = pd.DataFrame(lines)
        blast_df.columns = index_line
        numeric_cols = [
            "%_identity",
            "alignment_length",
            "mismatches",
            "gap_opens",
            "q._start",
            "q._end",
            "s._start",
            "s._end",
            "evalue",
            "bit_score",
        ]
        for col in numeric_cols:
            blast_df[col] = pd.to_numeric(blast_df[col], errors="coerce") # 20241218: transform cols to numeric

    return blast_df


def get_PP_position_Ecoli_K12(
    f_pri,
    r_pri,
    db="Model_data/Ecoli_K12/Ecoli_K12",
):
    # get the position of primer in Ecoli_K12 genome
    rand_int = random.randint(9931, 99419)
    temp_dir = f"temp_align_{rand_int}"
    os.makedirs(temp_dir, exist_ok=True)
    records = [
        SeqRecord(seq=Seq(f_pri), id="forward", description=""),
        SeqRecord(seq=Seq(r_pri), id="reverse", description=""),
    ]

    with open(f"./{temp_dir}/temp_al.fna", "w") as f:
        SeqIO.write(records, f, "fasta")

    try:
        blast_cmd(
            f"./{temp_dir}/temp_al.fna", db=db, out_path=f"./{temp_dir}/temp_blast.txt"
        )
        blast_df = parse_blastTxt(blast_txt=f"./{temp_dir}/temp_blast.txt")

        info_t = {}
        for pri_ty in ["forward", "reverse"]:
            df_ = blast_df[blast_df["query_acc.ver"] == pri_ty].reset_index(drop=True)
            # forward primer needs to be in the forward strand, and reverse primer needs to be in the reverse strand
            if pri_ty == "forward":
                df_ = df_[df_["s._end"] > df_["s._start"]]
            else:
                df_ = df_[df_["s._end"] < df_["s._start"]]
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


# check the off-target amplification of primer
def offTarget_amplicon_check(
    f_pri,
    r_pri,
    offTarget_fasta="Model_data/OffTarget_amplicon_check/offTarget_reference_seqs.fasta",
    permitted_mismatch=1,
    stringent_mode=False,  # if True, the any one of the primer pair has off-target amplification, then remove this primer pair
    verbose=False,
):
    # make blast db
    offTarget_db = offTarget_fasta.replace(".fasta", "")
    cmd = f"makeblastdb -in {offTarget_fasta} -dbtype nucl -out {offTarget_db}"
    if not os.path.exists(offTarget_db + ".nin"):
        os.system(cmd)

    # create a temp dir
    rand_int = random.randint(9931, 99419)
    temp_dir = f"temp_offTarget_{rand_int}"
    os.makedirs(temp_dir, exist_ok=True)
    records = [
        SeqRecord(seq=Seq(f_pri), id="forward", description=""),
        SeqRecord(seq=Seq(r_pri), id="reverse", description=""),
    ]
    with open(f"./{temp_dir}/temp_al.fna", "w") as f:
        SeqIO.write(records, f, "fasta")
    # get the off-target seqs and description
    offTarget_df = SeqIO.parse(offTarget_fasta, "fasta")
    offTarget_df = pd.DataFrame(
        [[x.id, str(x.seq), x.description] for x in offTarget_df],
        columns=["seq_id", "ATGC", "description"],
    )
    offTarget_df.set_index("seq_id", inplace=True)

    def get_ref_atgc(x):
        # todo corner case: maybe (start_ - 1) < 0 or (end_ - 1) > len(full_ref) ?
        full_ref = offTarget_df.loc[x["subject_acc.ver"], "ATGC"]
        # according to the forward or reverse annotation:
        fr_label = x["query_acc.ver"]
        start_, end_ = int(x["s._start"]), int(x["s._end"])
        start_q, end_q = int(x["q._start"]), int(x["q._end"])
        if end_ > start_:
            start_ = start_ - (start_q - 1)
            end_ = end_ + (len(f_pri) - end_q)
            ref_seq_ = full_ref[start_ - 1 : end_]
            ref_seq_ = "".join([degenerate_base_table[s][0] for s in list(ref_seq_)])
            if verbose:
                print(len(f_pri), f_pri, len(ref_seq_), ref_seq_)
            return ref_seq_
        else:
            start_ = start_ + (start_q - 1)
            end_ = end_ - (len(r_pri) - end_q)
            ref_seq_ = full_ref[start_ - 1 : end_ - 2 : -1]
            ref_seq_ = "".join(
                [
                    (atgc_to_complement[degenerate_base_table[s][0]])
                    for s in list(ref_seq_)
                ]
            )
            return ref_seq_

    offTarget_info_t = {"forward": [], "reverse": []}
    try:
        blast_cmd(
            f"./{temp_dir}/temp_al.fna",
            db=offTarget_db,
            out_path=f"./{temp_dir}/temp_blast.txt",
            evalue=100,
        )
        blast_df = parse_blastTxt(blast_txt=f"./{temp_dir}/temp_blast.txt")
        numeric_cols = [
            "%_identity",
            "alignment_length",
            "mismatches",
            "gap_opens",
            "q._start",
            "q._end",
            "s._start",
            "s._end",
            "evalue",
            "bit_score",
        ]
        for col in numeric_cols:
            blast_df[col] = pd.to_numeric(blast_df[col], errors="coerce")
        blast_df = blast_df[blast_df["alignment_length"] > 9].reset_index(drop=True)
        blast_df = blast_df[
            (
                (blast_df["query_acc.ver"] == "forward")
                & (blast_df["s._end"] > blast_df["s._start"])
            )
            | (
                (blast_df["query_acc.ver"] == "reverse")
                & (blast_df["s._end"] < blast_df["s._start"])
            )
        ].reset_index(drop=True)

        blast_df = blast_df.loc[
            blast_df.groupby(["query_acc.ver", "subject_acc.ver"])["evalue"].idxmin(),
        ].reset_index(
            drop=True
        )  # de-replicated
        blast_df["ref_seq"] = blast_df.apply(lambda x: get_ref_atgc(x), axis=1)
        if verbose:
            print(
                "blast_df: ",
                blast_df[
                    [
                        "query_acc.ver",
                        "ref_seq",
                        "q._start",
                        "q._end",
                        "s._start",
                        "s._end",
                    ]
                ],
            )

        for pri_ty in ["forward", "reverse"]:
            df_ = blast_df[blast_df["query_acc.ver"] == pri_ty].reset_index(drop=True)
            df_["pri_seq"] = f_pri if pri_ty == "forward" else r_pri
            if len(df_) == 0:
                continue
            df_["bind_prob"] = df_.apply(
                lambda x: primer_binding_probability(
                    x["pri_seq"], x["ref_seq"], K=permitted_mismatch
                ),
                axis=1,
            )
            df_binded = df_[df_["bind_prob"] > 0]
            offTarget_info_t[pri_ty] += df_binded["subject_acc.ver"].tolist()
    except:
        if verbose:
            print("Warning in off-target amplification check!")

    if verbose:
        print("off-target seqs of each primer: ", offTarget_info_t)
    if stringent_mode:
        offTarget_info_final = list(
            set(offTarget_info_t["forward"] + offTarget_info_t["reverse"])
        )
    else:
        offTarget_info_final = [
            x for x in offTarget_info_t["forward"] if x in offTarget_info_t["reverse"]
        ]

    os.system(f"rm -r ./{temp_dir}")
    return offTarget_info_final, offTarget_info_t


if __name__ == "__main__":
    # test of off-target amplification check function
    
    off_target_test_df = []
    # offTarget_fasta="Model_data/OffTarget_amplicon_check/MITOBANK_6w.fasta" # 6w human mitochondria seqs
    # offTarget_fasta = "Model_data/OffTarget_amplicon_check/offTarget_olive_seqs.fasta"  # plant chloroplast seqs
    offTarget_fasta="Model_data/OffTarget_amplicon_check/offTarget_reference_seqs.fasta" # human mitochondria seqs
    
    uni_pris = pd.read_excel(
        "Model_data/OffTarget_amplicon_check/primers_offTarget_test.xlsx"
    )
    for pri_i in range(0, uni_pris.shape[0]):
        pri_f, pri_r = uni_pris.loc[pri_i, ["forward_seq", "reverse_seq"]]
        pri_i_nm = uni_pris.loc[pri_i, "pri_nm"]

        verbose = False
        offTarget_res, offTarget_detail = offTarget_amplicon_check(
            pri_f,
            pri_r,
            permitted_mismatch=0, # 5 for mitochondria, 0 for chloroplast
            stringent_mode=False, # True for mitochondria, False for chloroplast
            offTarget_fasta=offTarget_fasta,
            verbose=verbose,
        )
        off_target_test_df.append(
            {
                "pri_nm": uni_pris.loc[pri_i, "pri_nm"],
                "off_target_seqs": offTarget_res,
                "off_target_seqs_forward": offTarget_detail["forward"],
                "off_target_seqs_reverse": offTarget_detail["reverse"],
            }
        )
    off_target_test_df = pd.DataFrame(off_target_test_df)
    offTarget_host_nm = os.path.basename(offTarget_fasta).split(".")[0]
    off_target_test_df.to_csv(
        f"Model_data/OffTarget_amplicon_check/offTarget_test_{offTarget_host_nm}.csv",
        index=False,
    )
