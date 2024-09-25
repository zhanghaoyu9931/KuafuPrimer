from torch.functional import tensordot
from torch.utils import data
from torchvision import datasets, transforms
from base import BaseDataLoader
from torch.utils.data import Dataset
import numpy as np
import pandas as pd
import torch
import random
from tqdm import tqdm
from utils.util import *

nucl_ambiguous = {
    "Y": ["C", "T"],
    "R": ["A", "G"],
    "W": ["A", "T"],
    "S": ["G", "C"],
    "K": ["T", "G"],
    "M": ["C", "A"],
    "N": ["A", "T", "G", "C"],
    "H": ["A", "T", "C"],
    "V": ["A", "G", "C"],
    "B": ["T", "G", "C"],
    "D": ["A", "G", "T"],
}

class RnaSegDataset10Class(Dataset):
    """
    @description:signal compare dataset class
    """

    nucl_table = ["unknown", "A", "T", "G", "C"]

    def __init__(
        self,
        data_dir,
        batch_size,
        seq_lens=1600,
        shuffle=True,
        validation_split=0.0,
        num_workers=1,
        training=True,
        *args,
        **kwargs,
    ):
        data_df = pd.read_csv(data_dir)
        if "Unnamed" in data_df.columns[0]:
            data_df = data_df.iloc[:, 1:]
        
        if validation_split < 0.001:
            ## only for test
            data_df = data_df.iloc[:, :]
        self.validation_split = validation_split
        seqs = list(data_df["16s_rna"])
        lens = list(data_df["lens"])

        _data = RnaSegDataset10Class._parse_RnaSeq(seqs, seq_lens)
        _label = RnaSegDataset10Class._parse_Label(
            data_df, seq_lens=seq_lens, validation_split=self.validation_split
        )

        self.y = _label
        self.datas = _data  ## size: [*, seq_steps, feature_dim]
        self.seq_lens = lens

    def __len__(self):
        return len(self.datas)

    def __getitem__(self, idx):
        X = self.datas[idx]
        y = self.y[idx]
        len_t = self.seq_lens[idx]

        if self.validation_split < 0.001:
            sample = (X, y, len_t)
        else:
            sample = (X, y)
        return sample

    @staticmethod
    def _parse_RnaSeq(seqs, seq_lens):
        _data = np.zeros(shape=(len(seqs), seq_lens), dtype=np.int32)

        for i, seq in enumerate(tqdm(seqs)):
            for j, nucl in enumerate(seq):
                if j >= seq_lens:
                    break
                if nucl in nucl_ambiguous.keys():
                    # method of parse ambiguous base: need to recheck
                    nucl_choice = nucl_ambiguous[nucl]
                    nucl = random.sample(nucl_choice, 1)[0]
                nucl_index = RnaSegDataset10Class.nucl_table.index(nucl)
                _data[i, j] = nucl_index

        _data = torch.Tensor(_data)
        _data = torch.nn.functional.one_hot(_data.to(torch.int64), 5)
        return _data

    @staticmethod
    def _parse_Label(df, seq_lens, validation_split):
        # modify
        df_t = df.iloc[:, 3:]
        _label = np.zeros(shape=(df.shape[0], seq_lens), dtype=np.int32)

        # 7.19: 为test也添加label
        if df_t.shape[1] < 9:
            # for output of the module, without the label
            _label = torch.Tensor(_label)
            _label = _label.to(torch.float).unsqueeze(2)
            return _label

        for i in tqdm(range(df_t.shape[0])):
            hvrs = df_t.iloc[i, :]
            for vj, hvr in enumerate(list(hvrs)):
                hvr = hvr.strip('"[').strip(']"')
                hvr = hvr.split(", ")
                hvr = [int(x) for x in hvr]

                if -1 in hvr:
                    print('Not clear data!')
                    continue
            
                # modify to a multi-class output
                _label[i][hvr[0] : hvr[1]] = vj + 1

        _label = torch.Tensor(_label)
        _label = _label.to(torch.float).unsqueeze(2)
        return _label

class RnaSegDataLoader(BaseDataLoader):
    """
    16sRNA seg data loading demo using BaseDataLoader
    """

    def __init__(
        self,
        data_dir,
        batch_size,
        seq_lens=16000,
        shuffle=True,
        validation_split=0.0,
        num_workers=1,
        training=True,
        data_type="ATGC",
    ):
        trsfm = transforms.Compose(
            [
                transforms.ToTensor(),
            ]
        )
        self.data_dir = data_dir
        if data_type == "ATGC_10":
            self.dataset = RnaSegDataset10Class(
                data_dir,
                batch_size,
                seq_lens,
                shuffle,
                validation_split,
                num_workers,
                training,
            )
        else:
            self.dataset = None # config wrong, please reset your data_type parameter!

        super().__init__(
            self.dataset, batch_size, shuffle, validation_split, num_workers
        )


if __name__ == "__main__":
    test = RnaSegClassDataset(
        data_dir="/data3/hyzhang/ont/16s_RNA_seg/parsed_data/taxonomic_classify/genus_more_100/samples.csv",
        batch_size=10,
        validation_split=0.1,
    )
    x, y = test[10]
    print(x.size(), y)
