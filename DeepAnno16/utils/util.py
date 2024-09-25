import json
import torch
import pandas as pd
from pathlib import Path
from itertools import repeat
from collections import OrderedDict
import pickle


def ensure_dir(dirname):
    dirname = Path(dirname)
    if not dirname.is_dir():
        dirname.mkdir(parents=True, exist_ok=False)


def read_json(fname):
    fname = Path(fname)
    with fname.open("rt") as handle:
        return json.load(handle, object_hook=OrderedDict)


def write_json(content, fname):
    fname = Path(fname)
    with fname.open("wt") as handle:
        json.dump(content, handle, indent=4, sort_keys=False)


def inf_loop(data_loader):
    """wrapper function for endless data loader."""
    for loader in repeat(data_loader):
        yield from loader


def prepare_device(n_gpu_use):
    """
    setup GPU device if available. get gpu device indices which are used for DataParallel
    """
    n_gpu = torch.cuda.device_count()
    if n_gpu_use > 0 and n_gpu == 0:
        print(
            "Warning: There's no GPU available on this machine,"
            "training will be performed on CPU."
        )
        n_gpu_use = 0
    if n_gpu_use > n_gpu:
        print(
            f"Warning: The number of GPU's configured to use is {n_gpu_use}, but only {n_gpu} are "
            "available on this machine."
        )
        n_gpu_use = n_gpu
    device = torch.device("cuda:0" if n_gpu_use > 0 else "cpu")
    list_ids = list(range(n_gpu_use))
    return device, list_ids


class MetricTracker:
    def __init__(self, *keys, writer=None):
        self.writer = writer
        self._data = pd.DataFrame(index=keys, columns=["total", "counts", "average"])
        self.reset()

    def reset(self):
        for col in self._data.columns:
            self._data[col].values[:] = 0

    def update(self, key, value, n=1):
        if self.writer is not None:
            self.writer.add_scalar(key, value)
        self._data.total[key] += value * n
        self._data.counts[key] += n
        self._data.average[key] = self._data.total[key] / self._data.counts[key]

    def avg(self, key):
        return self._data.average[key]

    def result(self):
        return dict(self._data.average)


## add by hy: date--11.30
def dump_pkl(obj, pkl_name):
    with open(pkl_name, "wb") as dest:  # 以二进制模式打开文件
        pickle.dump(obj, dest)  # 将文件 以utf-8的格式编码成二进制并写入文件


def load_pkl(pkl_name):
    with open(pkl_name, "rb") as src:  # 定义了 以二进制的方式读取文件, src可以理解为一个数据流管道
        res = pickle.load(src)
    return res


def parse_pkl_ont_data(
    pkl_file="/data3/sfwu/Project/bacteria_pneumoniae/seg_method/signal_new/AB001445.1.1538_signal_dic.pkl",
    sim_id="2",
):
    ont_dict = load_pkl(pkl_file)
    ont_signal = ont_dict[sim_id]
    ont_len = len(ont_signal)

    return ont_signal, ont_len


if __name__ == "__main__":
    parse_pkl_ont_data()
