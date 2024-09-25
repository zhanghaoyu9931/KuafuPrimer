import collections
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import itertools
from sklearn import metrics

import warnings;warnings.filterwarnings('ignore')

def str_list(str='"[1204, 1224]"'):
    str = str.strip('"[').strip(']"')
    str = str.split(", ")
    str = [int(x) for x in str]
    return str

def parse_single_seg_v2(
    seg_np=np.array([0, 0, 1, 0, 2, 2, 0, 0, 3, 0, 0], dtype=np.float32), lens=1600
):
    ## 在test时用于从seg的结果（概率矩阵）预测出具体的HVR信息的函数
    # 先将float的seq利用模糊化处理获得0 1二值化后的seq
    regions = []
    t_pd = 0.6

    def parse_single_v(seg_v, fix_step = 3):
        mx_lens = 0
        mx_region = str([-1, -1])
        
        tol_step = 0
        start = -1
        for i in range(len(seg_v)):
            ## 目前处于v区域
            if seg_v[i] > t_pd:
                if start == -1:
                    start = i
                continue
            if start == -1:
                continue
            if tol_step >= fix_step:
                # 超过最大容忍步数
                t_region = str([start, i - fix_step])
                t_lens = i - fix_step - start
                tol_step = 0
                
                start = -1
                if t_lens > mx_lens:
                    mx_lens = t_lens
                    mx_region = t_region
                    if mx_lens > 50:
                        # lens is long enough to stop
                        break
            else:
                # 容忍1步
                tol_step += 1
                
        if start != -1:
            if lens >= 1600:
                t_region = str([start, 1599])
                t_lens = 1599 - start

                if t_lens > mx_lens:
                    mx_lens = t_lens
                    mx_region = t_region

        if mx_lens < 5:
            # print('Small one.')
            pass

        return mx_region

    for vi in range(9):
        seg_v = np.zeros_like(seg_np)
        seg_v[seg_np == (vi + 1)] = 1
        region_t = parse_single_v(seg_v)
        # if vi == 8:
        #     # retrieve v9 annotaion
        #     if '-1' in region_t and ('-1' not in regions[-1]):
        #         v8_reg = regions[-1]
        #         v8_reg = str_list(v8_reg)
        #         if lens > v8_reg[-1] + 21:
        #             retri_region = [v8_reg[-1] + 21, lens]
        #             region_t = str(retri_region)
        #             print('Save v9.')
                
        regions.append(region_t)
        
    return regions