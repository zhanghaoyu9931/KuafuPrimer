import os
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
os.environ["CUDA_LAUNCH_BLOCKING"] = "1"

import argparse
import pandas as pd
import numpy as np
import torch
torch.backends.cudnn.deterministic = True # set to True to avoid cuda error
torch.backends.cudnn.benchmark = True # set to True to avoid cuda error

from tqdm import tqdm
import data_loader.data_loaders as module_data
import model.model as module_arch
from parse_config import ConfigParser
from utils.res_util import *
import torch.nn as nn
import torch.nn.functional as F
from utils import prepare_device
import time


def main(config, output_df_to="", input_csv="", plot_auc=""):
    start = time.time()
    
    fold_dir = os.path.dirname(output_df_to)
    os.makedirs(fold_dir, exist_ok=True)
    logger = config.get_logger("test")

    # setup data_loader instances
    data_loader = getattr(module_data, config["data_loader"]["type"])(
        input_csv,
        batch_size=256,
        shuffle=False,
        validation_split=0.0,
        training=False,
        num_workers=8,
        data_type=config["data_loader"]["args"]["data_type"],
        seq_lens=config["data_loader"]["args"]["seq_lens"],
    )

    # build model architecture
    model = config.init_obj("arch", module_arch)
    # logger.info(model)
    
    # prepare model for testing
    device, device_ids = prepare_device(config["n_gpu"])
    
    model = model.to(device)
    if config["n_gpu"] > 1:
        model = torch.nn.DataParallel(model)
        
    logger.info("Loading checkpoint: {} ...".format(config.resume))
    checkpoint = torch.load(config.resume)
    state_dict = checkpoint["state_dict"]
    model.load_state_dict(state_dict)
    model.eval()

    seg_info = []
    # for auc and confusion matrix figs
    if plot_auc is not None:
        position_wise_pred = []
        position_wise_true = []

    with torch.no_grad():
        # run the 16sDeepSeg module to annotate 16s rRNA genes
        for i, (data, target, lens) in enumerate(tqdm(data_loader)):
            data, target = data.to(device), target.to(device)
            data = data.permute(0, 2, 1).to(torch.float)
            output = model(data)
            output = output.permute(0, 2, 1)
            # the last layer of the model was not a sigmoid layer
            output = nn.Sigmoid()(output)
            
            # plot the auc curve
            if plot_auc is not None:
                position_wise_true.append(target.detach().cpu().contiguous().numpy().flatten())
                position_wise_pred.append(output.detach().cpu().contiguous().view(-1, 10).numpy()) # 10 for 9 V-reigons and conserved region
            
            for i in range(output.size()[0]):
                id_out = torch.argmax(output[i, :, :], dim=1)
                id_out =  np.array(id_out.cpu(), dtype = np.int)
                len_out = float(lens[i].cpu())
                t = parse_single_seg_v2(id_out, len_out)

                seg_info.append(t)
            # save sample images, or do something with output here

    # record the running time
    end = time.time()
    print(f'Running time of this procedure: {end - start:.4f} s.')
    
    with open(os.path.join(fold_dir, 'DeepAnno16_running_time.txt'), 'w') as f:
        f.write(f'Running time of this procedure: {end - start:.4f} s.')
        
    # save the segmentation results
    input_df = pd.read_csv(input_csv)
    seg_cols = [f"v{i+1}" for i in range(9)]
    
    df_seg_output = pd.DataFrame(seg_info)
    df_seg_output.columns = seg_cols
    df_seg_output['id'] = list(input_df['id'])
    df_seg_output = pd.DataFrame(df_seg_output, columns=['id'] + seg_cols)
    df_seg_output.to_csv(output_df_to, index=False)

    # save the position-wise prediction for auc and confusion matrix figs
    if plot_auc is not None:
        position_wise_pred = np.concatenate(position_wise_pred, axis=0)
        position_wise_true = np.concatenate(position_wise_true, axis=0)
        np.save(f'{fold_dir}/position_wise_pred.npy', position_wise_pred)
        np.save(f'{fold_dir}/position_wise_true.npy', position_wise_true)
        
if __name__ == "__main__":
    args = argparse.ArgumentParser(description="PyTorch Template")
    args.add_argument(
        "-p",
        "--plot_auc",
        default=None,
        type=str,
        help="plot position-wise auc and con(default: None)",
    )
    args.add_argument(
        "-c",
        "--config",
        default=None,
        type=str,
        help="config file path (default: None)",
    )
    args.add_argument(
        "-r",
        "--resume",
        default=None,
        type=str,
        help="path to latest checkpoint (default: None)",
    )
    args.add_argument(
        "-d",
        "--device",
        default=None,
        type=str,
        help="indices of GPUs to enable (default: all)",
    )
    args.add_argument(
        "-o",
        "--output",
        default="/data1/hyzhang/Projects/ONT_proj/16sRNA_seg/res/paper_fig/incomplete_pred_512.csv",
        type=str,
        help="output csv file to.",
    )
    args.add_argument(
        "-i",
        "--input",
        default="/data1/hyzhang/Projects/ONT_proj/16sRNA_seg/application2_phylogenetic_tree/res/test_random_E.coli/ori_16s_rna.csv",
        type=str,
        help="input csv file to.",
    )

    config = ConfigParser.from_args(args)
    output = args.parse_args().output
    input = args.parse_args().input
    plot_auc = args.parse_args().plot_auc
    main(config, output, input, plot_auc)
