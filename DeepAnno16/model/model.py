import torch.nn as nn
import torch.nn.functional as F
import torch
import sys

if __name__ == "__main__":
    sys.path.append("/data1/hyzhang/Projects/ONT_proj/16sRNA_seg/src")
    from metric import *
from base import BaseModel


# 1D Unet
class conbr_block(nn.Module):
    def __init__(self, in_layer, out_layer, kernel_size, stride, dilation):
        super(conbr_block, self).__init__()

        self.conv1 = nn.Conv1d(
            in_layer,
            out_layer,
            kernel_size=kernel_size,
            stride=stride,
            dilation=dilation,
            padding=3,
            bias=True,
        )
        self.bn = nn.BatchNorm1d(out_layer)
        self.relu = nn.ReLU()

    def forward(self, x):
        x = self.conv1(x)
        x = self.bn(x)
        out = self.relu(x)

        return out


class se_block(nn.Module):
    def __init__(self, in_layer, out_layer):
        # in_layer and out_layer are the same
        super(se_block, self).__init__()

        self.conv1 = nn.Conv1d(in_layer, out_layer // 8, kernel_size=1, padding=0)
        self.conv2 = nn.Conv1d(out_layer // 8, in_layer, kernel_size=1, padding=0)
        self.fc = nn.Linear(1, out_layer // 8)
        self.fc2 = nn.Linear(out_layer // 8, out_layer)
        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):

        x_se = nn.functional.adaptive_avg_pool1d(x, 1)
        x_se = self.conv1(x_se)
        x_se = self.relu(x_se)
        x_se = self.conv2(x_se)
        x_se = self.sigmoid(x_se)

        # x_out = torch.add(x, x_se)
        x_out = x * x_se
        return x_out


class re_block(nn.Module):
    def __init__(self, in_layer, out_layer, kernel_size, dilation):
        super(re_block, self).__init__()

        self.cbr1 = conbr_block(in_layer, out_layer, kernel_size, 1, dilation)
        self.cbr2 = conbr_block(out_layer, out_layer, kernel_size, 1, dilation)
        self.seblock = se_block(out_layer, out_layer)

    def forward(self, x):
        x_re = self.cbr1(x)
        x_re = self.cbr2(x_re)
        x_re = self.seblock(x_re)
        ## simple implement
        residual = x[:, :, : x_re.size()[2]]
        x_out = torch.add(residual, x_re)
        return x_out


class UNET_1D(nn.Module):
    def __init__(self, input_dim, layer_n, kernel_size, depth, output_class = 1):
        super(UNET_1D, self).__init__()
        self.input_dim = input_dim
        self.layer_n = layer_n # 128
        self.kernel_size = kernel_size # 7
        self.depth = depth # 3

        self.AvgPool1D1 = nn.AvgPool1d(input_dim, stride=5)
        self.AvgPool1D2 = nn.AvgPool1d(input_dim, stride=25)

        self.layer1 = self.down_layer(
            self.input_dim, self.layer_n, self.kernel_size, 1, 2
        )
        self.layer2 = self.down_layer(
            self.layer_n, int(self.layer_n * 2), self.kernel_size, 5, 2
        )
        self.layer3 = self.down_layer(
            int(self.layer_n * 2) + int(self.input_dim),
            int(self.layer_n * 3),
            self.kernel_size,
            5,
            2,
        )
        self.layer4 = self.down_layer(
            int(self.layer_n * 3) + int(self.input_dim),
            int(self.layer_n * 4),
            self.kernel_size,
            5,
            2,
        )
        self.layer5 = self.down_layer(
            int(self.layer_n * 4) + int(self.input_dim),
            int(self.layer_n * 5),
            self.kernel_size,
            4,
            2,
        )
        
        # upsample blocks
        self.cbr_same = conbr_block(
            int(self.layer_n * 4), int(self.layer_n * 4), self.kernel_size, 1, 1
        )

        self.cbr_up1 = conbr_block(
            int(self.layer_n * 7), int(self.layer_n * 3), self.kernel_size, 1, 1
        )
        self.cbr_up2 = conbr_block(
            int(self.layer_n * 5), int(self.layer_n * 2), self.kernel_size, 1, 1
        )
        self.cbr_up3 = conbr_block(
            int(self.layer_n * 3), self.layer_n, self.kernel_size, 1, 1
        )
        self.upsample = nn.Upsample(scale_factor=5, mode="nearest")

        self.outcov = nn.Conv1d(
            self.layer_n,
            output_class, # 721: modify to 10
            kernel_size=self.kernel_size,
            stride=1,
            padding=self.kernel_size // 2,
        )

    def down_layer(self, input_layer, out_layer, kernel, stride, depth):
        block = []
        block.append(conbr_block(input_layer, out_layer, kernel, stride, 1))
        for i in range(depth):
            block.append(re_block(out_layer, out_layer, kernel, 1))
        return nn.Sequential(*block)

    def forward(self, x):
        ori_lens = x.size()[2]

        pool_x1 = self.AvgPool1D1(x)
        pool_x2 = self.AvgPool1D2(x)

        #############Encoder#####################

        out_0 = self.layer1(x)
        out_1 = self.layer2(out_0)

        x = torch.cat([out_1, pool_x1[:, :, : out_1.size()[2]]], 1)
        out_2 = self.layer3(x)

        x = torch.cat([out_2, pool_x2[:, :, : out_2.size()[2]]], 1)
        x = self.layer4(x)
        
        # CBR the same size
        x = self.cbr_same(x)
        
        #############Decoder####################
        def padd1d(up, out):
            padding = (
                0,
                out.size()[2]
                - up.size()[2],  # 前面填充1个单位，后面填充两个单位，输入的最后一个维度则增加1+2个单位，成为8
            )
            up = F.pad(up, padding)
            return up

        up = self.upsample(x)
        up = padd1d(up, out_2)
        up = torch.cat([up, out_2], 1)
        up = self.cbr_up1(up)

        up = self.upsample(up)
        up = padd1d(up, out_1)
        up = torch.cat([up, out_1], 1)
        up = self.cbr_up2(up)
        
        up = self.upsample(up)
        up = padd1d(up, out_0)
        up = torch.cat([up, out_0], 1)
        up = self.cbr_up3(up)

        ## pad to original lens
        padding = (
            0,
            ori_lens - up.size()[2],  # 前面填充1个单位，后面填充两个单位，输入的最后一个维度则增加1+2个单位，成为8
        )
        up = F.pad(up, padding)
        out = self.outcov(up)

        # out = nn.functional.softmax(out,dim=2)

        return out

if __name__ == "__main__":
    # Have a pick at the model
    unet = UNET_1D(5,128,7,3)
    x = torch.rand(1000, 5, 1600)
    print(unet(x).size())
