import torch
import numpy as np

def accuracy(output, target):
    with torch.no_grad():
        pred = torch.argmax(output, dim=1)
        assert pred.shape[0] == len(target)
        correct = 0
        correct += torch.sum(pred == target).item()
    return correct / len(target)


def top_k_acc(output, target, k=3):
    with torch.no_grad():
        pred = torch.topk(output, k, dim=1)[1]
        assert pred.shape[0] == len(target)
        correct = 0
        for i in range(k):
            correct += torch.sum(pred[:, i] == target).item()
    return correct / len(target)

def unet_acc(
    output,
    target,
):
    pred = torch.sigmoid(output)
    pred = (pred > 0.5).float()

    num_correct = (pred == target).sum()
    num_pixels = torch.numel(pred)

    return num_correct / num_pixels

def unet_acc_v2(
    output,
    target,
):
    true = target.contiguous().view(-1, ).detach().cpu()
    pred = output.contiguous().view(-1, 10).detach().cpu()
    pred = torch.sigmoid(pred)

    accracy = np.mean((torch.argmax(pred,1)==true).numpy())
    # num_correct = (pred == target).sum()
    # num_pixels = torch.numel(pred)

    return accracy
