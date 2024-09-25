import torch.nn.functional as F
import torch


def nll_loss(output, target):
    return F.nll_loss(output, target)


def cross_entropy(output, target):
    output_ = output.contiguous().view(-1, 10)
    target_ = target.contiguous().view(-1, )
    return F.cross_entropy(output_, target_.long())


def bce_loss(output, target):
    return F.binary_cross_entropy_with_logits(output, target=target)


def focal_loss(output, target):
    alpha = 1.0
    gamma = 3.0

    ce_loss = F.cross_entropy(
        output, target, reduction="none"
    )  # important to add reduction='none' to keep per-batch-item loss
    pt = torch.exp(-ce_loss)  # ce loss 实际是log(softmax)所以直接exp就可以了
    Focal_loss = (alpha * (1 - pt) ** gamma * ce_loss).mean()  # mean over the batch

    return Focal_loss
