import random         # 导入 random 模块，用于随机数生成
import torch          # 导入 PyTorch 模块，用于深度学习任务
import torchvision
import numpy as np    # 导入 numpy 模块，用于数值计算
from torch import nn  # 从 PyTorch 中导入神经网络模块
from sklearn import datasets  # 从sklearn引入数据集
from sklearn.model_selection import train_test_split  # 导入 sklearn 库中的 train_test_split 函数，用于数据划分
from sklearn.preprocessing import StandardScaler     # 导入 sklearn 库中的 StandardScaler 类，用于数据标准化
from sklearn.ensemble import VotingClassifier
from sklearn.metrics import roc_curve, auc,precision_recall_curve
import pandas as pd
from torch.optim.lr_scheduler import StepLR,ReduceLROnPlateau
import os
from torch.nn import functional as F
from sklearn.preprocessing import label_binarize
import matplotlib.pyplot as plt
os.environ['CUDA_VISIBLE_DEVICES'] = "2"
#from sklearn.datasets import load_wine
"""
import rpy2.robjects as robjects
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
"""
from torch.utils.data import Dataset,DataLoader

import joblib
import math
from resnet import ResNet18
para_mean=np.zeros(6)
para_std=np.zeros(6)
para_max=np.zeros(6)
para_min=np.zeros(6)

para_mean_b=np.zeros(6)
para_std_b=np.zeros(6)
para_max_b=np.zeros(6)
para_min_b=np.zeros(6)

class FocalLoss(nn.Module):
    def __init__(self, alpha=0.25, gamma=2, reduce=True):
        super(FocalLoss, self).__init__()
        self.alpha = alpha
        self.gamma = gamma
        self.reduce = reduce

    def forward(self, inputs, targets):
        CE_loss = F.cross_entropy(inputs, targets, reduction='none')
        pt = torch.exp(-CE_loss)
        F_loss = self.alpha * (1-pt)**self.gamma * CE_loss

        if self.reduce:
            return torch.mean(F_loss)
        else:
            return F_loss
def normalize_train(x,std=False,time=False,avg_basequality=False,base_num_fraction=False,g_num=False):
    i=0
    if std: i=1
    if time: i=2
    if avg_basequality: i=3
    if base_num_fraction: i=4
    if g_num: i=5

    global para_mean
    para_mean[i] = np.mean(x) 
    global para_std
    para_std[i] = np.std(x) 
    x=(x - para_mean[i]) / (para_std[i] * 1.0)
    global para_max
    para_max[i]=np.max(x)
    global para_min
    para_min[i]=np.min(x)
    print(para_max[i],para_min[i],para_std[i],para_mean[i])
    return (x-para_min[i])/(para_max[i]-para_min[i]) 
def normalize_train_b(x,std=False,time=False,avg_basequality=False,base_num_fraction=False,g_num=False):
    i=0
    if std: i=1
    if time: i=2
    if avg_basequality: i=3
    if base_num_fraction: i=4
    if g_num: i=5

    global para_mean_b
    para_mean_b[i] = np.mean(x) 
    global para_std_b
    para_std_b[i] = np.std(x) 
    x=(x - para_mean_b[i]) / (para_std_b[i] * 1.0)
    global para_max_b
    para_max_b[i]=np.max(x)
    global para_min_b
    para_min_b[i]=np.min(x)
    print(para_max_b[i],para_min_b[i],para_std_b[i],para_mean_b[i])
    return (x-para_min_b[i])/(para_max_b[i]-para_min_b[i]) 
def normalize_test(x,std=False,time=False,avg_basequality=False,base_num_fraction=False,g_num=False):
    i=0
    if std: i=1
    if time: i=2
    if avg_basequality: i=3
    if base_num_fraction: i=4
    if g_num: i=5
    x=(x - para_mean[i]) / (para_std[i] * 1.0)
    #print(para_max[i],para_min[i],para_std[i],para_mean[i])
    return (x-para_min[i])/(para_max[i]-para_min[i]) 
def normalize_test_b(x,std=False,time=False,avg_basequality=False,base_num_fraction=False,g_num=False):
    i=0
    if std: i=1
    if time: i=2
    if avg_basequality: i=3
    if base_num_fraction: i=4
    if g_num: i=5
    x=(x - para_mean_b[i]) / (para_std_b[i] * 1.0)
    #print(para_max[i],para_min[i],para_std[i],para_mean[i])
    return (x-para_min_b[i])/(para_max_b[i]-para_min_b[i]) 
class MotifNet3(nn.Module):
    def __init__(self, input_size=22, hidden_size=128, num_classes=3):
        super(MotifNet3, self).__init__()
        self.layers1 = nn.Sequential(nn.Linear(12,32),nn.ReLU(),nn.Linear(32, 21))
    def forward(self,mean,base_num_fraction,g_num,avg_basequality,std,time,t_test_pval):
        mean_1 = self.layers1(mean)
        return mean_1
class MotifNet2(nn.Module):
    def __init__(self, input_size=22, hidden_size=128, num_classes=21):
        super(MotifNet2, self).__init__()
        mid_channels = 8
        mid_len = 10
        kernel_size=3
        self.layers1 = nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=kernel_size, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),Flatten(),nn.Linear(mid_channels*mid_len, num_classes))
        """
        self.layers2 = nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=kernel_size, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),nn.Conv1d(mid_channels, mid_channels, kernel_size=kernel_size, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),Flatten(),nn.Linear(mid_channels*mid_len, num_classes))
        self.layers3 = nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=kernel_size, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),nn.Conv1d(mid_channels, mid_channels, kernel_size=kernel_size, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),Flatten(),nn.Linear(mid_channels*mid_len, num_classes))
        self.layers4 = nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=kernel_size, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),nn.Conv1d(mid_channels, mid_channels, kernel_size=kernel_size, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),Flatten(),nn.Linear(mid_channels*mid_len, num_classes))
        self.layers5 = nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=kernel_size, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),nn.Conv1d(mid_channels, mid_channels, kernel_size=kernel_size, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),Flatten(),nn.Linear(mid_channels*mid_len, num_classes))
        
        self.weights = torch.Tensor([0.643180383
        ,0.058239536
        ,0.030066228
        ,0.100720973
        ,0.006191121
        #,0.038674039
        ])
        """
    def forward(self,mean,base_num_fraction,g_num,avg_basequality,std,time,t_test_pval):

        mean_1 = self.layers1(mean)
        """
        base_num_fraction_1 = self.layers2(base_num_fraction)
        avg_basequality_1 = self.layers3(avg_basequality)
        std_1 = self.layers4(std)
        time_1 = self.layers5(time)
        """
        return mean_1#+self.weights[1]*base_num_fraction_1+self.weights[2]*avg_basequality_1+self.weights[3]*std_1+self.weights[4]*time_1##++self.weights[5]*g_num result1+result2+result3
class Flatten(nn.Module):
    def forward(self, input):
        return input.view(input.size(0), -1)
class MyConvBlock(nn.Module):
    def __init__(self, in_channels, out_channels, kernel_size, stride=1):
        super(MyConvBlock, self).__init__()
        self.conv1 = nn.Conv1d(in_channels, in_channels, kernel_size, stride=stride, groups=in_channels)
        self.bn1 = nn.BatchNorm1d(in_channels)
        self.relu1 = nn.ReLU()
        #self.dp1 = nn.Dropout(p=0.5)
        self.conv2 = nn.Conv1d(in_channels, out_channels, 1)
        self.bn2 = nn.BatchNorm1d(out_channels)
        self.relu2 = nn.ReLU()
        #self.dp2 = nn.Dropout(p=0.5)

    def forward(self, x):
        x = self.conv1(x)
        x = self.bn1(x)
        x = self.relu1(x)
        x = self.conv2(x)
        x = self.bn2(x)
        x = self.relu2(x)
        return x
class MotifNet(nn.Module):
    def __init__(self, input_size=22, hidden_size=128, num_classes=3):
        super(MotifNet, self).__init__()
        mid_channels = 4
        mid_len = 6#9
        kernel_size = 4
        mid_len2 = 4#8
        kernel_size2 = 5
        mid_len3 = 8#10
        kernel_size3 = 3
        """
        self.layers1 = nn.ModuleList([
            nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=2, stride=1), nn.BatchNorm1d(mid_channels), nn.ReLU(),Flatten(), nn.Linear(mid_channels*11, num_classes)),
            nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=3, stride=1), nn.BatchNorm1d(mid_channels), nn.ReLU(),Flatten(), nn.Linear(mid_channels*10, num_classes)),
            nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=4, stride=1), nn.BatchNorm1d(mid_channels), nn.ReLU(),Flatten(), nn.Linear(mid_channels*9, num_classes))
        ])
        self.layers2 = nn.ModuleList([
            nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=2, stride=1), nn.BatchNorm1d(mid_channels), nn.ReLU(),Flatten(), nn.Linear(mid_channels*11, num_classes)),
            nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=3, stride=1), nn.BatchNorm1d(mid_channels), nn.ReLU(),Flatten(), nn.Linear(mid_channels*10, num_classes)),
            nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=4, stride=1), nn.BatchNorm1d(mid_channels), nn.ReLU(),Flatten(), nn.Linear(mid_channels*9, num_classes))
        ])
        self.layers3 = nn.ModuleList([
            nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=2, stride=1), nn.BatchNorm1d(mid_channels), nn.ReLU(),Flatten(), nn.Linear(mid_channels*11, num_classes)),
            nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=3, stride=1), nn.BatchNorm1d(mid_channels), nn.ReLU(),Flatten(), nn.Linear(mid_channels*10, num_classes)),
            nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=4, stride=1), nn.BatchNorm1d(mid_channels), nn.ReLU(),Flatten(), nn.Linear(mid_channels*9, num_classes))
        ])
        self.layers4 =nn.ModuleList([
            nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=2, stride=1), nn.BatchNorm1d(mid_channels), nn.ReLU(),Flatten(), nn.Linear(mid_channels*11, num_classes)),
            nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=3, stride=1), nn.BatchNorm1d(mid_channels), nn.ReLU(),Flatten(), nn.Linear(mid_channels*10, num_classes)),
            nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=4, stride=1), nn.BatchNorm1d(mid_channels), nn.ReLU(),Flatten(), nn.Linear(mid_channels*9, num_classes))
        ])
        self.layers5 = nn.ModuleList([
            nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=2, stride=1), nn.BatchNorm1d(mid_channels), nn.ReLU(),Flatten(), nn.Linear(mid_channels*11, num_classes)),
            nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=3, stride=1), nn.BatchNorm1d(mid_channels), nn.ReLU(),Flatten(), nn.Linear(mid_channels*10, num_classes)),
            nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=4, stride=1), nn.BatchNorm1d(mid_channels), nn.ReLU(),Flatten(), nn.Linear(mid_channels*9, num_classes))
        ]) 
        """
        c=4
        i=1
        self.layers1 = nn.Sequential(MyConvBlock(i,c,4),nn.BatchNorm1d(c),nn.ReLU(),)
        self.layers11 = nn.Sequential(MyConvBlock(i,c,3),nn.BatchNorm1d(c),nn.ReLU(),nn.Linear(10, 9),nn.BatchNorm1d(c),nn.ReLU())
        self.layers111 = nn.Sequential(MyConvBlock(i,c,5),nn.BatchNorm1d(c),nn.ReLU(),nn.Linear(8, 9),nn.BatchNorm1d(c),nn.ReLU())
        self.layers1111 = nn.Sequential(nn.Conv1d(3*c, c, kernel_size=1),nn.BatchNorm1d(c))
        self.cls = nn.Sequential(nn.ReLU(),Flatten(),nn.Linear(9*c, num_classes))
        #self.shortcut = nn.Sequential(MyConvBlock(i,c,5),nn.BatchNorm1d(c))
        """
        c=4
        self.layers1 = nn.Sequential(MyConvBlock(1,c,3),nn.BatchNorm1d(c),nn.ReLU())
        self.layers11 = nn.Sequential(MyConvBlock(1,c,2),nn.BatchNorm1d(c),nn.ReLU(),nn.Linear(11, 10),nn.BatchNorm1d(c),nn.ReLU())
        self.layers111 = nn.Sequential(MyConvBlock(1,c,4),nn.BatchNorm1d(c),nn.ReLU(),nn.Linear(9, 10),nn.BatchNorm1d(c),nn.ReLU())
        self.layers1111 = nn.Sequential(nn.Conv1d(3*c, c, kernel_size=1),nn.BatchNorm1d(c))
        self.cls = nn.Sequential(nn.ReLU(),Flatten(),nn.Linear(10*c, num_classes))

        self.layers1 = nn.Sequential(MyConvBlock(5,4,3),nn.BatchNorm1d(4),nn.ReLU())
        self.layers11 = nn.Sequential(MyConvBlock(5,4,2),nn.BatchNorm1d(4),nn.ReLU(),nn.Linear(11, 10),nn.BatchNorm1d(4),nn.ReLU())
        self.layers111 = nn.Sequential(MyConvBlock(5,4,4),nn.BatchNorm1d(4),nn.ReLU(),nn.Linear(9, 10),nn.BatchNorm1d(4),nn.ReLU())
        self.layers1111 = nn.Sequential(nn.Conv1d(3*4, 4, kernel_size=1),nn.BatchNorm1d(4),nn.ReLU(),Flatten(),nn.Linear(4*10, num_classes))

        self.layers1 = nn.Sequential(MyConvBlock(5,2,4))
        self.layers11 = nn.Sequential(MyConvBlock(5,2,3),MyConvBlock(2,2,2))
        self.layers111 = nn.Sequential(MyConvBlock(5,2,2),MyConvBlock(2,2,2),MyConvBlock(2,2,2))
        self.layers1111 = nn.Sequential(nn.Conv1d(6, 2, kernel_size=1),nn.BatchNorm1d(2),nn.ReLU(),Flatten(),nn.Linear(9*2, num_classes))
        
        self.layers1 = nn.Sequential(nn.Conv1d(5, mid_channels, kernel_size=kernel_size, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),nn.Conv1d(mid_channels, mid_channels, kernel_size=kernel_size, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU())
        self.layers11 = nn.Sequential(nn.Conv1d(5, mid_channels-1, kernel_size=kernel_size2, stride=1),nn.BatchNorm1d(mid_channels-1),nn.ReLU(),nn.Conv1d(mid_channels-1, mid_channels-1, kernel_size=kernel_size2, stride=1),nn.BatchNorm1d(mid_channels-1),nn.ReLU(),nn.Linear(mid_len2, mid_len),nn.BatchNorm1d(mid_channels-1),nn.ReLU())
        self.layers111 = nn.Sequential(nn.Conv1d(5, mid_channels+1, kernel_size=kernel_size3, stride=1),nn.BatchNorm1d(mid_channels+1),nn.ReLU(),nn.Conv1d(mid_channels+1, mid_channels+1, kernel_size=kernel_size3, stride=1),nn.BatchNorm1d(mid_channels+1),nn.ReLU(),nn.Linear(mid_len3, mid_len),nn.BatchNorm1d(mid_channels+1),nn.ReLU())
        self.layers1111 = nn.Sequential(nn.Conv1d(3*mid_channels, mid_channels, kernel_size=1, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),Flatten(),nn.Linear(mid_len*mid_channels, num_classes))
        """
        """
        self.layers2 = nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=kernel_size, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU())
        self.layers22 = nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=kernel_size2, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),nn.Linear(mid_len2, mid_len),nn.BatchNorm1d(mid_channels),nn.ReLU())
        self.layers222 = nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=kernel_size3, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),nn.Linear(mid_len3, mid_len),nn.BatchNorm1d(mid_channels),nn.ReLU())
        self.layers2222 = nn.Sequential(nn.Conv1d(3*mid_channels, mid_channels, kernel_size=1, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),Flatten(),nn.Linear(mid_len*mid_channels, num_classes))

        self.layers3 = nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=kernel_size, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU())
        self.layers33 = nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=kernel_size2, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),nn.Linear(mid_len2, mid_len),nn.BatchNorm1d(mid_channels),nn.ReLU())
        self.layers333 = nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=kernel_size3, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),nn.Linear(mid_len3, mid_len),nn.BatchNorm1d(mid_channels),nn.ReLU())
        self.layers3333 = nn.Sequential(nn.Conv1d(3*mid_channels, mid_channels, kernel_size=1, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),Flatten(),nn.Linear(mid_len*mid_channels, num_classes))

        self.layers4 = nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=kernel_size, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU())
        self.layers44 = nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=kernel_size2, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),nn.Linear(mid_len2, mid_len),nn.BatchNorm1d(mid_channels),nn.ReLU())
        self.layers444 = nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=kernel_size3, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),nn.Linear(mid_len3, mid_len),nn.BatchNorm1d(mid_channels),nn.ReLU())
        self.layers4444 = nn.Sequential(nn.Conv1d(3*mid_channels, mid_channels, kernel_size=1, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),Flatten(),nn.Linear(mid_len*mid_channels, num_classes))

        self.layers5 = nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=kernel_size, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU())
        self.layers55 = nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=kernel_size2, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),nn.Linear(mid_len2, mid_len),nn.BatchNorm1d(mid_channels),nn.ReLU())
        self.layers555 = nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=kernel_size3, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),nn.Linear(mid_len3, mid_len),nn.BatchNorm1d(mid_channels),nn.ReLU())
        self.layers5555 = nn.Sequential(nn.Conv1d(3*mid_channels, mid_channels, kernel_size=1, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),Flatten(),nn.Linear(mid_len*mid_channels, num_classes))
        """
        #self.layers6 = nn.Sequential(nn.Linear(input_size,hidden_size),nn.ReLU(),nn.Linear(hidden_size, num_classes))
        #self.weights = nn.Parameter(torch.Tensor([100,1,0.1,0.5,0.3]))
        #self.weights = torch.Tensor([0.643180383,0.058239536,0.030066228,0.100720973,0.006191121#,0.038674039])
        #self.weights = nn.Parameter(torch.Tensor([0.2,0.2,0.2,0.2,0.2]))
        #nn.init.kaiming_normal_(self.weights, nonlinearity='relu')

    def forward(self,mean,base_num_fraction,g_num,avg_basequality,std,time,t_test_pval):
        """
        mean_0 = self.layers1[0](mean)
        mean_1 = self.layers1[1](mean)
        mean_2 = self.layers1[2](mean)
        mean_3 = mean_0+mean_1+mean_2

        base_num_fraction_0 = self.layers1[0](base_num_fraction)
        base_num_fraction_1 = self.layers2[1](base_num_fraction)
        base_num_fraction_2 = self.layers2[2](base_num_fraction)
        base_num_fraction_3 = base_num_fraction_0+base_num_fraction_1+base_num_fraction_2

        avg_basequality_0 = self.layers1[0](avg_basequality)
        avg_basequality_1 = self.layers3[1](avg_basequality)
        avg_basequality_2 = self.layers3[2](avg_basequality)
        avg_basequality_3 = avg_basequality_0+avg_basequality_1+avg_basequality_2

        std_0 = self.layers1[0](std)
        std_1 = self.layers4[1](std)
        std_2 = self.layers4[2](std)
        std_3 = std_0+std_1+std_2
        
        time_0 = self.layers1[0](time)
        time_1 = self.layers5[1](time)
        time_2 = self.layers5[2](time)
        time_3 = time_0+time_1+time_2

        mean_1 = self.layers1(mean)
        base_num_fraction_1 = self.layers2(base_num_fraction)
        avg_basequality_1 = self.layers3(avg_basequality)
        std_1 = self.layers4(std)
        time_1 = self.layers5(time)

        mean_2 = self.layers11(mean)
        base_num_fraction_2 = self.layers22(base_num_fraction)
        avg_basequality_2 = self.layers33(avg_basequality)
        std_2 = self.layers44(std)
        time_2 = self.layers55(time)

        mean_3 = self.layers111(mean)
        base_num_fraction_3 = self.layers222(base_num_fraction)
        avg_basequality_3 = self.layers333(avg_basequality)
        std_3 = self.layers444(std)
        time_3 = self.layers555(time)

        result1=self.weights[0]*mean_1+self.weights[1]*base_num_fraction_1+self.weights[2]*avg_basequality_1+self.weights[3]*std_1+self.weights[4]*time_1
        result2=self.weights[0]*mean_2+self.weights[1]*base_num_fraction_2+self.weights[2]*avg_basequality_2+self.weights[3]*std_2+self.weights[4]*time_2
        result3=self.weights[0]*mean_3+self.weights[1]*base_num_fraction_3+self.weights[2]*avg_basequality_3+self.weights[3]*std_3+self.weights[4]*time_3
        """
        mean_1 = self.layers1(mean)
        mean_2 = self.layers11(mean)
        mean_3 = self.layers111(mean)
        mean4 = self.layers1111(torch.cat((mean_1,mean_2,mean_3), dim=1))
        mean4+=mean_1
        mean4=self.cls(mean4)
        """
        base_num_fraction_1 = self.layers2(base_num_fraction)
        base_num_fraction_2 = self.layers22(base_num_fraction)
        base_num_fraction_3 = self.layers222(base_num_fraction)
        base_num_fraction4 = self.layers2222(torch.cat((base_num_fraction_1,base_num_fraction_2,base_num_fraction_3), dim=1))


        avg_basequality_1 = self.layers3(avg_basequality)
        avg_basequality_2 = self.layers33(avg_basequality)
        avg_basequality_3 = self.layers333(avg_basequality)
        avg_basequality4 = self.layers3333(torch.cat((avg_basequality_1,avg_basequality_2,avg_basequality_3), dim=1))


        std_1 = self.layers4(std)
        std_2 = self.layers44(std)
        std_3 = self.layers444(std)
        std4 = self.layers4444(torch.cat((std_1,std_2,std_3), dim=1))

        time_1 = self.layers5(time)
        time_2 = self.layers55(time)
        time_3 = self.layers555(time)
        time4 = self.layers5555(torch.cat((time_1,time_2,time_3), dim=1))
        """
        return mean4 #+self.weights[1]*base_num_fraction4+self.weights[2]*avg_basequality4+self.weights[3]*std4+self.weights[4]*time4#  
class MotifNet3(nn.Module):
    def __init__(self, input_size=22, hidden_size=128, num_classes=3):
        super(MotifNet3, self).__init__()
        mid_channels = 4
        mid_len = 6#9
        kernel_size = 4
        mid_len2 = 4#8
        kernel_size2 = 5
        mid_len3 = 8#10
        kernel_size3 = 3

        c=32
        i=5
        self.layers1 = nn.Sequential(MyConvBlock(i,c,3),nn.BatchNorm1d(c),nn.ReLU(),MyConvBlock(c,c,2),nn.BatchNorm1d(c),nn.ReLU(),nn.Dropout())
        self.layers11 = nn.Sequential(MyConvBlock(i,c,2),nn.BatchNorm1d(c),nn.ReLU(),MyConvBlock(c,c,2),nn.BatchNorm1d(c),nn.ReLU(),MyConvBlock(c,c,2),nn.BatchNorm1d(c),nn.ReLU(),nn.Dropout())
        self.layers111 = nn.Sequential(MyConvBlock(i,c,4),nn.BatchNorm1d(c),nn.ReLU(),nn.Dropout())
        self.layers1111 = nn.Sequential(nn.Conv1d(3*c, c, kernel_size=1),nn.BatchNorm1d(c),nn.Dropout())
        self.cls = nn.Sequential(nn.ReLU(),Flatten(),nn.Linear(9*c, num_classes))
        #self.cls = nn.Sequential(Flatten(),nn.Linear(9*c, num_classes))
        """
        self.layers1 = nn.Sequential(MyConvBlock(i,c,4),nn.BatchNorm1d(c),nn.ReLU(),MyConvBlock(c,c,3),nn.BatchNorm1d(c),nn.ReLU())
        self.layers11 = nn.Sequential(MyConvBlock(i,c,3),nn.BatchNorm1d(c),nn.ReLU(),MyConvBlock(c,c,3),nn.BatchNorm1d(c),nn.ReLU(),MyConvBlock(c,c,2),nn.BatchNorm1d(c),nn.ReLU())
        self.layers111 = nn.Sequential(MyConvBlock(i,c,6),nn.BatchNorm1d(c),nn.ReLU())
        self.layers1111 = nn.Sequential(nn.Conv1d(3*c, c, kernel_size=1),nn.BatchNorm1d(c))
        self.cls = nn.Sequential(nn.ReLU(),Flatten(),nn.Linear(7*c, num_classes))
        """


    def forward(self,mean):

        mean_1 = self.layers1(mean)
        mean_2 = self.layers11(mean)
        mean_3 = self.layers111(mean)
        mean4 = self.layers1111(torch.cat((mean_1,mean_2,mean_3), dim=1))
        mean4+=(mean_1+mean_2+mean_3)
        mean4=self.cls(mean4)

        return mean4
class CNN_basic(nn.Module):
    def __init__(self,num_classes=21):
        super(CNN_basic, self).__init__()
        mid_channels = 32
        mid_len = 10
        kernel_size = 3 
        self.embedding = nn.Embedding(4, 4)
        self.layers1 = nn.Sequential(nn.Conv1d(9, mid_channels, kernel_size=kernel_size, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),Flatten(),nn.Linear(mid_channels*mid_len, num_classes))
    def forward(self,basic,input):
        basic=self.embedding(basic)
        basic = basic.permute(0, 2, 1)
        input=torch.cat((basic,input),dim=1)
        input=self.layers1(input)

        return input

class origin_dataset_train_basic(Dataset):
    def __init__(self, transform,valid):
        super(origin_dataset_train_basic, self).__init__()
        if valid==False:     
            path_data = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_mean_12_code_filter_iso_offset_upsample_train.pkl'
            path_data_base_num_fraction = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_base_num_fraction_12_code_filter_iso_offset_upsample_train.pkl"
            path_data_std_diff = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_std_diff_12_code_filter_iso_offset_upsample_train.pkl"
            path_data_time_diff = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_time_diff_12_code_filter_iso_offset_upsample_train.pkl"
            path_test_t_test_pval = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_t_test_pval_12_code_filter_iso_offset_upsample_train.pkl"
            path_test_u_test_pval = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_u_test_pval_12_code_filter_iso_offset_upsample_train.pkl"
            path_test_g_num = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_g_num_12_code_filter_iso_offset_upsample_train.pkl"
            path_test_g_value = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_g_value_12_code_filter_iso_offset_upsample_train.pkl"
            path_test_avg_basequality = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_avg_basequality_12_code_filter_iso_offset_upsample_train.pkl"
            path_y_basic_group = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_basic_group_12_code_filter_iso_offset_upsample_train.pkl'
            path_y_offset = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_offset_12_code_filter_iso_offset_upsample_train.pkl'
        else:
            path_data = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_mean_12_code_filter_iso_offset_upsample_valid.pkl'
            path_data_base_num_fraction = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_base_num_fraction_12_code_filter_iso_offset_upsample_valid.pkl"
            path_data_std_diff = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_std_diff_12_code_filter_iso_offset_upsample_valid.pkl"
            path_data_time_diff = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_time_diff_12_code_filter_iso_offset_upsample_valid.pkl"
            path_test_t_test_pval = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_t_test_pval_12_code_filter_iso_offset_upsample_valid.pkl"
            path_test_u_test_pval = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_u_test_pval_12_code_filter_iso_offset_upsample_valid.pkl"
            path_test_g_num = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_g_num_12_code_filter_iso_offset_upsample_valid.pkl"
            path_test_g_value = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_g_value_12_code_filter_iso_offset_upsample_valid.pkl"
            path_test_avg_basequality = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_avg_basequality_12_code_filter_iso_offset_upsample_valid.pkl"
            path_y_basic_group = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_basic_group_12_code_filter_iso_offset_upsample_valid.pkl'
            path_y_offset = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_offset_12_code_filter_iso_offset_upsample_valid.pkl'

        self.data=joblib.load(path_data) 
        self.data_base_num_fraction=joblib.load(path_data_base_num_fraction) 
        self.data_std_diff=joblib.load(path_data_std_diff)
        self.data_time_diff=joblib.load(path_data_time_diff)
        self.data_t_test_pval=joblib.load(path_test_t_test_pval)
        self.data_g_num=joblib.load(path_test_g_num)
        self.data_g_value=joblib.load(path_test_g_value)
        self.data_avg_basequality=joblib.load(path_test_avg_basequality)
        self.label_basic_group=joblib.load(path_y_basic_group)  
        self.label_offset=joblib.load(path_y_offset)

        if valid==False:
            self.data = normalize_train_b(self.data)
            self.data_base_num_fraction = normalize_train_b(self.data_base_num_fraction,base_num_fraction=True)
            self.data_std_diff = normalize_train_b(self.data_std_diff,std=True)
            self.data_time_diff = normalize_train_b(self.data_time_diff,time=True)
            self.data_avg_basequality = normalize_train_b(self.data_avg_basequality,avg_basequality=True)
            self.data_g_num = normalize_train_b(self.data_g_num,g_num=True)
        else:
            self.data = normalize_test_b(self.data)
            self.data_base_num_fraction = normalize_test_b(self.data_base_num_fraction,base_num_fraction=True)
            self.data_std_diff = normalize_test_b(self.data_std_diff,std=True)
            self.data_time_diff = normalize_test_b(self.data_time_diff,time=True)
            self.data_avg_basequality = normalize_test_b(self.data_avg_basequality,avg_basequality=True)
            self.data_g_num = normalize_test_b(self.data_g_num,g_num=True)
        #print(self.data_basic_group.shape)
        print(self.data.shape)

        self.len_data=self.data.shape[0]
        self.transform = transform

        self.data=torch.FloatTensor(self.data)
        self.data_base_num_fraction=torch.FloatTensor(self.data_base_num_fraction)
        self.data_std_diff=torch.FloatTensor(self.data_std_diff)
        self.data_time_diff=torch.FloatTensor(self.data_time_diff)
        self.data_t_test_pval=torch.FloatTensor(self.data_t_test_pval)
        self.data_avg_basequality=torch.FloatTensor(self.data_avg_basequality)
        self.data_g_num=torch.FloatTensor(self.data_g_num)
        self.data_g_value=torch.FloatTensor(self.data_g_value)
        self.label_basic_group=torch.LongTensor(self.label_basic_group)
        self.label_offset=torch.LongTensor(self.label_offset)

    def __len__(self):
        return self.len_data
 
    def __getitem__(self, index):
        #data = torch.FloatTensor(data)
        return  self.data[index],self.data_base_num_fraction[index],self.data_std_diff[index],self.data_time_diff[index],self.data_t_test_pval[index],self.data_avg_basequality[index],self.data_g_num[index],self.data_g_value[index], self.label_basic_group[index],self.label_offset[index]
# Create a mapping from letters to numbers
mapping = {'A': 0, 'C': 1, 'T': 2, 'G': 3}

# Create a function to convert a string to a list of numbers
def convert_string(s):
    return [mapping[char] for char in s]

# Vectorize the function
vfunc = np.vectorize(convert_string, otypes=[list])
class origin_dataset_train(Dataset):
    def __init__(self, transform,valid):
        super(origin_dataset_train, self).__init__()
        if valid==False:     
            path_data = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_mean_12_code_filter_iso_offset_train2.pkl'
            path_data_basic_group = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_basic_group_12_code_filter_iso_offset_train.pkl"
            path_data_base_num_fraction = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_base_num_fraction_12_code_filter_iso_offset_train2.pkl"
            path_data_std_diff = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_std_diff_12_code_filter_iso_offset_train2.pkl"
            path_data_time_diff = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_time_diff_12_code_filter_iso_offset_train2.pkl"
            path_test_t_test_pval = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_t_test_pval_12_code_filter_iso_offset_train.pkl"
            path_test_u_test_pval = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_u_test_pval_12_code_filter_iso_offset_train.pkl"
            path_test_g_num = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_g_num_12_code_filter_iso_offset_train.pkl"
            path_test_g_value = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_g_value_12_code_filter_iso_offset_train.pkl"
            path_test_avg_basequality = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_avg_basequality_12_code_filter_iso_offset_train2.pkl"
            path_y_basic_group = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_basic_group_12_code_filter_iso_offset_train2.pkl'
            path_y_offset = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_offset_12_code_filter_iso_offset_train2.pkl'
        else:
            path_data = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_mean_12_code_filter_iso_offset_valid2.pkl'
            path_data_basic_group = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_basic_group_12_code_filter_iso_offset_valid.pkl"
            path_data_base_num_fraction = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_base_num_fraction_12_code_filter_iso_offset_valid2.pkl"
            path_data_std_diff = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_std_diff_12_code_filter_iso_offset_valid2.pkl"
            path_data_time_diff = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_time_diff_12_code_filter_iso_offset_valid2.pkl"
            path_test_t_test_pval = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_t_test_pval_12_code_filter_iso_offset_valid.pkl"
            path_test_u_test_pval = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_u_test_pval_12_code_filter_iso_offset_valid.pkl"
            path_test_g_num = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_g_num_12_code_filter_iso_offset_valid.pkl"
            path_test_g_value = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_g_value_12_code_filter_iso_offset_valid.pkl"
            path_test_avg_basequality = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_avg_basequality_12_code_filter_iso_offset_valid2.pkl"
            path_y_basic_group = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_basic_group_12_code_filter_iso_offset_valid2.pkl'
            path_y_offset = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_offset_12_code_filter_iso_offset_valid2.pkl'
            """
            path_data = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_mean_12_code_filter_iso_offset.pkl'
            path_data_base_num_fraction = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_base_num_fraction_12_code_filter_iso_offset.pkl"
            path_data_std_diff = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_std_diff_12_code_filter_iso_offset.pkl"
            path_data_time_diff = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_time_diff_12_code_filter_iso_offset.pkl"
            path_test_t_test_pval = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_t_test_pval_12_code_filter_iso_offset.pkl"
            path_test_u_test_pval = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_u_test_pval_12_code_filter_iso_offset.pkl"
            path_test_g_num = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_g_num_12_code_filter_iso_offset.pkl"
            path_test_g_value = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_g_value_12_code_filter_iso_offset.pkl"
            path_test_avg_basequality = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_avg_basequality_12_code_filter_iso_offset.pkl"
            path_y_basic_group = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_basic_group_12_code_filter_iso_offset.pkl'
            path_y_offset = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_offset_12_code_filter_iso_offset.pkl'
            """
        self.data=joblib.load(path_data) 
        self.data_basic_group=joblib.load(path_data_basic_group) 
        self.data_base_num_fraction=joblib.load(path_data_base_num_fraction) 
        self.data_std_diff=joblib.load(path_data_std_diff)
        self.data_time_diff=joblib.load(path_data_time_diff)
        self.data_t_test_pval=joblib.load(path_test_t_test_pval)
        self.data_g_num=joblib.load(path_test_g_num)
        self.data_g_value=joblib.load(path_test_g_value)
        self.data_avg_basequality=joblib.load(path_test_avg_basequality)
        self.label_basic_group=joblib.load(path_y_basic_group)  
        self.label_offset=joblib.load(path_y_offset)
        """
        self.data1=joblib.load(path_data1) 
        self.data_base_num_fraction1=joblib.load(path_data_base_num_fraction1) 
        self.label_basic_group1=joblib.load(path_y_basic_group1)  
        self.label_offset1=joblib.load(path_y_offset1)
        self.data=np.concatenate((self.data, self.data1), axis=0)
        self.label_basic_group=np.concatenate((self.label_basic_group, self.label_basic_group1), axis=0)
        self.label_offset=np.concatenate((self.label_offset, self.label_offset1), axis=0)
        """
        """
        self.data_basic_group=joblib.load(path_data_basic_group) 
        self.data_std_diff=joblib.load(path_data_std) 
        self.data_time_diff=joblib.load(path_data_time)   
        self.data_t_test_pval=joblib.load(path_data_t_test_pval) 
        self.data_u_test_pval=joblib.load(path_data_u_test_pval) 
        self.data_avg_basequality=joblib.load(path_data_avg_basequality) 
        self.data_g_num=joblib.load(path_data_g_num) 
        self.data_g_value=joblib.load(path_data_g_value) 
     
        indices = np.where(self.label_basic_group == 1)[0]
        self.data = self.data[indices]
        self.data_base_num_fraction = self.data_base_num_fraction[indices]
        self.label_offset = self.label_offset[indices]
        """
        self.data_basic_group = vfunc(self.data_basic_group)
        # Convert the list of lists to a 2D array
        self.data_basic_group = np.stack(self.data_basic_group)
        if valid==False:
            self.data = normalize_train(self.data)
            self.data_base_num_fraction = normalize_train(self.data_base_num_fraction,base_num_fraction=True)
            self.data_std_diff = normalize_train(self.data_std_diff,std=True)
            self.data_time_diff = normalize_train(self.data_time_diff,time=True)
            self.data_avg_basequality = normalize_train(self.data_avg_basequality,avg_basequality=True)
            self.data_g_num = normalize_train(self.data_g_num,g_num=True)
        else:
            self.data = normalize_test(self.data)
            self.data_base_num_fraction = normalize_test(self.data_base_num_fraction,base_num_fraction=True)
            self.data_std_diff = normalize_test(self.data_std_diff,std=True)
            self.data_time_diff = normalize_test(self.data_time_diff,time=True)
            self.data_avg_basequality = normalize_test(self.data_avg_basequality,avg_basequality=True)
            self.data_g_num = normalize_test(self.data_g_num,g_num=True)
        #print(self.data_basic_group.shape)
        print(self.data.shape)

        self.len_data=self.data.shape[0]
        self.transform = transform

        self.data_basic_group=torch.LongTensor(self.data_basic_group)
        self.data=torch.FloatTensor(self.data)
        self.data_base_num_fraction=torch.FloatTensor(self.data_base_num_fraction)
        self.data_std_diff=torch.FloatTensor(self.data_std_diff)
        self.data_time_diff=torch.FloatTensor(self.data_time_diff)
        self.data_t_test_pval=torch.FloatTensor(self.data_t_test_pval)
        self.data_avg_basequality=torch.FloatTensor(self.data_avg_basequality)
        self.data_g_num=torch.FloatTensor(self.data_g_num)
        self.data_g_value=torch.FloatTensor(self.data_g_value)
        self.label_basic_group=torch.LongTensor(self.label_basic_group)
        self.label_offset=torch.LongTensor(self.label_offset)

        #self.data_t_test_pval = rollmean(-torch.log10(self.data_t_test_pval))  
        #self.data_g_value = rollmean(-torch.log10(self.data_g_value))  
        """
        self.data_basic_group=torch.LongTensor(self.data_basic_group)
        self.data_std_diff=torch.FloatTensor(self.data_std_diff)
        self.data_time_diff=torch.FloatTensor(self.data_time_diff)
        self.data_t_test_pval=torch.FloatTensor(self.data_t_test_pval)
        self.data_u_test_pval=torch.FloatTensor(self.data_u_test_pval)
        self.data_avg_basequality=torch.FloatTensor(self.data_avg_basequality)
        self.data_g_num=torch.FloatTensor(self.data_g_num)
        self.data_g_value=torch.FloatTensor(self.data_g_value)
        """
        self.data = self.data.unsqueeze(1)
        self.data_base_num_fraction = self.data_base_num_fraction.unsqueeze(1)
        self.data_std_diff = self.data_std_diff.unsqueeze(1)
        self.data_time_diff = self.data_time_diff.unsqueeze(1)
        self.data_avg_basequality = self.data_avg_basequality.unsqueeze(1)
        #self.data_nat =  self.data_nat.unsqueeze(1)
        #self.data_wga =  self.data_wga.unsqueeze(1)
    def __len__(self):
        return self.len_data
 
    def __getitem__(self, index):
        #data = torch.FloatTensor(data)
        return  self.data[index],self.data_base_num_fraction[index],self.data_std_diff[index],self.data_time_diff[index],self.data_t_test_pval[index],self.data_avg_basequality[index],self.data_g_num[index],self.data_g_value[index], self.data_basic_group[index], self.label_basic_group[index],self.label_offset[index]
        #return self.data_basic_group[index],self.label[index]

class origin_dataset_test(Dataset):
    def __init__(self, transform,motif,species):
        super(origin_dataset_test, self).__init__()
        #species="TP"
        path_test = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/{}_pd_from_r_df_test_12_code_filter_iso.pkl".format(species)
        path_test_base_num_fraction = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/{}_pd_from_r_df_test_base_num_fraction_12_code_filter_iso.pkl".format(species)
        path_test_std_diff = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/{}_pd_from_r_df_test_std_diff_12_code_filter_iso.pkl".format(species)
        path_test_time_diff = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/{}_pd_from_r_df_test_time_diff_12_code_filter_iso.pkl".format(species)
        path_test_basic_group = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/{}_pd_from_r_df_test_basic_group_12_code_filter_iso.pkl".format(species)
        path_test_t_test_pval = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/{}_pd_from_r_df_test_t_test_pval_12_code_filter_iso.pkl".format(species)
        path_test_u_test_pval = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/{}_pd_from_r_df_test_u_test_pval_12_code_filter_iso.pkl".format(species)
        path_test_avg_basequality = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/{}_pd_from_r_df_test_avg_basequality_12_code_filter_iso.pkl".format(species)
        path_test_g_num = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/{}_pd_from_r_df_test_g_num_12_code_filter_iso.pkl".format(species)
        path_test_g_value = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/{}_pd_from_r_df_test_g_value_12_code_filter_iso.pkl".format(species)

        pd_from_r_df_test=joblib.load(path_test)   
        pd_from_r_df_test_base_num_fraction=joblib.load(path_test_base_num_fraction)
        pd_from_r_df_test_std_diff=joblib.load(path_test_std_diff)
        pd_from_r_df_test_time_diff=joblib.load(path_test_time_diff)
        pd_from_r_df_test_basic_group=joblib.load(path_test_basic_group)
        pd_from_r_df_test_t_test_pval=joblib.load(path_test_t_test_pval)
        pd_from_r_df_test_avg_basequality=joblib.load(path_test_avg_basequality)
        pd_from_r_df_test_g_num=joblib.load(path_test_g_num)
        pd_from_r_df_test_g_value=joblib.load(path_test_g_value)

        subet_test_data=pd_from_r_df_test[pd_from_r_df_test["label"]==motif]
        subet_test_data_base_num_fraction=pd_from_r_df_test_base_num_fraction[pd_from_r_df_test_base_num_fraction["label"]==motif]
        subet_test_data_std_diff=pd_from_r_df_test_std_diff[pd_from_r_df_test_std_diff["label"]==motif]
        subet_test_data_time_diff=pd_from_r_df_test_time_diff[pd_from_r_df_test_time_diff["label"]==motif]
        subet_test_data_basic_group=pd_from_r_df_test_basic_group[pd_from_r_df_test_basic_group["label"]==motif]
        subet_test_data_t_test_pval=pd_from_r_df_test_t_test_pval[pd_from_r_df_test_t_test_pval["label"]==motif]
        subet_test_data_avg_basequality=pd_from_r_df_test_avg_basequality[pd_from_r_df_test_avg_basequality["label"]==motif]
        subet_test_data_g_num=pd_from_r_df_test_g_num[pd_from_r_df_test_g_num["label"]==motif]
        subet_test_data_g_value=pd_from_r_df_test_g_value[pd_from_r_df_test_g_value["label"]==motif]

        self.data=subet_test_data.values[:,6:18] #6:18
        self.data=self.data.astype(float) 
        self.data = normalize_test(self.data)
        self.data_b=subet_test_data.values[:,1:23]
        self.data_b=self.data_b.astype(float)
        self.data_b = normalize_test_b(self.data_b)

        self.data_base_num_fraction=subet_test_data_base_num_fraction.values[:,6:18] #6:18
        self.data_base_num_fraction=self.data_base_num_fraction.astype(float) 
        
        self.data_base_num_fraction = normalize_test(self.data_base_num_fraction,base_num_fraction=True)
        self.data_base_num_fraction_b=subet_test_data_base_num_fraction.values[:,1:23] #6:18
        self.data_base_num_fraction_b=self.data_base_num_fraction_b.astype(float) 
        #if motif=="CTCGAG_-14":
        #    joblib.dump(self.data_base_num_fraction_b,"/home/nipeng/chenzheng/mulan-methyl/data_pkl/CTCGAG_base_num_fraction.pkl")
        self.data_base_num_fraction_b = normalize_test_b(self.data_base_num_fraction_b,base_num_fraction=True)   

        self.data_std_diff=subet_test_data_std_diff.values[:,6:18] #6:18
        self.data_std_diff=self.data_std_diff.astype(float) 
        self.data_std_diff = normalize_test(self.data_std_diff,std=True)
        self.data_std_diff_b=subet_test_data_std_diff.values[:,1:23] #6:18
        self.data_std_diff_b=self.data_std_diff_b.astype(float)
        self.data_std_diff_b = normalize_test_b(self.data_std_diff_b,std=True)

        self.data_time_diff=subet_test_data_time_diff.values[:,6:18] #6:18
        self.data_time_diff=self.data_time_diff.astype(float)
        self.data_time_diff = normalize_test(self.data_time_diff,time=True)
        self.data_time_diff_b=subet_test_data_time_diff.values[:,1:23] #6:18
        self.data_time_diff_b=self.data_time_diff_b.astype(float)
        self.data_time_diff_b = normalize_test_b(self.data_time_diff_b,time=True)

        self.data_avg_basequality=subet_test_data_avg_basequality.values[:,6:18] #6:18
        self.data_avg_basequality=self.data_avg_basequality.astype(float)
        self.data_avg_basequality = normalize_test(self.data_avg_basequality,avg_basequality=True)
        self.data_avg_basequality_b=subet_test_data_avg_basequality.values[:,1:23]
        self.data_avg_basequality_b=self.data_avg_basequality_b.astype(float)
        self.data_avg_basequality_b = normalize_test_b(self.data_avg_basequality_b,avg_basequality=True)

        self.data_g_num=subet_test_data_g_num.values[:,6:18] #6:18
        self.data_g_num=self.data_g_num.astype(float)
        self.data_g_num = normalize_test(self.data_g_num,g_num=True)
        self.data_g_num_b=subet_test_data_g_num.values[:,1:23]
        self.data_g_num_b=self.data_g_num_b.astype(float)
        self.data_g_num_b = normalize_test_b(self.data_g_num_b,g_num=True)

        #self.data = umap_model.transform(self.data)
        self.data_basic_group=subet_test_data_basic_group.values[:,6:18] #1:23
        self.data_basic_group[self.data_basic_group=='A']=0
        self.data_basic_group[self.data_basic_group=='C']=1
        self.data_basic_group[self.data_basic_group=='T']=2
        self.data_basic_group[self.data_basic_group=='G']=3
        self.data_basic_group = self.data_basic_group.astype(int) 
        self.data_t_test_pval=subet_test_data_t_test_pval.values[:,6:18] #1:23
        self.data_t_test_pval=self.data_t_test_pval.astype(float)
        self.data_g_value=subet_test_data_g_value.values[:,6:18] #1:23
        self.data_g_value=self.data_g_value.astype(float)
        #self.data_basic_group=np.eye(4)[self.data_basic_group] 
        self.len_data=self.data.shape[0]
        self.transform = transform
        self.data=torch.FloatTensor(self.data)
        self.data_std_diff=torch.FloatTensor(self.data_std_diff)
        self.data_time_diff=torch.FloatTensor(self.data_time_diff)
        self.data_basic_group=torch.LongTensor(self.data_basic_group)
        #self.data_t_test_pval=torch.FloatTensor(self.data_t_test_pval)
        self.data_avg_basequality=torch.FloatTensor(self.data_avg_basequality)
        self.data_base_num_fraction=torch.FloatTensor(self.data_base_num_fraction)
        #self.data_g_value=torch.FloatTensor(self.data_g_value)
        self.data_g_num=torch.FloatTensor(self.data_g_num)

        self.data_b=torch.FloatTensor(self.data_b)
        self.data_std_diff_b=torch.FloatTensor(self.data_std_diff_b)
        self.data_time_diff_b=torch.FloatTensor(self.data_time_diff_b)
        self.data_avg_basequality_b=torch.FloatTensor(self.data_avg_basequality_b)
        self.data_base_num_fraction_b=torch.FloatTensor(self.data_base_num_fraction_b)
        self.data_g_num_b=torch.FloatTensor(self.data_g_num_b)

        self.data = self.data.unsqueeze(1)
        self.data_base_num_fraction = self.data_base_num_fraction.unsqueeze(1)
        self.data_std_diff = self.data_std_diff.unsqueeze(1)
        self.data_time_diff = self.data_time_diff.unsqueeze(1)
        self.data_avg_basequality = self.data_avg_basequality.unsqueeze(1)
        #self.data = self.data.unsqueeze(1)
        #self.data=torch.FloatTensor(self.data)
        #self.data_t_test_pval = rollmean(-torch.log10(self.data_t_test_pval))  
        #self.data_g_value = rollmean(-torch.log10(self.data_g_value)) 
    def __len__(self):
        return self.len_data
 
    def __getitem__(self, index):
        #data = torch.FloatTensor(data)
        #return  torch.cat((self.data[index],self.data_basic_group[index]),dim=1)
         return  (self.data[index],self.data_base_num_fraction[index],self.data_std_diff[index],self.data_time_diff[index],self.data_t_test_pval[index],self.data_avg_basequality[index],self.data_g_num[index],self.data_g_value[index], self.data_basic_group[index]),(self.data_b[index],self.data_base_num_fraction_b[index],self.data_avg_basequality_b[index],self.data_std_diff_b[index],self.data_time_diff_b[index])
def sort_and_return_index(lst):
    return np.array([i[0] for i in sorted(enumerate(lst), key=lambda x:x[1], reverse=True)])

# usage
lst = [2, 0, 1, 5, 3, 4]
print(sort_and_return_index(lst))  # Output: [1, 2, 0, 4, 5, 3]
if __name__ == '__main__':

    # 设置随机种子，让模型每次输出的结果都一样
    seed_value = 42
    random.seed(seed_value)                         # 设置 random 模块的随机种子
    np.random.seed(seed_value)                      # 设置 numpy 模块的随机种子
    torch.manual_seed(seed_value)                   # 设置 PyTorch 中 CPU 的随机种子
    #tf.random.set_seed(seed_value)                 # 设置 Tensorflow 中随机种子
    if torch.cuda.is_available():                   # 如果可以使用 CUDA，设置随机种子
        torch.cuda.manual_seed(seed_value)          # 设置 PyTorch 中 GPU 的随机种子
        torch.backends.cudnn.deterministic = True   # 使用确定性算法，使每次运行结果一样
        torch.backends.cudnn.benchmark = False      # 不使用自动寻找最优算法加速运算
    # 检测GPU是否可用
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    batch_size=256 #512
    tol=1e-3
    n_iter_no_change=5
    loss_not_improved_count = 0
    best_loss = float('inf')
    #scaler = StandardScaler() 
    dataset_train=origin_dataset_train(transform=None,valid=False)
    data_loader_train = DataLoader(
        dataset_train,
        batch_size=batch_size,
        num_workers=4,
        shuffle=True,
        pin_memory=True
    )

    dataset_valid=origin_dataset_train(transform=None,valid=True)
    data_loader_valid = DataLoader(
        dataset_valid,
        batch_size=batch_size,
        num_workers=4,
        shuffle=True,
        pin_memory=True
    )
    """
    dataset_train_basic=origin_dataset_train_basic(transform=None,valid=False)
    data_loader_train_basic = DataLoader(
        dataset_train_basic,
        batch_size=batch_size,
        num_workers=4,
        shuffle=True,
        pin_memory=True
    )
    """
    """
    #model_basic_group = CNN_basic(num_classes=21).to(device)
    model_basic_group= MotifNet3(input_size=12, hidden_size=64, num_classes=21).to(device) 
    #model_basic_group = ResNet18().to(device)#ComplexResidualMLP(layer_sizes=[4, 32, 36, 64, 104, 3], residual_layers=[1, 3]).to(device) #ResNet18()  
    #model_offset = TransformerModel().to(device) 
    #model = LSTMModel().to(device)  
    criterion = FocalLoss() #nn.CrossEntropyLoss()   
    #criterion_SupCon = SupConLoss2()
    optimizer_basic_group = torch.optim.Adam(model_basic_group.parameters()) 
    #optimizer_offset = torch.optim.Adam(model_offset.parameters(),lr=0.001) 
    #schedulr steplr 2-1 epoch *0.1 
    scheduler_basic_group = StepLR(optimizer_basic_group, step_size=2, gamma=0.1)
    #scheduler_basic_group = ReduceLROnPlateau(optimizer_basic_group, 'min')



    #scheduler_offset = StepLR(optimizer_offset, step_size=2, gamma=0.1)
    #scheduler = ReduceLROnPlateau(optimizer, mode='max', factor=args.lr_decay, patience=args.lr_patience, verbose=True)
    num_epochs = 50   
    loss_each_epoch=[]   
    for epoch in range(num_epochs):
        #if epoch%10==0 : print(epoch)
        model_basic_group.train()
        total_train_loss = []
        total_train_loss_SupCon = []
        for index,data in enumerate(data_loader_train):
            mean,base_num_fraction,std,time,t_test_pval,avg_basequality,g_num,g_value,basic,y_basic_group,y_offset=data
            mean=mean.to(device)
            std=std.to(device)
            time=time.to(device)
            t_test_pval=t_test_pval.to(device)
            avg_basequality=avg_basequality.to(device)
            g_num=g_num.to(device)
            g_value=g_value.to(device)
            basic=basic.to(device)
            base_num_fraction=base_num_fraction.to(device)
            y_basic_group=y_basic_group.to(device)
            y_offset=y_offset.to(device)
            y=y_basic_group*7+y_offset
            #outputs_basic_group=model_basic_group(mean,basic,sd,nat,wga)
            input = torch.cat((mean,std,time,base_num_fraction,avg_basequality), dim=1)
            #outputs_basic_group=model_basic_group(basic,input)
            #mean= mean.squeeze(1)
            outputs_basic_group=model_basic_group(input) #mean,base_num_fraction,g_num,avg_basequality,std,time
            #print(outputs_basic_group.shape,y.shape)
            #print(y)
            #print(torch.max(outputs.data, 1))
            #print(outputs.shape,y.shape)
            #loss_basic_group=criterion(outputs_basic_group, y_basic_group)  
            loss_basic_group=criterion(outputs_basic_group, y) 
            #loss_SupCon = criterion_SupCon(features=mid_mean, labels=y)
            loss=loss_basic_group#+0.1*loss_SupCon
            loss_value=loss.item()
            total_train_loss.append(loss_value)
            #total_train_loss_SupCon.append(loss_SupCon.item())
            #print(loss.item())
            optimizer_basic_group.zero_grad()  
            loss.backward()        
            optimizer_basic_group.step() 
        scheduler_basic_group.step()
     
        model_basic_group.eval()
        total_valid_loss=[]
        total_valid_loss_SupCon=[]
        with torch.no_grad():
            for index,data in enumerate(data_loader_valid):
                mean,base_num_fraction,std,time,t_test_pval,avg_basequality,g_num,g_value,basic,y_basic_group,y_offset=data
                mean=mean.to(device)
                std=std.to(device)
                time=time.to(device)
                t_test_pval=t_test_pval.to(device)
                avg_basequality=avg_basequality.to(device)
                g_num=g_num.to(device)
                g_value=g_value.to(device)
                basic=basic.to(device)
                base_num_fraction=base_num_fraction.to(device)
                #x=x.to(device)
                y_basic_group=y_basic_group.to(device)
                y_offset=y_offset.to(device)
                y=y_basic_group*7+y_offset
                input = torch.cat((mean,std,time,base_num_fraction,avg_basequality), dim=1)
                #outputs_basic_group=model_basic_group(basic,input)
                #mean= mean.squeeze(1)
                outputs_basic_group=model_basic_group(input)
                #print(outputs)
                #print(y)
                #print(torch.max(outputs.data, 1))
                #loss_basic_group=criterion(outputs_basic_group, y_basic_group)  
                loss_basic_group=criterion(outputs_basic_group,y) 
                #loss_SupCon = criterion_SupCon(features=mid_mean, labels=y)
                loss_valid=loss_basic_group#+0.1*loss_SupCon
                loss_value=loss_valid.item()
                total_valid_loss.append(loss_value)
                #total_valid_loss_SupCon.append(loss_SupCon.item())
        #scheduler_basic_group.step(np.mean(total_valid_loss))
        if loss_valid.item() < best_loss - tol:
            best_loss = loss_valid.item()
            loss_not_improved_count = 0
        else:
            loss_not_improved_count += 1

        if loss_not_improved_count >= n_iter_no_change:
            print(f"Training stopped at epoch {epoch+1}")
            torch.save(model_basic_group.state_dict(), f"/home/nipeng/chenzheng/mulan-methyl/my_model/mlp_{epoch+1}.pt")
            break

        if (epoch + 1) % 5 == 0:
            print(f'Epoch [{epoch + 1}/{num_epochs}],Train Loss: {np.mean(total_train_loss):.6f},Test Loss: {np.mean(total_valid_loss):.6f}')#,Train supCon Loss: {np.mean(total_train_loss_SupCon):.6f},Test supCon Loss: {np.mean(total_valid_loss_SupCon):.6f}')
            torch.save(model_basic_group.state_dict(), f"/home/nipeng/chenzheng/mulan-methyl/my_model/mlp_{epoch+1}.pt")
        #break
        if (epoch + 1) % 1 == 0:
            #print(f'Epoch [{epoch + 1}/{num_epochs}],Train Loss: {np.mean(total_train_loss):.6f}')#,Test Loss: {np.mean(total_valid_loss):.6f})
            torch.save(model_basic_group.state_dict(), f"/home/nipeng/chenzheng/mulan-methyl/my_model/mlp_{epoch+1}.pt")
        #break
    """
    #model = ResNet18().to(device)  
    #model = CNN_basic(num_classes=21).to(device)
    model = MotifNet3(input_size=12, hidden_size=64, num_classes=21).to(device)
    #model = joblib.load('/home/nipeng/chenzheng/mulan-methyl/saved_model/MLPClassifier.pkl')
    #model_basic = MLP_basic(input_size=22, hidden_size=22, num_classes=3).to(device)
    model.load_state_dict(torch.load("/home/nipeng/chenzheng/mulan-methyl/my_model/mlp_7.pt")) #mlp_29_data_mean_12_code_filter_iso_offset_21_avg.pt 上采样:mlp_30_data_mean_12_code_filter_iso_offset_21_upsample 降采样:mlp_29_data_mean_12_code_filter_iso_offset_21_avg.pt mlp_46_data_mean_12_code_filter_iso_offset_21_avg.pt
    #model_basic.load_state_dict(torch.load("/home/nipeng/chenzheng/mulan-methyl/my_model/mlp_23_data_mean_12_code_filter_iso_offset_3.pt")) #mlp_15_data_mean_12_code_filter_iso_offset_3_avg #NG:21,MH:27
    #for name, param in model.state_dict().items():
    #    print(name, param)
    #    break
    TP_motif=["CGCG_-14","CCTCC_-14","CAGAAA_-14","CCCRAG_-14","CTACT_-14","GATC_-14","GGNCC_-14","RAACTC_-14"]
    #TP_motif=["CGCG_-14","CCTCC_-14","CTACT_-14","GGNCC_-14","RAACTC_-14"]
    NO_motif=["CAANNNNNNNCTGG_-14","CCAGNNNNNNNTTG_-14","CTCGAG_-14","GCGGCCGC_-14"]
    NG_motif=["CCGCGG_-14","GCCGGC_-14","GAGNNNNNTAC_-14","GCANNNNNNNNTGC_-14","GGCC_-14","GGNNCC_-14","GGTGA_-14","GTANNNNNCTC_-14","RGCGCY_-14"]
    MH_motif=["CTNAG_-14","AGCT_-14","CCACGK_-14","GATC_-14","GCYYGAT_-14","GTAC_-14"]
    HP_motif=["CCGG_-14","ATTAAT_-14","CATG_-14","CRTANNNNNNNWC_-14","CSAG_-14","CTRYAG_-14","CYANNNNNNTTC_-14","GCGC_-14","GAGG_-14","GANNNNNNNTAYG_-14","GAANNNNNNTRG_-14","GAATTC_-14","GGCC_-14","GMRGA_-14","GTAC_-14","GTNNAC_-14","TCTTC_-14","TCGA_-14","TCNNGA_-14","TGCA_-14"]
    EC_motif=["AACNNNNNNGTGC_-14","CCWGG_-14","GATC_-14","GCACNNNNNNGTT_-14"]
    CP_motif=["CCGG_-14","CACNNNNNRTAAA_-14","GATC_-14","GGWCC_-14","GTATAC_-14","TTTAYNNNNNGTG_-14","VGACAT_-14"]
    BF_motif=["AACNNNNNNGTGC_-14","CCWGG_-14","GATC_-14","GCACNNNNNNGTT_-14"]
    true_label_TP=[10,11,17,18,10,17,10,10]
    true_label_TP_offset=[3,4,3,4,3,3,3,3]
    #true_label_TP_offset=[3,4,3,3,3]
    true_label_TP_basic_group=[1,1,2,2,1,2,1,1]
    true_label_NO=[17,17,4,10]
    true_label_NO_basic_group=[2,2,0,1]
    true_label_NG_basic_group=[0,0,2,2,0,0,2,2,0]
    true_label_MH_basic_group=[1,1,1,2,2,1]
    true_label_HP_basic_group=[1,2,2,2,2,2,2,0,2,2,2,2,0,2,2,2,1,2,2,2]
    true_label_EC_basic_group=[2,0,2,2]
    true_label_CP_basic_group=[0,2,0,0,2,2,2]
    MEAN_total=0
    MEAN=0
    for i,motif_item in enumerate(NO_motif):
        y_score=[]
        y_basic_score=[]
        predict_result=[]
        predict_basic_result=[]
        dataset_test=origin_dataset_test(transform=None,motif=motif_item,species="NO")
        data_loader_test = DataLoader(
            dataset_test,
            batch_size=32,
            num_workers=1,
            shuffle=False,
            pin_memory=True,
            drop_last=False,
        )
        model.eval()
        #model_basic.eval()
        with torch.no_grad():
            for index,Data in enumerate(data_loader_test):
                data,data_b=Data
                mean,base_num_fraction,std,time,t_test_pval,avg_basequality,g_num,g_value,basic=data
                mean_b,base_num_fraction_b,avg_basequality_b,std_b,time_b=data_b
                
                base_num_fraction=base_num_fraction.to(device)
                std=std.to(device)
                time=time.to(device)
                avg_basequality=avg_basequality.to(device)
                g_num=g_num.to(device)
                t_test_pval=t_test_pval.to(device)
                basic=basic.to(device)

                mean_b=mean_b.to(device)
                avg_basequality_b=avg_basequality_b.to(device)
                base_num_fraction_b=base_num_fraction_b.to(device)
                std_b=std_b.to(device)
                time_b=time_b.to(device)
                #outputs=model.predict_proba(mean)#,base_num_fraction,g_num,avg_basequality,std)
                mean=mean.to(device)
                input = torch.cat((mean,std,time,base_num_fraction,avg_basequality), dim=1)
                #utputs=model(basic,input)
                #mean= mean.squeeze(1)
                outputs=model(input)
                #outputs_basic=model_basic(mean,base_num_fraction,g_num,avg_basequality,std,time,t_test_pval)
                #outputs_basic=model_basic(mean_b,base_num_fraction_b,g_num,avg_basequality_b,std_b,time_b)
                _, predicted = torch.max(outputs.data, 1)
                #_, predicted_basic = torch.max(outputs_basic.data, 1)
                predict_prob=outputs.data
                #predict_basic_prob=outputs_basic.data

                #_, predicted = torch.max(outputs.data, 1)
                predicted=predicted.cpu().numpy().tolist()
                predict_prob=predict_prob.cpu().numpy().tolist()
                #predict_basic_prob=predict_basic_prob.cpu().numpy().tolist()
                #predicted_basic=predicted_basic.cpu().numpy().tolist()
                #predicted=outputs.tolist()
                predict_result+=predicted
                #predict_basic_result+=predicted_basic
                y_score+=predict_prob
                #y_basic_score+=predict_basic_prob

        #y_test=np.array([true_label_NO[i]]*len(dataset_test))
        #y_test = label_binarize(y_test, classes=range(21))
        y_score=np.array(y_score)
        #y_basic_score=np.array(y_basic_score)
        #predict_basic_result=np.array(predict_basic_result)
        predict_result=np.array(predict_result)
    
        result = pd.Series(predict_result).value_counts()
        #counts = np.bincount(predict_result)
        #most_num = np.argmax(counts)
        most_num=result.idxmax()
        accuracy=int(round(result[true_label_NO[i]]/len(dataset_test),2)*100)
        MEAN+=accuracy
        MEAN_total+=accuracy
        print(accuracy)

    print(round(MEAN/len(true_label_NO),2))

    MEAN=0
    for i,motif_item in enumerate(TP_motif):
        y_score=[]
        y_basic_score=[]
        predict_result=[]
        predict_basic_result=[]
        dataset_test=origin_dataset_test(transform=None,motif=motif_item,species="TP")
        data_loader_test = DataLoader(
            dataset_test,
            batch_size=32,
            num_workers=1,
            shuffle=False,
            pin_memory=True,
            drop_last=False,
        )
        model.eval()
        #model_basic.eval()
        with torch.no_grad():
            for index,Data in enumerate(data_loader_test):
                data,data_b=Data
                mean,base_num_fraction,std,time,t_test_pval,avg_basequality,g_num,g_value,basic=data
                mean_b,base_num_fraction_b,avg_basequality_b,std_b,time_b=data_b
                
                base_num_fraction=base_num_fraction.to(device)
                std=std.to(device)
                time=time.to(device)
                avg_basequality=avg_basequality.to(device)
                g_num=g_num.to(device)
                t_test_pval=t_test_pval.to(device)
                basic=basic.to(device)

                mean_b=mean_b.to(device)
                avg_basequality_b=avg_basequality_b.to(device)
                base_num_fraction_b=base_num_fraction_b.to(device)
                std_b=std_b.to(device)
                time_b=time_b.to(device)
                #outputs=model.predict_proba(mean)#,base_num_fraction,g_num,avg_basequality,std)
                mean=mean.to(device)
                input = torch.cat((mean,std,time,base_num_fraction,avg_basequality), dim=1)
                #utputs=model(basic,input)
                #mean= mean.squeeze(1)
                outputs=model(input)
                #outputs_basic=model_basic(mean,base_num_fraction,g_num,avg_basequality,std,time,t_test_pval)
                #outputs_basic=model_basic(mean_b,base_num_fraction_b,g_num,avg_basequality_b,std_b,time_b)
                _, predicted = torch.max(outputs.data, 1)
                #_, predicted_basic = torch.max(outputs_basic.data, 1)
                predict_prob=outputs.data
                #predict_basic_prob=outputs_basic.data

                #_, predicted = torch.max(outputs.data, 1)
                predicted=predicted.cpu().numpy().tolist()
                predict_prob=predict_prob.cpu().numpy().tolist()
                #predict_basic_prob=predict_basic_prob.cpu().numpy().tolist()
                #predicted_basic=predicted_basic.cpu().numpy().tolist()
                #predicted=outputs.tolist()
                predict_result+=predicted
                #predict_basic_result+=predicted_basic
                y_score+=predict_prob
                #y_basic_score+=predict_basic_prob

        #y_test=np.array([true_label_NO[i]]*len(dataset_test))
        #y_test = label_binarize(y_test, classes=range(21))
        y_score=np.array(y_score)
        #y_basic_score=np.array(y_basic_score)
        #predict_basic_result=np.array(predict_basic_result)
        predict_result=np.array(predict_result)
    
        result = pd.Series(predict_result).value_counts()
        #counts = np.bincount(predict_result)
        #most_num = np.argmax(counts)
        most_num=result.idxmax()
        accuracy=int(round(result[true_label_TP[i]]/len(dataset_test),2)*100)
        MEAN+=accuracy
        MEAN_total+=accuracy
        print(accuracy)

    #print(MEAN/len(true_label_TP))
    print(round(MEAN/len(true_label_TP),2))
    print(round(MEAN_total/12,2))
    #"""