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
from torch.optim.lr_scheduler import StepLR
import os
from torch.nn import functional as F
from sklearn.preprocessing import label_binarize
import matplotlib.pyplot as plt
os.environ['CUDA_VISIBLE_DEVICES'] = "0"
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
        self.po = nn.AdaptiveMaxPool1d(1)
        """
        self.layers1 = nn.Sequential(MyConvBlock(i,c,3),nn.BatchNorm1d(c),nn.ReLU(),MyConvBlock(c,c,2),nn.BatchNorm1d(c),nn.ReLU())
        self.layers11 = nn.Sequential(MyConvBlock(i,c,2),nn.BatchNorm1d(c),nn.ReLU(),MyConvBlock(c,c,2),nn.BatchNorm1d(c),nn.ReLU(),MyConvBlock(c,c,2),nn.BatchNorm1d(c),nn.ReLU())
        self.layers111 = nn.Sequential(MyConvBlock(i,c,4),nn.BatchNorm1d(c),nn.ReLU())
        self.layers1111 = nn.Sequential(nn.Conv1d(3*c, c, kernel_size=1),nn.BatchNorm1d(c))
        self.cls = nn.Sequential(nn.ReLU(),Flatten(),nn.Linear(9*c, num_classes))
        """

    def forward(self,mean,base_num_fraction,g_num,avg_basequality,std,time):

        mean_1 = self.layers1(mean)
        mean_2 = self.layers11(mean)
        mean_3 = self.layers111(mean)
        mean_cat=torch.cat((mean_1,mean_2,mean_3), dim=1)
        mean_cat=mean_cat*F.softmax(self.po(mean_cat),dim=1)
        mean4 = self.layers1111(mean_cat)
        mean4=mean4+(mean_1+mean_2+mean_3)
        mean4=self.cls(mean4)

        return mean4
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
    def forward(self,mean,base_num_fraction,g_num,avg_basequality,std,time):

        mean_1 = self.layers1(mean)
        """
        base_num_fraction_1 = self.layers2(base_num_fraction)
        avg_basequality_1 = self.layers3(avg_basequality)
        std_1 = self.layers4(std)
        time_1 = self.layers5(time)
        """
        return mean_1
class Flatten(nn.Module):
    def forward(self, input):
        return input.view(input.size(0), -1)
class MyConvBlock(nn.Module):
    def __init__(self, in_channels, out_channels, kernel_size, stride=1):
        super(MyConvBlock, self).__init__()
        self.conv1 = nn.Conv1d(in_channels, in_channels, kernel_size, stride=stride, groups=in_channels)
        self.bn1 = nn.BatchNorm1d(in_channels)
        self.relu1 = nn.ReLU()
        self.conv2 = nn.Conv1d(in_channels, out_channels, 1)
        self.bn2 = nn.BatchNorm1d(out_channels)
        self.relu2 = nn.ReLU()

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

        c=4
        i=5
        self.layers1 = nn.Sequential(MyConvBlock(i,c,6),nn.BatchNorm1d(c),nn.ReLU())
        self.layers11 = nn.Sequential(MyConvBlock(i,c,5),nn.BatchNorm1d(c),nn.ReLU(),nn.Linear(8, 7),nn.BatchNorm1d(c),nn.ReLU())
        self.layers111 = nn.Sequential(MyConvBlock(i,c,7),nn.BatchNorm1d(c),nn.ReLU(),nn.Linear(6, 7),nn.BatchNorm1d(c),nn.ReLU())
        self.layers1111 = nn.Sequential(nn.Conv1d(3*c, c, kernel_size=1),nn.BatchNorm1d(c))
        self.cls = nn.Sequential(nn.ReLU(),Flatten(),nn.Linear(7*c, num_classes))


    def forward(self,mean,base_num_fraction,g_num,avg_basequality,std,time):

        mean_1 = self.layers1(mean)
        mean_2 = self.layers11(mean)
        mean_3 = self.layers111(mean)
        mean4 = self.layers1111(torch.cat((mean_1,mean_2,mean_3), dim=1))
        mean4+=mean_1
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
    def __init__(self, transform,valid,motif):
        super(origin_dataset_train, self).__init__()
        if valid==False:     
            path_data = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_mean_{}_train.pkl'.format(motif)
            path_data_base_num_fraction = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_base_num_fraction_{}_train.pkl".format(motif)
            path_data_std_diff = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_std_diff_{}_train.pkl".format(motif)
            path_data_time_diff = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_time_diff_{}_train.pkl".format(motif)
            path_test_g_num = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_g_num_{}_train.pkl".format(motif)
            path_test_avg_basequality = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_avg_basequality_{}_train.pkl".format(motif)
            path_y = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_{}_train.pkl'.format(motif)
        else:
            path_data = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_mean_{}_valid.pkl'.format(motif)
            path_data_base_num_fraction = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_base_num_fraction_{}_valid.pkl".format(motif)
            path_data_std_diff = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_std_diff_{}_valid.pkl".format(motif)
            path_data_time_diff = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_time_diff_{}_valid.pkl".format(motif)
            path_test_g_num = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_g_num_{}_valid.pkl".format(motif)
            path_test_avg_basequality = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_avg_basequality_{}_valid.pkl".format(motif)
            path_y = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_{}_valid.pkl'.format(motif)

        self.data=joblib.load(path_data) 
        self.data_base_num_fraction=joblib.load(path_data_base_num_fraction) 
        self.data_std_diff=joblib.load(path_data_std_diff)
        self.data_time_diff=joblib.load(path_data_time_diff)
        self.data_g_num=joblib.load(path_test_g_num)
        self.data_avg_basequality=joblib.load(path_test_avg_basequality)
        self.label=joblib.load(path_y)  

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
        self.data=torch.FloatTensor(self.data)
        self.data_base_num_fraction=torch.FloatTensor(self.data_base_num_fraction)
        self.data_std_diff=torch.FloatTensor(self.data_std_diff)
        self.data_time_diff=torch.FloatTensor(self.data_time_diff)
        self.data_avg_basequality=torch.FloatTensor(self.data_avg_basequality)
        self.data_g_num=torch.FloatTensor(self.data_g_num)
        self.label=torch.LongTensor(self.label)

        self.data = self.data.unsqueeze(1)
        self.data_base_num_fraction = self.data_base_num_fraction.unsqueeze(1)
        self.data_std_diff = self.data_std_diff.unsqueeze(1)
        self.data_time_diff = self.data_time_diff.unsqueeze(1)
        self.data_avg_basequality = self.data_avg_basequality.unsqueeze(1)
    def __len__(self):
        return self.len_data
 
    def __getitem__(self, index):
        #data = torch.FloatTensor(data)
        return  self.data[index],self.data_base_num_fraction[index],self.data_std_diff[index],self.data_time_diff[index],self.data_avg_basequality[index],self.data_g_num[index],self.label[index]
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
        #path_test_g_value = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/{}_pd_from_r_df_test_g_value_12_code_filter_iso.pkl".format(species)

        pd_from_r_df_test=joblib.load(path_test)   
        pd_from_r_df_test_base_num_fraction=joblib.load(path_test_base_num_fraction)
        pd_from_r_df_test_std_diff=joblib.load(path_test_std_diff)
        pd_from_r_df_test_time_diff=joblib.load(path_test_time_diff)
        pd_from_r_df_test_basic_group=joblib.load(path_test_basic_group)
        pd_from_r_df_test_t_test_pval=joblib.load(path_test_t_test_pval)
        pd_from_r_df_test_avg_basequality=joblib.load(path_test_avg_basequality)
        pd_from_r_df_test_g_num=joblib.load(path_test_g_num)
        #pd_from_r_df_test_g_value=joblib.load(path_test_g_value)

        subet_test_data=pd_from_r_df_test[pd_from_r_df_test["label"]==motif]
        subet_test_data_base_num_fraction=pd_from_r_df_test_base_num_fraction[pd_from_r_df_test_base_num_fraction["label"]==motif]
        subet_test_data_std_diff=pd_from_r_df_test_std_diff[pd_from_r_df_test_std_diff["label"]==motif]
        subet_test_data_time_diff=pd_from_r_df_test_time_diff[pd_from_r_df_test_time_diff["label"]==motif]
        subet_test_data_basic_group=pd_from_r_df_test_basic_group[pd_from_r_df_test_basic_group["label"]==motif]
        subet_test_data_t_test_pval=pd_from_r_df_test_t_test_pval[pd_from_r_df_test_t_test_pval["label"]==motif]
        subet_test_data_avg_basequality=pd_from_r_df_test_avg_basequality[pd_from_r_df_test_avg_basequality["label"]==motif]
        subet_test_data_g_num=pd_from_r_df_test_g_num[pd_from_r_df_test_g_num["label"]==motif]
        #subet_test_data_g_value=pd_from_r_df_test_g_value[pd_from_r_df_test_g_value["label"]==motif]

        self.data=subet_test_data.values[:,6:18] #6:18
        self.data=self.data.astype(float) 
        self.data = normalize_test(self.data)

        self.data_base_num_fraction=subet_test_data_base_num_fraction.values[:,6:18] #6:18
        self.data_base_num_fraction=self.data_base_num_fraction.astype(float) 
        self.data_base_num_fraction = normalize_test(self.data_base_num_fraction,base_num_fraction=True)

        self.data_std_diff=subet_test_data_std_diff.values[:,6:18] #6:18
        self.data_std_diff=self.data_std_diff.astype(float) 
        self.data_std_diff = normalize_test(self.data_std_diff,std=True)

        self.data_time_diff=subet_test_data_time_diff.values[:,6:18] #6:18
        self.data_time_diff=self.data_time_diff.astype(float)
        self.data_time_diff = normalize_test(self.data_time_diff,time=True)

        self.data_avg_basequality=subet_test_data_avg_basequality.values[:,6:18] #6:18
        self.data_avg_basequality=self.data_avg_basequality.astype(float)
        self.data_avg_basequality = normalize_test(self.data_avg_basequality,avg_basequality=True)

        self.data_g_num=subet_test_data_g_num.values[:,6:18] #6:18
        self.data_g_num=self.data_g_num.astype(float)
        self.data_g_num = normalize_test(self.data_g_num,g_num=True)

        #self.data = umap_model.transform(self.data)
        self.data_basic_group=subet_test_data_basic_group.values[:,6:18] #1:23
        self.data_basic_group[self.data_basic_group=='A']=0
        self.data_basic_group[self.data_basic_group=='C']=1
        self.data_basic_group[self.data_basic_group=='T']=2
        self.data_basic_group[self.data_basic_group=='G']=3
        self.data_basic_group = self.data_basic_group.astype(int) 
        self.data_t_test_pval=subet_test_data_t_test_pval.values[:,6:18] #1:23
        self.data_t_test_pval=self.data_t_test_pval.astype(float)
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
         return  self.data[index],self.data_base_num_fraction[index],self.data_std_diff[index],self.data_time_diff[index],self.data_t_test_pval[index],self.data_avg_basequality[index],self.data_g_num[index],self.data_basic_group[index]
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
    tol=1e-4
    n_iter_no_change=10
    loss_not_improved_count = 0
    best_loss = float('inf')
    #scaler = StandardScaler() 
    motifs_test= {}
    """
    motifs_test["TP"]=["CGCG_-14","CCTCC_-14","CAGAAA_-14","CCCRAG_-14","CTACT_-14","GATC_-14","GGNCC_-14","RAACTC_-14"]
    motifs_test["NO"]=["CAANNNNNNNCTGG_-14","CCAGNNNNNNNTTG_-14","CTCGAG_-14","GCGGCCGC_-14"]
    motifs_test["NG"]=["CCGCGG_-14","GCCGGC_-14","GAGNNNNNTAC_-14","GCANNNNNNNNTGC_-14","GGCC_-14","GGNNCC_-14","GGTGA_-14","GTANNNNNCTC_-14","RGCGCY_-14"]
    motifs_test["MH"]=["CTNAG_-14","AGCT_-14","CCACGK_-14","GATC_-14","GCYYGAT_-14","GTAC_-14"]
    motifs_test["HP"]=["CCGG_-14","ATTAAT_-14","CATG_-14","CRTANNNNNNNWC_-14","CSAG_-14","CTRYAG_-14","CYANNNNNNTTC_-14","GCGC_-14","GAGG_-14","GANNNNNNNTAYG_-14","GAANNNNNNTRG_-14","GAATTC_-14","GGCC_-14","GMRGA_-14","GTAC_-14","GTNNAC_-14","TCTTC_-14","TCGA_-14","TCNNGA_-14","TGCA_-14"]
    motifs_test["EC"]=["AACNNNNNNGTGC_-14","CCWGG_-14","GATC_-14","GCACNNNNNNGTT_-14"]
    motifs_test["CP"]=["CCGG_-14","CACNNNNNRTAAA_-14","GATC_-14","GGWCC_-14","GTATAC_-14","TTTAYNNNNNGTG_-14","VGACAT_-14"]
    motifs_test["BF"]=["GATC_-14"]
    motifs_test["BA"]=["GCWGC_-14","GGATCC_-14"]
    """
    motifs_test["TP"]=["CGCG_-14","CCTCC_-14","CAGAAA_-14","CCCRAG_-14","CTACT_-14","GATC_-14","GGNCC_-14","RAACTC_-14"]
    motifs_test["NO"]=["CAANNNNNNNCTGG_-14","CCAGNNNNNNNTTG_-14","CTCGAG_-14","GCGGCCGC_-14"]

    motifs_test["NG"]=["C5mCGCGG","G5mCCGGC","G6mAGNNNNNTAC","GC6mANNNNNNNNTGC","GG5mCC","GGNN5mCC","GGTG6mA","GT6mANNNNNCTC","RG5mCGCY"]
    motifs_test["MH"]=["4mCTNAG","AG4mCT","CCA4mCGK","G6mATC","GCYYG6mAT","GTA4mC"]
    motifs_test["HP"]=["4mCCGG","ATTA6mAT","C6mATG","CRT6mANNNNNNNWC","CS6mAG","CTRY6mAG","CY6mANNNNNNTTC","G5mCGC","G6mAGG","G6mANNNNNNNTAYG","GA6mANNNNNNTRG","GA6mATTC","GG5mCC","GMRG6mA","GT6mAC","GTNN6mAC","T4mCTTC","TCG6mA","TCNNG6mA","TGC6mA"]
    motifs_test["EC"]=["A6mACNNNNNNGTGC","C5mCWGG","G6mATC","GC6mACNNNNNNGTT"]
    motifs_test["CP"]=["5mCCGG", "C6mACNNNNNRTAAA","GAT5mC","GGW5mCC","GTAT6mAC","TTT6mAYNNNNNGTG","VGAC6mAT"]
    motifs_test["BF"]=["GAT5mC"]
    motifs_test["BA"]=["G5mCWGC", "GGAT4mCC"]
    #map_list={"GATC_-14":"GAT5mC","CCGG_-14":"5mCCGG","CACNNNNNRTAAA_-14":"C6mACNNNNNRTAAA","GGWCC_-14":"GGW5mCC","GTATAC_-14":"GTAT6mAC","TTTAYNNNNNGTG_-14":"TTT6mAYNNNNNGTG","VGACAT_-14":"VGAC6mAT"}
    map_list={"G5mCWGC":"GCWGC_-14", "GGAT4mCC":"GGATCC_-14", 
                "GAT5mC":"GATC_-14", 
                "5mCCGG":"CCGG_-14","C6mACNNNNNRTAAA":"CACNNNNNRTAAA_-14", "GGW5mCC":"GGWCC_-14","GTAT6mAC":"GTATAC_-14","TTT6mAYNNNNNGTG":"TTTAYNNNNNGTG_-14","VGAC6mAT":"VGACAT_-14",
                "A6mACNNNNNNGTGC":"AACNNNNNNGTGC_-14","C5mCWGG":"CCWGG_-14","G6mATC":"GATC_-14","GC6mACNNNNNNGTT":"GCACNNNNNNGTT_-14",
                "4mCCGG":"CCGG_-14","ATTA6mAT":"ATTAAT_-14","C6mATG":"CATG_-14","CRT6mANNNNNNNWC":"CRTANNNNNNNWC_-14","CS6mAG":"CSAG_-14","CTRY6mAG":"CTRYAG_-14","CY6mANNNNNNTTC":"CYANNNNNNTTC_-14","G5mCGC":"GCGC_-14","G6mAGG":"GAGG_-14","G6mANNNNNNNTAYG":"GANNNNNNNTAYG_-14","GA6mANNNNNNTRG":"GAANNNNNNTRG_-14","GA6mATTC":"GAATTC_-14","GMRG6mA":"GMRGA_-14","GT6mAC":"GTAC_-14","GTNN6mAC":"GTNNAC_-14","T4mCTTC":"TCTTC_-14","TCG6mA":"TCNNGA_-14","TCNNG6mA":"TCNNGA_-14","TGC6mA":"TGCA_-14",
                "4mCTNAG":"CTNAG_-14","AG4mCT":"AGCT_-14","CCA4mCGK":"CCACGK_-14","GCYYG6mAT":"GCYYGAT_-14" ,"GTA4mC":"GTAC_-14",
                "C5mCGCGG":"CCGCGG_-14","G5mCCGGC":"GCCGGC_-14","G6mAGNNNNNTAC":"GAGNNNNNTAC_-14","GC6mANNNNNNNNTGC":"GCANNNNNNNNTGC_-14","GG5mCC":"GGCC_-14","GGNN5mCC":"GGNNCC_-14","GGTG6mA":"GGTGA_-14","GT6mANNNNNCTC":"GTANNNNNCTC_-14","RG5mCGCY":"RGCGCY_-14"}
    #true_label_TP_offset=[3,4,3,4,3,3,3,3]
    #true_label_TP_basic_group=[1,1,2,2,1,2,1,1]
    #true_label_NO_basic_group=[2,2,0,1]
    true_label={}
    true_label["TP"]=[10,11,17,18,10,17,10,10]
    true_label["NO"]=[17,17,4,10]
    true_label["NG"]=[4,3,18,17,3,3,17,17,3]
    true_label["MH"]=[10,10,11,17,18,10]
    true_label["HP"]=[10,18,18,17,17,17,17,4,17,17,17,18,3,17,17,17,11,17,17,17]
    true_label["EC"]=[17,4,17,17]
    true_label["CP"]=[3,17,3,4,17,18,18]
    true_label["BF"]=[3]
    true_label["BA"]=[4,8]
    """
    true_label_basic_group={}
    true_label_offset={}
    true_label_basic_group["NG"]=[0,0,2,2,0,0,2,2,0]
    true_label_basic_group["MH"]=[1,1,1,2,2,1]
    true_label_basic_group["HP"]=[1,2,2,2,2,2,2,0,2,2,2,2,0,2,2,2,1,2,2,2]
    true_label_basic_group["EC"]=[2,0,2,2]
    true_label_basic_group["CP"]=[0,2,0,0,2,2,2]
    true_label_basic_group["BF"]=[0]
    true_label_basic_group["BA"]=[0,1]

    true_label_offset["NG"]=[4,3,4,3,3,3,3,3,3]
    true_label_offset["MH"]=[3,3,4,3,4,3]
    true_label_offset["HP"]=[3,4,4,3,3,3,3,4,3,3,3,4,3,3,3,3,4,3,3,3]
    true_label_offset["EC"]=[3,4,3,3]
    true_label_offset["CP"]=[3,3,3,4,3,4,4]
    true_label_offset["BF"]=[3]
    true_label_offset["BA"]=[4,1]
    """
    bacterials_set=["BA","BF","CP","EC","HP","MH","NG"] #["BA","BF","CP","EC","HP","MH","NG"]
    #bacterials_set=["HP"]
    #motifs_test["HP"]=["G6mAGG"]
    #motifs_test["HP"]=["T4mCTTC"]
    #true_label["HP"]=[17]
    #true_label["HP"]=[11]
    model = MotifNet3(input_size=12, hidden_size=64, num_classes=21).to(device)
    MEAN_total=0
    motif_num={}
    motif_acc={}
    motif_acc_final={}
    motif_num["GAT5mC"]=[17300,7664]
    motif_num["G6mATC"]=[28187,27440]
    motif_num["GG5mCC"]=[829,6508]
    motif_acc["GAT5mC"]=[]
    motif_acc["G6mATC"]=[]
    motif_acc["GG5mCC"]=[]
    with open('log_test_new.txt', 'a') as f:
        for bacterial in bacterials_set:
            MEAN=0
            #print(bacterial,file=f)
            for i,motif_item in enumerate(motifs_test[bacterial]):
                #if motif_item=="GAT5mC":
                print(motif_item,file=f)
                model.load_state_dict(torch.load(f"/home/nipeng/chenzheng/mulan-methyl/my_model/mlp_mean_{motif_item}.pt"))
                y_score=[]
                y_basic_score=[]
                predict_result=[]
                predict_basic_result=[]
                dataset_train=origin_dataset_train(transform=None,valid=False,motif=motif_item)
                data_loader_train = DataLoader(
                    dataset_train,
                    batch_size=batch_size,
                    num_workers=4,
                    shuffle=True,
                    pin_memory=True
                )
                dataset_test=origin_dataset_test(transform=None,motif=map_list[motif_item],species=bacterial)
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
                        mean,base_num_fraction,std,time,t_test_pval,avg_basequality,g_num,basic=Data
                        
                        base_num_fraction=base_num_fraction.to(device)
                        std=std.to(device)
                        time=time.to(device)
                        avg_basequality=avg_basequality.to(device)
                        g_num=g_num.to(device)
                        t_test_pval=t_test_pval.to(device)
                        basic=basic.to(device)
                        mean=mean.to(device)

                        input = torch.cat((mean,std,time,base_num_fraction,avg_basequality), dim=1)
                        outputs=model(input,base_num_fraction,g_num,avg_basequality,std,time)
                        _, predicted = torch.max(outputs.data, 1)
                        #_, predicted_basic = torch.max(outputs_basic.data, 1)
                        predict_prob=outputs.data
                        predicted=predicted.cpu().numpy().tolist()
                        predict_prob=predict_prob.cpu().numpy().tolist()
                        predict_result+=predicted
                        #predict_basic_result+=predicted_basic
                        y_score+=predict_prob
                        #y_basic_score+=predict_basic_prob

                y_score=np.array(y_score)
                predict_result=np.array(predict_result)
            
                result = pd.Series(predict_result).value_counts()
                print(result,file=f)
                most_num=result.idxmax()
                accuracy=int(round(result[true_label[bacterial][i]]/len(dataset_test),2)*100)
                MEAN+=accuracy
                #print(accuracy,file=f)
                if (bacterial=="BF" or bacterial=="CP") and motif_item=="GAT5mC":
                    motif_acc["GAT5mC"].append(accuracy)
                elif (bacterial=="EC" or bacterial=="MH") and motif_item=="G6mATC":
                    motif_acc["G6mATC"].append(accuracy)
                elif (bacterial=="HP" or bacterial=="NG") and motif_item=="GG5mCC":
                    motif_acc["GG5mCC"].append(accuracy)
                else:
                    MEAN_total+=accuracy
                print(map_list[motif_item]," ",accuracy," ",len(dataset_test),file=f)

            #print(round(MEAN/len(true_label[bacterial])),file=f)
        #print(round(MEAN_total/49,2),file=f)
        motif_twice=["GAT5mC","G6mATC","GG5mCC"]
        print(motif_acc)
        for item in motif_twice:
            motif_acc_final[item]=(motif_acc[item][0]*motif_num[item][0]+motif_acc[item][1]*motif_num[item][1])/(motif_num[item][0]+motif_num[item][1])
            MEAN_total+=motif_acc_final[item]
        print(round(MEAN_total/46,2),file=f)
        print(motif_acc_final,file=f)

