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
class SELayer(nn.Module):
    def __init__(self, channel, reduction=16):
        super(SELayer, self).__init__()
        self.avg_pool = nn.AdaptiveAvgPool2d(1)
        self.fc = nn.Sequential(
            nn.Linear(channel, channel // reduction, bias=False),
            nn.ReLU(inplace=True),
            nn.Linear(channel, channel, bias=False),
            nn.Sigmoid()
        )

    def forward(self, x):
        b, c, _, _ = x.size()
        y = self.avg_pool(x).view(b, c)
        y = self.fc(y).view(b, c, 1, 1)
        return x * y.expand_as(x)

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
    #scaler = StandardScaler() 
    motifs_set1 = ["G5mCWGC","GGAT4mCC","GAT5mC","5mCCGG", "C6mACNNNNNRTAAA","GGW5mCC","GTAT6mAC","TTT6mAYNNNNNGTG","VGAC6mAT","A6mACNNNNNNGTGC"]
    motifs_set2 = ["C5mCWGG","G6mATC","GC6mACNNNNNNGTT","4mCCGG","ATTA6mAT","C6mATG","CRT6mANNNNNNNWC","CS6mAG","CTRY6mAG","CY6mANNNNNNTTC","G5mCGC","G6mAGG","G6mANNNNNNNTAYG"]
    motifs_set3 = ["GA6mANNNNNNTRG","GA6mATTC","GMRG6mA","GT6mAC","GTNN6mAC","T4mCTTC","TCG6mA","TCNNG6mA","TGC6mA","4mCTNAG","AG4mCT","CCA4mCGK","GCYYG6mAT","GTA4mC"]
    motifs_set4 = ["C5mCGCGG","G5mCCGGC","G6mAGNNNNNTAC","GC6mANNNNNNNNTGC","GG5mCC","GGNN5mCC","GGTG6mA","GT6mANNNNNCTC","RG5mCGCY"]
    #done_set = ["G5mCWGC", "GGAT4mCC","GAT5mC","5mCCGG", "C6mACNNNNNRTAAA","GAT5mC","GGW5mCC","GTAT6mAC","TTT6mAYNNNNNGTG","VGAC6mAT","A6mACNNNNNNGTGC","C5mCWGG","G6mATC","GC6mACNNNNNNGTT",
    #                "4mCCGG","ATTA6mAT","C6mATG","CRT6mANNNNNNNWC","CS6mAG","CTRY6mAG","CY6mANNNNNNTTC","G5mCGC","G6mAGG","G6mANNNNNNNTAYG","GA6mANNNNNNTRG","GA6mATTC","GG5mCC"]
    #motifs_set=["GAT5mC","GAT5mC"]
    print(len(motifs_set1)+len(motifs_set2)+len(motifs_set3)+len(motifs_set4))
    with open('log4.txt', 'a') as f:
        for motif in motifs_set4:
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
            dataset_train=origin_dataset_train(transform=None,valid=False,motif=motif)
            data_loader_train = DataLoader(
                dataset_train,
                batch_size=batch_size,
                num_workers=4,
                shuffle=True,
                pin_memory=True
            )

            dataset_valid=origin_dataset_train(transform=None,valid=True,motif=motif)
            data_loader_valid = DataLoader(
                dataset_valid,
                batch_size=batch_size,
                num_workers=4,
                shuffle=True,
                pin_memory=True
            )
            print(motif,file=f)
            model_basic_group= MotifNet3(input_size=12, hidden_size=64, num_classes=21).to(device)   
            criterion = FocalLoss()  #nn.CrossEntropyLoss() 
            optimizer_basic_group = torch.optim.Adam(model_basic_group.parameters()) 
            scheduler_basic_group = StepLR(optimizer_basic_group, step_size=2, gamma=0.1)
            num_epochs = 50  
            loss_each_epoch=[]   
            for epoch in range(num_epochs):
                #if epoch%10==0 : print(epoch)
                model_basic_group.train()
                total_train_loss = []
                total_train_loss_SupCon = []
                for index,data in enumerate(data_loader_train):
                    mean,base_num_fraction,std,time,avg_basequality,g_num,y=data
                    mean=mean.to(device)
                    std=std.to(device)
                    time=time.to(device)
                    avg_basequality=avg_basequality.to(device)
                    g_num=g_num.to(device)
                    base_num_fraction=base_num_fraction.to(device)
                    y=y.to(device)
                    input = torch.cat((mean,std,time,base_num_fraction,avg_basequality), dim=1)
                    outputs_basic_group=model_basic_group(input,base_num_fraction,g_num,avg_basequality,std,time) #mean,base_num_fraction,g_num,avg_basequality,std,time
                    loss_basic_group=criterion(outputs_basic_group, y) 
                    #loss_SupCon = criterion_SupCon(features=mid_mean, labels=y)
                    loss=loss_basic_group#+0.1*loss_SupCon
                    loss_value=loss.item()
                    total_train_loss.append(loss_value)
                    optimizer_basic_group.zero_grad()  
                    loss.backward()        
                    optimizer_basic_group.step() 
                scheduler_basic_group.step()
            
                model_basic_group.eval()
                total_valid_loss=[]
                total_valid_loss_SupCon=[]
                with torch.no_grad():
                    for index,data in enumerate(data_loader_valid):
                        mean,base_num_fraction,std,time,avg_basequality,g_num,y=data
                        mean=mean.to(device)
                        std=std.to(device)
                        time=time.to(device)
                        avg_basequality=avg_basequality.to(device)
                        g_num=g_num.to(device)
                        base_num_fraction=base_num_fraction.to(device)
                        y=y.to(device)
                        input = torch.cat((mean,std,time,base_num_fraction,avg_basequality), dim=1)
                        outputs_basic_group=model_basic_group(input,base_num_fraction,g_num,avg_basequality,std,time)
                        loss_basic_group=criterion(outputs_basic_group,y) 
                        loss_valid=loss_basic_group#+0.1*loss_SupCon
                        loss_value=loss_valid.item()
                        total_valid_loss.append(loss_value)
                if loss_valid.item() < best_loss - tol:
                    best_loss = loss_valid.item()
                    loss_not_improved_count = 0
                else:
                    loss_not_improved_count += 1

                if loss_not_improved_count >= n_iter_no_change:
                    print(f"Training stopped at epoch {epoch+1}",file=f)
                    torch.save(model_basic_group.state_dict(), f"/home/nipeng/chenzheng/mulan-methyl/my_model/mlp_mean_{motif}.pt")
                    break

                if (epoch + 1) % 5 == 0:
                    print(f'Epoch [{epoch + 1}/{num_epochs}],Train Loss: {np.mean(total_train_loss):.6f},Test Loss: {np.mean(total_valid_loss):.6f}',file=f)#,Train supCon Loss: {np.mean(total_train_loss_SupCon):.6f},Test supCon Loss: {np.mean(total_valid_loss_SupCon):.6f}')
                    #torch.save(model_basic_group.state_dict(), f"/home/nipeng/chenzheng/mulan-methyl/my_model/mlp_mean_{motif}.pt")
                #break
