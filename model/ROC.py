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
from sklearn.preprocessing import label_binarize
import matplotlib.pyplot as plt
os.environ['CUDA_VISIBLE_DEVICES'] = "3"
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
class ComplexResidualMLP(nn.Module):
    def __init__(self, layer_sizes, residual_layers):
        super(ComplexResidualMLP, self).__init__()
        self.layers = nn.ModuleList()
        for i in range(len(layer_sizes) - 1):
            self.layers.append(nn.Linear(layer_sizes[i], layer_sizes[i+1]))
            self.layers.append(nn.ReLU())
        self.residual_layers = residual_layers

    def forward(self,mean,basic,sd,mean_nat,mean_wga,sd_nat,sd_wga):
        residuals = []
        x=torch.cat((mean_nat,mean_wga,sd_nat,sd_wga),dim=2)
        for i, layer in enumerate(self.layers):
            if i in self.residual_layers:
                residuals.append(x)
            x = layer(x)
            if i in self.residual_layers:
                x = torch.cat([x] + residuals, dim=1)
        return x
#ComplexResidualMLP(layer_sizes=[4, 32, 36, 64, 104, 3], residual_layers=[1, 3])
class Flatten(nn.Module):
    def forward(self, input):
        return input.view(input.size(0), -1)
class MLP_basic(nn.Module):
    def __init__(self, input_size=22, hidden_size=128, num_classes=3):
        super(MLP_basic, self).__init__()
        mid_channels = 32
        mid_len = 10
        kernel_size = 3
        self.layers1 = nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=kernel_size, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),Flatten(),nn.Linear(mid_channels*mid_len, num_classes))
        self.layers2 = nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=kernel_size, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),Flatten(),nn.Linear(mid_channels*mid_len, num_classes))
        self.layers3 = nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=kernel_size, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),Flatten(),nn.Linear(mid_channels*mid_len, num_classes))
        self.layers4 = nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=kernel_size, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),Flatten(),nn.Linear(mid_channels*mid_len, num_classes))
        self.layers5 = nn.Sequential(nn.Conv1d(1, mid_channels, kernel_size=kernel_size, stride=1),nn.BatchNorm1d(mid_channels),nn.ReLU(),Flatten(),nn.Linear(mid_channels*mid_len, num_classes))
        #self.weights = nn.Parameter(torch.Tensor([100,1,0.1,0.5,0.3]))
        #"""
        self.weights = torch.Tensor([0.643180383
        ,0.058239536
        ,0.030066228
        ,0.100720973
        ,0.006191121
        ,0.038674039
        ])
        #"""

    def forward(self,mean,base_num_fraction,g_num,avg_basequality,std,time):
        mean = self.layers1(mean)
        base_num_fraction = self.layers2(base_num_fraction)
        avg_basequality = self.layers3(avg_basequality)
        std = self.layers4(std)
        time = self.layers5(time)
        return self.weights[0]*mean+self.weights[1]*base_num_fraction+self.weights[2]*avg_basequality+self.weights[3]*std+self.weights[4]*time
class MLP(nn.Module):
    def __init__(self, input_size=22, hidden_size=128, num_classes=3):
        super(MLP, self).__init__()
        #self.layers1 = nn.Sequential(nn.Linear(input_size,hidden_size),nn.ReLU(),nn.Linear(hidden_size, num_classes))
        self.layers1 = nn.Sequential(nn.Linear(input_size,hidden_size),nn.ReLU(),nn.Linear(hidden_size, num_classes))
        #self.layers2 = nn.Sequential(nn.Linear(input_size,num_classes))
        #self.layers3 = nn.Sequential(nn.Linear(input_size,num_classes))
        self.weights = nn.Parameter(torch.Tensor([100,10]))
        #print(self.weights)
        #self.hidden1 = nn.Linear(input_size,hidden_size)
        #self.output = nn.Linear(hidden_size, num_classes)
        #self.relu = nn.ReLU()

    def forward(self,mean,base_num_fraction,g_num,avg_basequality,std,time):
        mean = self.layers1(mean)
        #base_num_fraction = self.layers2(base_num_fraction)
        #std = self.layers3(std)
        #x=self.relu(self.hidden1(mean))
        #x = self.output(x)
        #return mean
        return mean#self.weights[0]*mean+self.weights[1]*base_num_fraction#+self.weights[2]*std
    
class CNN(nn.Module):
    def __init__(self):
        super(CNN, self).__init__()
        self.conv = nn.Sequential(
            nn.Conv2d(2, 32, kernel_size=3, stride=1, padding=1),
            nn.ReLU(inplace=True),
            nn.MaxPool2d(2),
            nn.Conv2d(32, 64, kernel_size=3, stride=1, padding=1),
            nn.ReLU(),
            nn.MaxPool2d(2),
        )
        self.fc = nn.Sequential(
           nn.Linear(64 * 3 * 32, 128),
           nn.Sigmoid(),
           nn.Dropout(0.5),
           nn.Linear(128,21)
        )

    def forward(self, x):
        #print(x.shape)
        x = self.conv(x)
        #print(x.shape)
        x = x.view(x.size()[0], -1)
        x = self.fc(x)
        return x
  
class TransformerModel(nn.Module):
    def __init__(self, input_size=12,feature_size=5,hidden_dim=64,num_classes=21):
        super(TransformerModel, self).__init__()
        # 构建Transformer编码层，参数包括输入维度、注意力头数
        # 其中d_model要和模型输入维度相同
        #self.positional_encoding = PositionalEncoding(feature_size)
        #self.embedding_mean = nn.Linear(1, 256)
        self.fc1 = nn.Linear(1, 63)
        self.fc2 = nn.Linear(63, 63)
        self.relu = nn.ReLU()
        self.embedding = nn.Embedding(4, 128) 
        self.positional_encoding = nn.Embedding(input_size, hidden_dim)
        self.encoder_layer = nn.TransformerEncoderLayer(d_model=hidden_dim, nhead=8)        
        self.encoder = nn.TransformerEncoder(self.encoder_layer,num_layers=4)                

        self.fc = nn.Linear(64,num_classes)
    def forward(self,mean,basic,sd,mean_nat,mean_wga,sd_nat,sd_wga):
        #mean,basic=x
        #basic = self.embedding(basic)
        #x=torch.cat((mean,mean_nat,mean_wga,sd_nat,sd_wga),dim=2)
        x=mean
        #"""
        x=self.fc1(x)
        x = self.relu(x)
        x = self.fc2(x)
        x = self.relu(x)
        x=torch.cat((x,mean),dim=2)
        #x=torch.cat((nat,wga),dim=1)
        #mean = self.embedding_mean(mean)
        #sd = self.embedding_sd(sd)
        #"""
        x=x.transpose(1,0)
        seq_length, batch_size, hidden_dim = x.shape
        position = torch.arange(seq_length).unsqueeze(1).repeat(1, batch_size).to(x.device)
        positional_encoding = self.positional_encoding(position)
        x = x + positional_encoding
        x = self.encoder(x)   
        x = torch.mean(x, dim=0)
        x = self.fc(x)        # 输入线性层进行分类预测

        return x

class LSTMModel(nn.Module):
    """
        Parameters:
        - input_size: feature size
        - hidden_size: number of hidden units
        - output_size: number of output
        - num_layers: layers of LSTM to stack
    """
    def __init__(self, input_size=128, hidden_size=512, output_size=21, num_layers=2,bidirectional=True):
        super().__init__()
        self.seq_len=12
        self.is_bidirectional=1
        if bidirectional: self.is_bidirectional=2
        self.embedding = nn.Embedding(4, input_size)
        #self.dropout1 = nn.Dropout(p=0.5)
        self.fc_after_add = nn.Linear(input_size,input_size)   
        self.lstm = nn.LSTM(input_size, hidden_size, num_layers,bidirectional=bidirectional,dropout=0.5) # utilize the LSTM model in torch.nn 
        #self.dropout2 = nn.Dropout(p=0.5)
        self.forwardCalculation = nn.Linear(hidden_size*self.seq_len*self.is_bidirectional, output_size)
 
    def forward(self, mean,basic,sd):
        #_x =  _x.unsqueeze(2)
        basic = self.embedding(basic)
        _x = basic+sd
        _x = self.fc_after_add(_x)
        _x = _x+mean
        #_x=self.dropout1(_x)
        x, _ = self.lstm(_x)  # _x is input, size (batch,seq_len, input_size)
        #print(x.shape)
        b, s, h = x.shape  # x is output, size (batch,seq_len, hidden_size)
        x = x.view(b,s*h)
        x = self.forwardCalculation(x)
        #x=self.dropout2(x)
        #x = x.view(s, b, -1)
        return x
def rollmean(t):
    t_rolling = torch.zeros_like(t)

    # 遍历每一行
    for i in range(t.shape[0]):
        # 遍历每一行的每一个元素
        for j in range(t.shape[1]):
            # 根据元素的位置来决定窗口的大小
            start = max(0, j - 2)
            end = min(t.shape[1], j + 3)
            # 计算窗口的平均值
            t_rolling[i, j] = t[i, start:end].mean()

    return t_rolling
class origin_dataset_train_basic(Dataset):
    def __init__(self, transform,valid):
        super(origin_dataset_train_basic, self).__init__()
        if valid==False:     
            path_data = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_mean_12_code_filter_iso.pkl'
            path_data_base_num_fraction = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_base_num_fraction_12_code_filter_iso.pkl"
            path_data_std_diff = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_std_diff_12_code_filter_iso.pkl"
            path_data_time_diff = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_time_diff_12_code_filter_iso.pkl"
            path_test_t_test_pval = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_t_test_pval_12_code_filter_iso.pkl"
            path_test_u_test_pval = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_u_test_pval_12_code_filter_iso.pkl"
            path_test_g_num = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_g_num_12_code_filter_iso.pkl"
            path_test_g_value = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_g_value_12_code_filter_iso.pkl"
            path_test_avg_basequality = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_avg_basequality_12_code_filter_iso.pkl"
            path_y_basic_group = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_basic_group_12_code_filter_iso.pkl'
            path_y_offset = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_offset_12_code_filter_iso.pkl'

        else:
            path_data = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_mean_12_code_filter_iso.pkl'
            path_data_base_num_fraction = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_base_num_fraction_12_code_filter_iso.pkl"
            path_data_std_diff = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_std_diff_12_code_filter_iso.pkl"
            path_data_time_diff = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_time_diff_12_code_filter_iso.pkl"
            path_test_t_test_pval = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_t_test_pval_12_code_filter_iso.pkl"
            path_test_u_test_pval = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_u_test_pval_12_code_filter_iso.pkl"
            path_test_g_num = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_g_num_12_code_filter_iso.pkl"
            path_test_g_value = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_g_value_12_code_filter_iso.pkl"
            path_test_avg_basequality = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_avg_basequality_12_code_filter_iso.pkl"
            path_y_basic_group = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_basic_group_12_code_filter_iso.pkl'
            path_y_offset = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_offset_12_code_filter_iso.pkl'

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

class origin_dataset_train(Dataset):
    def __init__(self, transform,valid):
        super(origin_dataset_train, self).__init__()
        if valid==False:     
            path_data = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_mean_12_code_filter_iso_offset_train.pkl'
            path_data_base_num_fraction = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_base_num_fraction_12_code_filter_iso_offset_train.pkl"
            path_data_std_diff = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_std_diff_12_code_filter_iso_offset_train.pkl"
            path_data_time_diff = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_time_diff_12_code_filter_iso_offset_train.pkl"
            path_test_t_test_pval = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_t_test_pval_12_code_filter_iso_offset_train.pkl"
            path_test_u_test_pval = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_u_test_pval_12_code_filter_iso_offset_train.pkl"
            path_test_g_num = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_g_num_12_code_filter_iso_offset_train.pkl"
            path_test_g_value = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_g_value_12_code_filter_iso_offset_train.pkl"
            path_test_avg_basequality = "/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_avg_basequality_12_code_filter_iso_offset_train.pkl"
            path_y_basic_group = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_basic_group_12_code_filter_iso_offset_train.pkl'
            path_y_offset = '/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_offset_12_code_filter_iso_offset_train.pkl'
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
        #self.data = self.data.unsqueeze(1)
        #self.data_nat =  self.data_nat.unsqueeze(1)
        #self.data_wga =  self.data_wga.unsqueeze(1)
    def __len__(self):
        return self.len_data
 
    def __getitem__(self, index):
        #data = torch.FloatTensor(data)
        return  self.data[index],self.data_base_num_fraction[index],self.data_std_diff[index],self.data_time_diff[index],self.data_t_test_pval[index],self.data_avg_basequality[index],self.data_g_num[index],self.data_g_value[index], self.label_basic_group[index],self.label_offset[index]
        #return self.data_basic_group[index],self.label[index]

class origin_dataset_test(Dataset):
    def __init__(self, transform,motif):
        super(origin_dataset_test, self).__init__()
        species="TP"
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
        self.data_basic_group=subet_test_data_basic_group.values[:,1:23] #1:23
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
    batch_size=256
    tol=1e-4
    n_iter_no_change=10
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
    #"""
    dataset_train_basic=origin_dataset_train_basic(transform=None,valid=False)
    data_loader_train_basic = DataLoader(
        dataset_train_basic,
        batch_size=batch_size,
        num_workers=4,
        shuffle=True,
        pin_memory=True
    )
    #"""
    """
    model_basic_group = MLP_basic(input_size=12, hidden_size=64, num_classes=21).to(device) 
    #model_basic_group = ResNet18().to(device)#ComplexResidualMLP(layer_sizes=[4, 32, 36, 64, 104, 3], residual_layers=[1, 3]).to(device) #ResNet18()  
    model_offset = TransformerModel().to(device) 
    #model = LSTMModel().to(device)  
    criterion = nn.CrossEntropyLoss()
    optimizer_basic_group = torch.optim.Adam(model_basic_group.parameters()) 
    optimizer_offset = torch.optim.Adam(model_offset.parameters(),lr=0.001) 
    #schedulr steplr 2-1 epoch *0.1 
    #scheduler_basic_group = StepLR(optimizer_basic_group, step_size=1, gamma=1)
    scheduler_offset = StepLR(optimizer_offset, step_size=2, gamma=0.1)
    #scheduler = ReduceLROnPlateau(optimizer, mode='max', factor=args.lr_decay, patience=args.lr_patience, verbose=True)
    num_epochs = 5000  
    loss_each_epoch=[]   
    for epoch in range(num_epochs):
        #if epoch%10==0 : print(epoch)
        model_basic_group.train()
        total_train_loss = []
        for index,data in enumerate(data_loader_train):
            mean,base_num_fraction,std,time,t_test_pval,avg_basequality,g_num,g_value,y_basic_group,y_offset=data
            mean=mean.to(device)
            std=std.to(device)
            time=time.to(device)
            t_test_pval=t_test_pval.to(device)
            avg_basequality=avg_basequality.to(device)
            g_num=g_num.to(device)
            g_value=g_value.to(device)
            base_num_fraction=base_num_fraction.to(device)
            y_basic_group=y_basic_group.to(device)
            y_offset=y_offset.to(device)
            y=y_basic_group*7+y_offset
            #outputs_basic_group=model_basic_group(mean,basic,sd,nat,wga)
            outputs_basic_group=model_basic_group(mean,base_num_fraction,g_num,avg_basequality,std,time) #mean,base_num_fraction,g_num,avg_basequality,std,time
            #print(outputs_basic_group.shape,y.shape)
            #print(y)
            #print(torch.max(outputs.data, 1))
            #print(outputs.shape,y.shape)
            #loss_basic_group=criterion(outputs_basic_group, y_basic_group)  
            loss_basic_group=criterion(outputs_basic_group, y) 
            loss=loss_basic_group
            loss_value=loss.item()
            total_train_loss.append(loss_value)
            #print(loss.item())
            optimizer_basic_group.zero_grad()  
            loss.backward()        
            optimizer_basic_group.step() 
        #scheduler_basic_group.step()
     
        model_basic_group.eval()
        total_valid_loss=[]
        with torch.no_grad():
            for index,data in enumerate(data_loader_valid):
                mean,base_num_fraction,std,time,t_test_pval,avg_basequality,g_num,g_value,y_basic_group,y_offset=data
                mean=mean.to(device)
                std=std.to(device)
                time=time.to(device)
                t_test_pval=t_test_pval.to(device)
                avg_basequality=avg_basequality.to(device)
                g_num=g_num.to(device)
                g_value=g_value.to(device)
                base_num_fraction=base_num_fraction.to(device)
                #x=x.to(device)
                y_basic_group=y_basic_group.to(device)
                y_offset=y_offset.to(device)
                y=y_basic_group*7+y_offset
                outputs_basic_group=model_basic_group(mean,base_num_fraction,g_num,avg_basequality,std,time)
                #print(outputs)
                #print(y)
                #print(torch.max(outputs.data, 1))
                #loss_basic_group=criterion(outputs_basic_group, y_basic_group)  
                loss_basic_group=criterion(outputs_basic_group, y) 
                loss=loss_basic_group
                loss_value=loss.item()
                total_valid_loss.append(loss_value)
      
        if loss.item() < best_loss - tol:
            best_loss = loss.item()
            loss_not_improved_count = 0
        else:
            loss_not_improved_count += 1

        if loss_not_improved_count >= n_iter_no_change:
            print(f"Training stopped at epoch {epoch+1}")
            torch.save(model_basic_group.state_dict(), f"/home/nipeng/chenzheng/mulan-methyl/my_model/mlp_{epoch+1}_data_mean_12_code_filter_iso_offset_21_avg.pt")
            break

        if (epoch + 1) % 5 == 0:
            print(f'Epoch [{epoch + 1}/{num_epochs}],Train Loss: {np.mean(total_train_loss):.6f}')#,Test Loss: {np.mean(total_valid_loss):.6f})
            #torch.save(model_basic_group.state_dict(), f"/home/nipeng/chenzheng/mulan-methyl/my_model/mlp_{epoch+1}_data_mean_12_code_filter_iso_offset_21.pt")
        #break
        if (epoch + 1) % 50 == 0:
            #print(f'Epoch [{epoch + 1}/{num_epochs}],Train Loss: {np.mean(total_train_loss):.6f}')#,Test Loss: {np.mean(total_valid_loss):.6f})
            torch.save(model_basic_group.state_dict(), f"/home/nipeng/chenzheng/mulan-methyl/my_model/mlp_{epoch+1}_data_mean_12_code_filter_iso_offset_21_avg.pt")
        #break
    """
    #model = ResNet18().to(device)  
    model = MLP_basic(input_size=12, hidden_size=64, num_classes=21).to(device)
    #model = joblib.load('/home/nipeng/chenzheng/mulan-methyl/saved_model/MLPClassifier.pkl')
    model_basic = MLP_basic(input_size=12, hidden_size=64, num_classes=21).to(device)
    model.load_state_dict(torch.load("/home/nipeng/chenzheng/mulan-methyl/my_model/mlp_43_data_mean_12_code_filter_iso_offset_21_conv.pt")) #mlp_46_data_mean_12_code_filter_iso_offset_21_avg.pt
    model_basic.load_state_dict(torch.load("/home/nipeng/chenzheng/mulan-methyl/my_model/mlp_43_data_mean_12_code_filter_iso_offset_21_conv.pt")) #mlp_1_data_mean_12_code_filter_iso_offset_3_avg
    #for name, param in model_basic.state_dict().items():
    #    print(name, param)
    #    break
    TP_motif=["CGCG_-14","CCTCC_-14","CAGAAA_-14","CCCRAG_-14","CTACT_-14","GATC_-14","GGNCC_-14","RAACTC_-14"]
    #TP_motif=["CGCG_-14","CCTCC_-14","CTACT_-14","GGNCC_-14","RAACTC_-14"]
    NO_motif=["CAANNNNNNNCTGG_-14","CCAGNNNNNNNTTG_-14","CTCGAG_-14","GCGGCCGC_-14"]
    true_label_TP=[10,11,17,18,10,17,10,10]
    true_label_TP_offset=[3,4,3,4,3,3,3,3]
    #true_label_TP_offset=[3,4,3,3,3]
    true_label_TP_basic_group=[1,1,2,2,1,2,1,1]
    true_label_NO=[17,17,4,10]
    MEAN=0
    y_4mC_prob_4mC=[]
    y_4mC_prob_6mA=[]
    y_6mA_prob_4mC=[]
    y_6mA_prob_6mA=[]
    #fig, axs = plt.subplots(8, 2, figsize=(16,10))
    for i,motif_item in enumerate(TP_motif):
        y_score=[]
        y_basic_score=[]
        predict_result=[]
        predict_basic_result=[]
        dataset_test=origin_dataset_test(transform=None,motif=motif_item)
        data_loader_test = DataLoader(
            dataset_test,
            batch_size=32,
            num_workers=1,
            shuffle=False,
            pin_memory=True,
            drop_last=False,
        )
        model.eval()
        model_basic.eval()
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

                mean_b=mean_b.to(device)
                avg_basequality_b=avg_basequality_b.to(device)
                base_num_fraction_b=base_num_fraction_b.to(device)
                std_b=std_b.to(device)
                time_b=time_b.to(device)
                #outputs=model.predict_proba(mean)#,base_num_fraction,g_num,avg_basequality,std)
                mean=mean.to(device)
                outputs=model(mean,base_num_fraction,g_num,avg_basequality,std,time)
                #outputs_basic=model_basic(mean_b,base_num_fraction_b,g_num,avg_basequality_b,std_b,time_b)
                _, predicted = torch.max(outputs.data, 1)
                #_, predicted_basic = torch.max(outputs_basic.data, 1)
                predict_prob=outputs.data
                #predict_basic_prob=outputs_basic.data
                #print(outputs.shape)
                """
                for b in range(outputs.shape[0]):
                    #predicted=sort_and_return_index(outputs.data[b])
                    predicted=sort_and_return_index(outputs[b])
                    predicted_offset=predicted%7
                    predicted_offset_index=0
                    for j in range(21):
                        #print(basic[b])
                        if basic[b][17-predicted_offset[j]]==1 or basic[b][17-predicted_offset[j]]==0: #A or C #8-predicted_offset[j] self.data_basic_group=subet_test_data_basic_group.values[:,6:18]
                            predicted_offset_index=j
                            break
                    final=predicted_basic[b]*7+predicted_offset[predicted_offset_index]
                    final=final.cpu().item()
                    #print(type(final))
                    predict_result.append(final)
                """
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
        #print(type(predict_result))
        lw=2  
        fpr = dict()
        tpr = dict()
        roc_auc = dict()
        y_test=np.array([true_label_TP[i]]*len(dataset_test))
        y_test = label_binarize(y_test, classes=range(21))
        y_score=np.array(y_score)
        #y_basic_score=np.array(y_basic_score)
        #predict_basic_result=np.array(predict_basic_result)
        predict_result=np.array(predict_result)
        """
        fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel())
        roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])
        plt.plot(fpr["micro"], tpr["micro"],
                label=motif_item+'(area = {0:0.2f})'.format(roc_auc["micro"]),
                 linestyle=':', linewidth=4)
        """
        precision, recall, _ = precision_recall_curve(y_test.ravel(), y_score.ravel())
        auc_score = auc(recall, precision)
        plt.plot(recall, precision, marker='.', label=motif_item+'(area = {0:0.2f})'.format(auc_score),
                 linestyle=':', linewidth=4)
        
        x=np.linspace(0,10,len(dataset_test))
        #"""


        result = pd.Series(predict_result).value_counts()
        #counts = np.bincount(predict_result)
        #most_num = np.argmax(counts)
        most_num=result.idxmax()
        MEAN+=result[true_label_TP[i]]/len(dataset_test)
        print(most_num,true_label_TP[i],":",result[true_label_TP[i]]/len(dataset_test))

    plt.plot([0, 1], [0, 1], 'k--', lw=lw)
    #plt.xlabel('False Positive Rate')
    #plt.ylabel('True Positive Rate')
    plt.legend(loc="lower right")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('PR-curve') #
    plt.savefig("/home/nipeng/chenzheng/mulan-methyl/TP_PR_conv_43.png")
