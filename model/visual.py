import umap
import joblib
import matplotlib.pyplot as plt
import numpy as np
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score
x_mean, x_std, x_min, x_max = 0, 0, 0, 0
def scatter_plot(cluster_X,y):
    plt.figure(figsize=(20, 20))
    plt.scatter(cluster_X[y=="CGCG_-9", 0], cluster_X[y=="CGCG_-9", 1],c='blue',label="3")
    plt.scatter(cluster_X[y=="CCTCC_-9", 0], cluster_X[y=="CCTCC_-9", 1],c='red',label="4")
    plt.scatter(cluster_X[y=="CAGAAA_-9", 0], cluster_X[y=="CAGAAA_-9", 1],c='blue',label="3")
    plt.scatter(cluster_X[y=="CCCRAG_-9", 0], cluster_X[y=="CCCRAG_-9", 1],c='red',label="4")
    plt.scatter(cluster_X[y=="CTACT_-9", 0], cluster_X[y=="CTACT_-9", 1],c='blue',label="3")
    plt.scatter(cluster_X[y=="GATC_-9", 0], cluster_X[y=="GATC_-9", 1],c='blue',label="3")
    plt.scatter(cluster_X[y=="GGNCC_-9", 0], cluster_X[y=="GGNCC_-9", 1],c='blue',label="3")
    plt.scatter(cluster_X[y=="RAACTC_-9", 0], cluster_X[y=="RAACTC_-9", 1],c='blue',label="3")
    plt.legend(loc='best', fontsize=32,markerscale=3)
    plt.show()
    plt.savefig("/home/nipeng/chenzheng/DNA_Graduation_project/offset.png")

def scatter_plot_train_data(cluster_X,y):
    plt.figure(figsize=(20, 20))
    plt.scatter(cluster_X[(y==0), 0], cluster_X[(y==0), 1],c='blue',label="5mC")
    plt.scatter(cluster_X[(y==1), 0], cluster_X[(y==1), 1],c='red',label="4mC")
    plt.scatter(cluster_X[(y==2), 0], cluster_X[(y==2), 1],c='yellow',label="6mA")
    #plt.scatter(x[:, 0], x[:, 1],c='green',label="X_5mC")
    """
    plt.scatter(X[(Y_offset==0), 0], X[(Y_offset==0), 1],c='cyan',label="X_offset_3")
    plt.scatter(X[(Y_offset==1), 0], X[(Y_offset==1), 1],c='orange',label="X_offset_2")
    plt.scatter(X[(Y_offset==2), 0], X[(Y_offset==2), 1],c='brown',label="X_offset_1")

    plt.scatter(X[(Y_offset==3), 0], X[(Y_offset==3), 1],c='olive',label="X_offset_0")
    plt.scatter(X[(Y_offset==4), 0], X[(Y_offset==4), 1],c='gray',label="X_offset_-1")
    plt.scatter(X[(Y_offset==5), 0], X[(Y_offset==5), 1],c='magenta',label="X_offset_-2")
    plt.scatter(X[(Y_offset==6), 0], X[(Y_offset==6), 1],c='white',label="X_offset_-3") 
    plt.scatter(cluster_X[(y==0), 0], cluster_X[(y==0), 1],c='blue',label="5mC")
    plt.scatter(cluster_X[(y==1), 0], cluster_X[(y==1), 1],c='red',label="4mC")
    plt.scatter(cluster_X[(y==2), 0], cluster_X[(y==2), 1],c='yellow',label="6mA")
    plt.scatter(X[(Y==0), 0],X[(Y==0), 1],c='black',label="X_5mC")
    plt.scatter(X[(Y==1), 0],X[(Y==1), 1],c='orange',label="X_4mC")
    plt.scatter(X[(Y==2), 0],X[(Y==2), 1],c='brown',label="X_6mA")

    plt.scatter(cluster_X[(y==0)&(y_offset==3), 0], cluster_X[(y==0)&(y_offset==3), 1],c='blue',label="5mC")
    plt.scatter(cluster_X[(y==1)&(y_offset==3), 0], cluster_X[(y==1)&(y_offset==3), 1],c='red',label="4mC")
    plt.scatter(cluster_X[(y==2)&(y_offset==3), 0], cluster_X[(y==2)&(y_offset==3), 1],c='yellow',label="6mA")

    plt.scatter(cluster_X[y==3, 0], cluster_X[y==3, 1],c='green',label="3")
    plt.scatter(cluster_X[y==4, 0], cluster_X[y==4, 1],c='black',label="4")
    plt.scatter(cluster_X[y==5, 0], cluster_X[y==5, 1],c='pink',label="5")
    plt.scatter(cluster_X[y==6, 0], cluster_X[y==6, 1],c='purple',label="6")
    plt.scatter(X[(Y==0)&(Y_offset==3), 0],X[(Y==0)&(Y_offset==3), 1],c='black',label="X_0")
    plt.scatter(X[(Y==1)&(Y_offset==3), 0], X[(Y==1)&(Y_offset==3), 1],c='orange',label="X_1")
    plt.scatter(X[(Y==2)&(Y_offset==3), 0],X[(Y==2)&(Y_offset==3), 1],c='brown',label="X_2")
    plt.scatter(X[(Y==0), 0],X[(Y==0), 1],c='black',label="X_0")
    plt.scatter(X[(Y==1), 0],X[(Y==1), 1],c='orange',label="X_1")
    plt.scatter(X[(Y==2), 0],X[(Y==2), 1],c='brown',label="X_2")
    """
    plt.legend(loc='best', fontsize=32,markerscale=3)
    plt.xticks([])
    plt.yticks([])
    plt.show()
    #plt.title("UMAP for 46 motifs")
    plt.savefig("/home/nipeng/chenzheng/DNA_Graduation_project/train_basic_avg_basequality.png")

def normlize_features(features,train=False):
    if train:
        global x_mean, x_std, x_min, x_max
        x_mean, x_std = features.mean(), features.std()
        features = (features - x_mean) / x_std 
        x_min, x_max = features.min(), features.max()
        print(x_mean, x_std, x_min, x_max)
    else:
        print(x_mean, x_std, x_min, x_max)
        features = (features - x_mean) / x_std
    features = (features - x_min) / (x_max - x_min)  
    return features

x=joblib.load("/home/nipeng/chenzheng/mulan-methyl/data_pkl/CTCGAG_mean.pkl") 
y=joblib.load("/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_basic_group_valid_1_290.pkl") 
y_offset=joblib.load("/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_offset_valid_1_290.pkl")  
#/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_mean_12_code_filter_iso_train.pkl
x_train=joblib.load("/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_avg_basequality_12_code_filter_iso_offset_train.pkl") 
y_train=joblib.load("/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_basic_group_12_code_filter_iso_offset_train.pkl")
y_train_offset=joblib.load("/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_offset_12_code_filter_iso_offset_train.pkl") 
y_train_21=y_train*7+y_train_offset
print(x_train.shape,y_train.shape)
""" 
x=joblib.load("/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_mean_valid.pkl") 
y=joblib.load("/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_basic_group_valid.pkl") 
y_offset=joblib.load("/home/nipeng/chenzheng/mulan-methyl/data_pkl/y_offset_valid.pkl")  

x_train=joblib.load("/home/nipeng/chenzheng/DNABERT_2/basic_emb_S.pkl") #/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_mean_12_origin_filter_iso_offset.pkl 
y_train=joblib.load("/home/nipeng/chenzheng/DNABERT_2/Y_basic_group_S.pkl")
y_train_offset=joblib.load("/home/nipeng/chenzheng/DNABERT_2/Y_offset_S.pkl") 
""" 
"""
indices = np.where(y_train == 2)
x_train=x_train[indices]
y_train=y_train[indices]
y_train_offset=y_train_offset[indices]

X=np.concatenate((x_train,x),axis=0)
X=normlize_features(X,train=True)
x_train=normlize_features(x_train)
x=normlize_features(x)

tsne = TSNE(n_components=2, random_state=0)
tsne_result = tsne.fit(X)
tsne_result = tsne.fit(x)
"""
#joblib.dump(tsne_result,"/home/nipeng/chenzheng/mulan-methyl/data_pkl/data_mean_22_tsne.pkl")
umap_model = umap.UMAP() 
#tsne = TSNE(n_components=2, random_state=0)
x_train=normlize_features(x_train,train=True)
umap_result = umap_model.fit_transform(x_train)
#x=normlize_features(x,train=False)
#umap_result2 = umap_model.transform(x)
"""
kmeans = KMeans(n_clusters=3, random_state=1).fit(umap_result)
ari = adjusted_rand_score(y_train, kmeans.labels_)
print(ari)
"""
#umap_result2 = umap_model.transform(x)

#umap_result = umap_model.transform(x_train)
#umap_result2 = umap_model.transform(x)

scatter_plot_train_data(umap_result,y_train)