import os
import pandas as pd
import numpy as np
from sklearn.metrics import auc
import matplotlib.pyplot as plt

def roc(array):
    #Select test matrix based on known lncRNA-protein associations
    lncRNA_index = []
    for i in range(len(train_rna_id)):
        if train_rna_id[i] in test_rna_id:
            lncRNA_index.append(i) 
    protein_index = []
    for j in range(len(train_protein_id)):
        if train_protein_id[j] in test_protein_id:
            protein_index.append(j)
    m = 0
    temparray = np.zeros((len(test_rna_id),len(train_protein_id)))
    for i in lncRNA_index:
        temparray[m] = array[i]
        m += 1
    n = 0
    data = np.zeros((len(test_rna_id),len(test_protein_id)))
    for j in protein_index:
        data[:,n] = temparray[:,j]
        n += 1
    #Arrange the thresholds from large to small
    TPR = []
    FPR = []
    RECALL= []
    PRECISION = []
    datalist = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            datalist.append(data[i][j])
    dataindex = list(set(datalist))
    dataindex.sort(reverse=True)
    for threshold in dataindex:
        tp = 0
        fn = 0
        fp = 0
        tn = 0
        for i in range(len(test_rna_id)):
            for j in range(len(test_protein_id)):      
                if data[i][j] >= threshold:
                    if test_protein_id[j] in know_interactions[i]:
                        tp += 1
                    else:
                        fp += 1
                else:
                    if test_protein_id[j] in know_interactions[i]:
                        fn += 1
                    else:
                        tn += 1
        tpr = tp / (tp+fn)
        fpr = fp / (tn+fp)
        recall = tp/(tp+fn)
        precision = tp/(tp+fp)
        TPR.append(tpr)
        FPR.append(fpr)
        RECALL.append(recall)
        PRECISION.append(precision)
    Auc = auc(FPR,TPR)
    aupr = auc(RECALL,PRECISION)
    return TPR,FPR,RECALL,PRECISION,Auc,aupr

path = 'C:\\Users\\zyk\\Desktop\\PAPER\\data\\10_results and curves\\results\\'
datapath = path + 'net-net.csv'
content = pd.read_csv(datapath,header=None,index_col=None)
data = np.array(content)

datapath1 = path + 'net-seq.csv'
content1 = pd.read_csv(datapath1,header=None,index_col=None)
data1 = np.array(content1)

datapath2 = path + 'net-feat.csv'
content2 = pd.read_csv(datapath2,header=None,index_col=None)
data2 = np.array(content2)

datapath3 = path + 'seq-net.csv'
content3 = pd.read_csv(datapath3,header=None,index_col=None)
data3 = np.array(content3)

datapath4 = path + 'seq-seq.csv'
content4 = pd.read_csv(datapath4,header=None,index_col=None)
data4 = np.array(content4)

datapath5 = path + 'seq-feat.csv'
content5 = pd.read_csv(datapath5,header=None,index_col=None)
data5 = np.array(content5)

datapath6 = path + 'feat-net.csv'
content6 = pd.read_csv(datapath6,header=None,index_col=None)
data6 = np.array(content6)

datapath7 = path + 'feat-seq.csv'
content7 = pd.read_csv(datapath7,header=None,index_col=None)
data7 = np.array(content7)

datapath8 = path + 'feat-feat.csv'
content8 = pd.read_csv(datapath8,header=None,index_col=None)
data8 = np.array(content8)

path1 = 'C:\\Users\\zyk\\Desktop\\PAPER\\data\\'
test_rna_content = pd.read_csv(os.path.join(path1,'6_select lncRNA-protein interactions\\lncRNA species.csv'),header=None,index_col=None)
test_rna_data = np.array(test_rna_content)
test_rna_id = test_rna_data[:,2].tolist()
#species of test proteins
test_protein_content = pd.read_csv(os.path.join(path1,'6_select lncRNA-protein interactions\\protein species.csv'),header=None,index_col=None)
test_protein_data = np.array(test_protein_content)
test_protein_id = test_protein_data[:,2].tolist()
#species of train lncRNAs
train_rna_content = pd.read_csv(os.path.join(path1,'5_select based on common miRNA\\lncRNA species.csv'),header=None,index_col=None)
train_rna_data = np.array(train_rna_content)
train_rna_id = train_rna_data[:,2].tolist()
#species of train proteins
train_protein_content = pd.read_csv(os.path.join(path1,'5_select based on common miRNA\\protein species.csv'),header=None,index_col=None)
train_protein_data = np.array(train_protein_content)
train_protein_id = train_protein_data[:,2].tolist()
#known lncRNA-protein interactions used for testing
know_interactions_content = pd.read_csv(os.path.join(path1,'8_known interactions correlation matrix\\known interactions.csv'),header=None,index_col=None)
know_interactions = np.array(know_interactions_content)

coordinate = []
tpr1,fpr1,recall1,precision1,AUC1,PRUC1 = roc(data)
coordinate.append(fpr1)
coordinate.append(tpr1)
coordinate.append(recall1)
coordinate.append(precision1)
print('1')
tpr2,fpr2,recall2,precision2,AUC2,PRUC2 = roc(data1)
coordinate.append(fpr2)
coordinate.append(tpr2)
coordinate.append(recall2)
coordinate.append(precision2)
print('2')
tpr3,fpr3,recall3,precision3,AUC3,PRUC3 = roc(data2)
coordinate.append(fpr3)
coordinate.append(tpr3)
coordinate.append(recall3)
coordinate.append(precision3)
print('3')
tpr4,fpr4,recall4,precision4,AUC4,PRUC4 = roc(data3)
coordinate.append(fpr4)
coordinate.append(tpr4)
coordinate.append(recall4)
coordinate.append(precision4)
print('4')
tpr5,fpr5,recall5,precision5,AUC5,PRUC5 = roc(data4)
coordinate.append(fpr5)
coordinate.append(tpr5)
coordinate.append(recall5)
coordinate.append(precision5)
print('5')
tpr6,fpr6,recall6,precision6,AUC6,PRUC6 = roc(data5)
coordinate.append(fpr6)
coordinate.append(tpr6)
coordinate.append(recall6)
coordinate.append(precision6)
print('6')
tpr7,fpr7,recall7,precision7,AUC7,PRUC7 = roc(data6)
coordinate.append(fpr7)
coordinate.append(tpr7)
coordinate.append(recall7)
coordinate.append(precision7)
print('7')
tpr8,fpr8,recall8,precision8,AUC8,PRUC8 = roc(data7)
coordinate.append(fpr8)
coordinate.append(tpr8)
coordinate.append(recall8)
coordinate.append(precision8)
print('8')
tpr9,fpr9,recall9,precision9,AUC9,PRUC9 = roc(data8)
coordinate.append(fpr9)
coordinate.append(tpr9)
coordinate.append(recall9)
coordinate.append(precision9)
print('9')
all_coordinate = np.array(coordinate)
all_coordinate = all_coordinate.T
np.savetxt(r'C:\Users\zyk\Desktop\PAPER\data\10_results and curves\curves\all_roc_coordinate.csv',all_coordinate,delimiter=',',fmt='%s')

plt.figure(figsize=(12,4), dpi=300,facecolor=(1, 1, 1))
plt.tight_layout()
plt.subplots_adjust(wspace =0.2, hspace =0.2)

plt.subplot(121)
plt.plot(fpr1,tpr1,label=r'Net + Net AUROC=%.3f'%AUC1)
plt.plot(fpr2,tpr2,label=r'Net + Seq AUROC=%.3f'%AUC2)
plt.plot(fpr3,tpr3,label=r'Net + Feat AUROC=%.3f'%AUC3)
plt.plot(fpr4,tpr4,label=r'Seq + Net AUROC=%.3f'%AUC4)
plt.plot(fpr5,tpr5,label=r'Seq + Seq AUROC=%.3f'%AUC5)
plt.plot(fpr6,tpr6,label=r'Seq + Feat AUROC=%.3f'%AUC6)
plt.plot(fpr7,tpr7,label=r'Feat + Net AUROC=%.3f'%AUC7)
plt.plot(fpr8,tpr8,label=r'Feat + Seq AUROCC=%.3f'%AUC8)
plt.plot(fpr9,tpr9,label=r'Feat + Feat AUROC=%.3f'%AUC9)
plt.xlim([-0.05, 1.05])
plt.xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
plt.ylim([-0.05, 1.05])
plt.yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
plt.xlabel('False Positive Rate',fontsize=11)
plt.ylabel('True Positive Rate',fontsize=11)
plt.title('ROC curve',fontsize=11)
plt.legend(loc="lower right",fontsize=7)

plt.subplot(122)
plt.plot(recall1,precision1,label=r'Net + Net AUPR=%.3f'%PRUC1)
plt.plot(recall2,precision2,label=r'Net + Seq AUPR=%.3f'%PRUC2)
plt.plot(recall3,precision3,label=r'Net + Feat AUPR=%.3f'%PRUC3)
plt.plot(recall4,precision4,label=r'Seq + Net AUPR=%.3f'%PRUC4)
plt.plot(recall5,precision5,label=r'Seq + Seq AUPR=%.3f'%PRUC5)
plt.plot(recall6,precision6,label=r'Seq + Fea AUPR=%.3f'%PRUC6)
plt.plot(recall7,precision7,label=r'Feat + Net AUPR=%.3f'%PRUC7)
plt.plot(recall8,precision8,label=r'Feat + Seq AUPR=%.3f'%PRUC8)
plt.plot(recall9,precision9,label=r'Feat + Feat AUPR=%.3f'%PRUC9)
plt.xlim([-0.05, 1.05])
plt.xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
plt.ylim([-0.05, 1.05])
plt.yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
plt.xlabel('Recall',fontsize=11)
plt.ylabel('Precision',fontsize=11)
plt.title('PR curve',fontsize=11)
plt.legend(loc="upper right",fontsize=7)

plt.savefig(r'C:\Users\zyk\Desktop\PAPER\data\10_results and curves\curves\all_roc.svg',bbox_inches='tight')
plt.savefig(r'C:\Users\zyk\Desktop\PAPER\data\10_results and curves\curves\all_roc.npg',bbox_inches='tight')