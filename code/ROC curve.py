import os
import pandas as pd
import numpy as np
from sklearn.metrics import auc
import matplotlib.pyplot as plt

def Normalize(data):
    temp = data
    for i in range(data.shape[0]):
        total = sum(data[i])
        for j in range(data.shape[1]):
            temp[i][j] = data[i][j] / total
    return temp

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
    temparray = np.zeros((266,103))
    for i in lncRNA_index:
        temparray[m] = array[i]
        m += 1
    n = 0
    data = np.zeros((266,58))
    for j in protein_index:
        data[:,n] = temparray[:,j]
        n += 1
    #
    TPR = []
    FPR = []
    RECALL= []
    PRECISION = []
    #Arrange the thresholds from large to small
    datalist = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            datalist.append(data[i][j])
    dataindex = list(set(datalist))
    dataindex.sort(reverse=True)
    #
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
    return TPR,FPR,Auc,RECALL,PRECISION,aupr

#final scoring matrix
path = 'C:\\Users\\zyk\\Desktop\\PAPER\\data\\'
#content = pd.read_csv(os.path.join(path,'10_results and curves\\final scoring matrix.csv'),header=None,index_col=None)
content = pd.read_csv(os.path.join(path,'11_prediction besed on all database\\final scoring matrix.csv'),header=None,index_col=None)
score_matrix = np.array(content)

#species of test lncRNAs
#test_rna_content = pd.read_csv(os.path.join(path,'6_select lncRNA-protein interactions\\lncRNA species.csv'),header=None,index_col=None)
test_rna_content = pd.read_csv(os.path.join(path,'6_select lncRNA-protein interactions\\lncRNA species.csv'),header=None,index_col=None)
test_rna_data = np.array(test_rna_content)
test_rna_id = test_rna_data[:,2].tolist()
#species of test proteins
#test_protein_content = pd.read_csv(os.path.join(path,'6_select lncRNA-protein interactions\\protein species.csv'),header=None,index_col=None)
test_protein_content = pd.read_csv(os.path.join(path,'6_select lncRNA-protein interactions\\protein species.csv'),header=None,index_col=None)
test_protein_data = np.array(test_protein_content)
test_protein_id = test_protein_data[:,2].tolist()

#species of train lncRNAs
train_rna_content = pd.read_csv(os.path.join(path,'5_select based on common miRNA\\lncRNA species.csv'),header=None,index_col=None)
train_rna_data = np.array(train_rna_content)
train_rna_id = train_rna_data[:,2].tolist()
#species of train proteins
train_protein_content = pd.read_csv(os.path.join(path,'5_select based on common miRNA\\protein species.csv'),header=None,index_col=None)
train_protein_data = np.array(train_protein_content)
train_protein_id = train_protein_data[:,2].tolist()
#known lncRNA-protein interactions used for testing
#know_interactions_content = pd.read_csv(os.path.join(path,'8_known interactions correlation matrix\\known interactions.csv'),header=None,index_col=None)
know_interactions_content = pd.read_csv(os.path.join(path,'8_known interactions and correlation matrix\\known interactions.csv'),header=None,index_col=None)
know_interactions = np.array(know_interactions_content)

tpr,fpr,roc_auc,recall,precision,aupr = roc(score_matrix)
coordinate = []
coordinate.append(fpr)
coordinate.append(tpr)
coordinate.append(recall)
coordinate.append(precision)
coordinate = np.array(coordinate)
coordinate = coordinate.T
plt.figure(figsize=(12,4), dpi=300,facecolor=(1, 1, 1))
plt.tight_layout()
plt.subplots_adjust(wspace =0.2, hspace =0.2)
plt.subplot(121)
plt.plot(fpr,tpr,label=r'AUROC=%.3f'%roc_auc)
plt.xlim([-0.05, 1.05])
plt.xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
plt.ylim([-0.05, 1.05])
plt.yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
plt.xlabel('False Positive Rate',fontsize=11)
plt.ylabel('True Positive Rate',fontsize=11)
plt.title('ROC curve',fontsize=11)
plt.legend(loc="lower right",fontsize=7)
plt.subplot(122)
plt.plot(recall,precision,label=r'AUPR=%.3f'%aupr)
plt.xlim([-0.05, 1.05])
plt.xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
plt.ylim([-0.05, 1.05])
plt.yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
plt.xlabel('Recall',fontsize=11)
plt.ylabel('Precision',fontsize=11)
plt.title('PR curve',fontsize=11)
plt.legend(loc="upper right",fontsize=7)
#plt.show()
plt.savefig(os.path.join(path,'10_results and curves\curves\\ROC curve.png'))
plt.savefig(os.path.join(path,'10_results and curves\curves\\ROC curve.svg'))
np.savetxt(os.path.join(path,'10_results and curves\curves\\coordinate.csv'),coordinate,delimiter=',',fmt='%s')