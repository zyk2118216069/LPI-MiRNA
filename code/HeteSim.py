import pandas as pd
import numpy as np
from sklearn.metrics import auc
#import matplotlib.pyplot as plt

def HeteSim(lncRNA_protein, protein_protein):
    #路径为LPP
    result = np.zeros((lncRNA_protein.shape[0],protein_protein.shape[1]))
    for i in range(lncRNA_protein.shape[0]):
        for j in range(protein_protein.shape[1]):
            result_i_j = np.dot(lncRNA_protein[i], protein_protein[:,j])
            result_i = np.linalg.norm(lncRNA_protein[i])
            result_j = np.linalg.norm(protein_protein[:,j]) 
            if result_i_j != 0:
                result[i][j] = result_i_j/(result_i*result_j)
    return result
    
def get_auc(score_matrix, real_label):
    score_list = []
    real_matrix = np.array(pd.read_csv(r'C:\Users\zyk\Desktop\PAPER\data\12_cmparison\interactions_matrix.csv',header=None,index_col=None,low_memory=False))
    for i in range(score_matrix.shape[0]):
        for j in range(score_matrix.shape[1]):
            score_list.append(score_matrix[i][j])
    score_list = list(set(score_list))
    score_list.sort(reverse = True)
    FPR = []
    TPR = []
    RECALL= []
    PRECISION = []
    for threshold in score_list:
        tp = 0
        fn = 0
        fp = 0
        tn = 0
        for i in range(score_matrix.shape[0]):
            for j in range(score_matrix.shape[1]):
                if score_matrix[i][j] >= threshold:
                    if real_matrix[i][j] == 1:
                        tp += 1
                    else:
                        fp += 1
                else:
                    if real_matrix[i][j] == 1:
                        fn += 1
                    else:
                        tn += 1
        fpr = fp / (tn+fp)
        tpr = tp / (tp+fn)
        recall = tp/(tp+fn)
        precision = tp/(tp+fp)
        FPR.append(fpr)
        TPR.append(tpr)
        RECALL.append(recall)
        PRECISION.append(precision)
    AUC = auc(FPR,TPR)
    aupr = auc(RECALL,PRECISION)
    return TPR,FPR,AUC,RECALL,PRECISION,aupr

if __name__ == '__main__':
    CV = 5
    lncRNA_protein = np.array(pd.read_csv(r'C:\Users\zyk\Desktop\PAPER\data\12_cmparison\interactions_matrix.csv',header=None,index_col=None,low_memory=False))
    protein_protein = np.array(pd.read_csv(r'C:\Users\zyk\Desktop\PAPER\data\12_cmparison\protein_interactions_matrix.csv',header=None,index_col=None,low_memory=False))
    final_matrix = np.zeros((lncRNA_protein.shape[0],lncRNA_protein.shape[1]))
    interactions_num = int(sum(sum(lncRNA_protein)))
    interactions_list = [i for i in range(interactions_num)]
    np.random.seed(1)
    np.random.shuffle(interactions_list)
    interactions_index = np.zeros((interactions_num,2))
    n = 0
    for i in range(lncRNA_protein.shape[0]):
        for j in range(lncRNA_protein.shape[1]):
            if lncRNA_protein [i][j] == 1:
                interactions_index[n][0] = i
                interactions_index[n][1] = j
                n += 1
    for k in range(CV):
        print(k+1)
        train_interactions = np.array(pd.read_csv(r'C:\Users\zyk\Desktop\PAPER\data\12_cmparison\interactions_matrix.csv',header=None,index_col=None,low_memory=False))
        select_interactions = interactions_list[int(interactions_num/CV)*k:int(interactions_num/CV)*(k+1)]
        for h in select_interactions:
            train_interactions[int(interactions_index[h][0])][int(interactions_index[h][1])] = 0
        result = HeteSim(train_interactions, protein_protein)
        #进行5折交叉验证，每次得到连接由1改为0这部分的关联分数，5次后可以得到所有已知关联的分数
        for h in select_interactions:
            final_matrix[int(interactions_index[h][0])][int(interactions_index[h][1])] = result[int(interactions_index[h][0])][int(interactions_index[h][1])]
        #每次将不存在已知关联部分的关联分数相加，得到5次之和
        for i in range(lncRNA_protein.shape[0]):
            for j in range(lncRNA_protein.shape[1]):
                if lncRNA_protein[i][j] == 0:
                    final_matrix[i][j] += result[i][j]
    #将不存在已知关联部分的关联分数除以5，取平均，得到最终的关联矩阵
    for i in range(lncRNA_protein.shape[0]):
            for j in range(lncRNA_protein.shape[1]):
                if lncRNA_protein[i][j] == 0:
                    final_matrix[i][j] /= CV
    tpr,fpr,roc_auc,recall,precision,aupr = get_auc(result,lncRNA_protein)
    coordinate = []
    coordinate.append(fpr)
    coordinate.append(tpr)
    coordinate.append(recall)
    coordinate.append(precision)
    coordinate = np.array(coordinate)
    coordinate = coordinate.T
    np.savetxt(r'C:\Users\zyk\Desktop\Hete_coordinate.csv',coordinate,delimiter=',',fmt='%s')