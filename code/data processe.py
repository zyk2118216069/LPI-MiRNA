import pandas as pd
import numpy as np
import math
import csv
import os

path = 'C:\\Users\\zyk\\Desktop\\PAPER\\data\\11_prediction besed on all database\\'
simpath = 'C:\\Users\\zyk\\Desktop\\PAPER\\data\\11_prediction besed on all database\\'
#lncRNA species
RNAcontent = pd.read_csv(os.path.join(path,'lncRNA species.csv'),header=None,index_col=None)
RNAdata = np.array(RNAcontent)
RNAID = RNAdata[:,2].tolist()

#protein species
PROcontent = pd.read_csv(os.path.join(path,'protein species.csv'),header=None,index_col=None)
PROdata = np.array(PROcontent)
PROID = PROdata[:,2].tolist()

#miRNA-lncRNA interaction
rnapath = path + '.csv'
rnacontent = pd.read_csv(os.path.join(path,'lncRNA-miRNA interactions.csv'),header=None,index_col=None)
rnadata = np.array(rnacontent)
rnarna = rnadata[:,7]
rnami = rnadata[:,3]

#protein-lncRNA interaction
procontent = pd.read_csv(os.path.join(path,'protein-miRNA interactions.csv'),header=None,index_col=None)
prodata = np.array(procontent)
propro = prodata[:,7]
promi = prodata[:,3]

def sumlist(data):
    temp = 0
    for i in data:
        temp += i
    return temp

def initialize(data,k=1):
    temp = data
    num = 0
    for i in range(data.shape[0]):
        num = sumlist(data[i])
        for j in range(data.shape[1]):
            if i == j and k == 1:
                temp[i][j] = 1
            else:
                if data[i][j] != 0:
                    temp[i][j] = data[i][j]/num
    return temp

def mirna_importance():
    #caluate the importance of miRNA in different network
    mirna = list(set(rnami))
    mirna.sort()
    rnamirna = []
    for i in mirna:
        n = 0
        for j in rnami:
            if i == j:
                n += 1
        rnamirna.append(n)    
    promirna = []
    for i in mirna:
        n = 0
        for j in promi:
            if i == j:
                n += 1
        promirna.append(n)
    m = []
    for i in rnamirna:
        m.append(-math.log(i/len(rnami)))
    n = []
    for j in rnamirna:
        n.append(-math.log(j/len(promi)))
    h = []
    for x in range(len(mirna)):
        num = rnamirna[x] + promirna[x]
        h.append(-math.log(num/(len(promi)+len(rnami))))
    savepath = r'C:\Users\zyk\Desktop\PAPER\data\11_prediction besed on all database\improtance of miRNAs.csv'
    miRNA = open(savepath,'w',newline='')
    miRNA_csv = csv.writer(miRNA)
    for i in range(len(mirna)):
        temp = []
        temp.append(mirna[i])
        temp.append(m[i])
        temp.append(n[i])
        temp.append(h[i])
        miRNA_csv.writerow(temp)
    miRNA.close()

def know_interaction():
    #process known lncRNA-protein interactions for calculation
    RNAcontent1 = pd.read_csv(os.path.join(path,'final lncRNA.csv'),header=None,index_col=None)
    RNAdata1 = np.array(RNAcontent1)
    RNAID1 = RNAdata1[:,2].tolist()

    knowcontent = pd.read_csv(os.path.join(path,'selected lncRNA-protein.csv'),header=None,index_col=None)
    knowdata = np.array(knowcontent)
    knowrna = knowdata[:,3]
    knowpro = knowdata[:,7]
    
    savepath = r'C:\Users\zyk\Desktop\PAPER\data\11_prediction besed on all database\know.csv'
    save = open(savepath,'w',newline='')
    save_csv = csv.writer(save)
    for i in range(len(RNAID1)):
        temp = []
        for j in range(len(knowrna)):
            if RNAID1[i] == knowrna[j]:
                temp.append(knowpro[j])
        save_csv.writerow(temp)
    save.close()
  
def whole_correlation_matrix():
    #get the whole network correlation matrix
    miRNApath = r'C:\Users\zyk\Desktop\PAPER\data\11_prediction besed on all database\improtance of miRNAs.csv'
    miRNAcontent = pd.read_csv(miRNApath,header=None,index_col=None)
    miRNAdata = np.array(miRNAcontent)
    miRNAID = miRNAdata[:,0].tolist()
    miRNA1 = miRNAdata[:,1].tolist()
    miRNA2 = miRNAdata[:,2].tolist()
    miRNA3 = miRNAdata[:,3].tolist()
    
    pro_pro = np.zeros((len(PROID),len(PROID)))
    for i in range(len(PROID)):
        for j in range(len(propro)):
            if PROID[i] == propro[j]:
                for h in range(len(promi)):
                    if promi[j] == promi[h] and j != h:
                        m = miRNAID.index(promi[h])
                        n = PROID.index(propro[h])
                        pro_pro[i][n] += miRNA1[m]                    
    rna_rna = np.zeros((len(RNAID),len(RNAID)))
    for i in range(len(RNAID)):
        for j in range(len(rnarna)):
            if RNAID[i] == rnarna[j]:
                for h in range(len(rnami)):
                    if rnami[j] == rnami[h] and j != h:
                        m = miRNAID.index(rnami[h])
                        n = RNAID.index(rnarna[h])
                        rna_rna[i][n] += miRNA2[m]             
    initial = np.zeros((len(RNAID),len(PROID)))
    for i in range(len(RNAID)):
        for j in range(len(rnarna)):
            if RNAID[i] == rnarna[j]:
                for h in range(len(promi)):
                    if rnami[j] == promi[h]:
                        n = PROID.index(propro[h])
                        m = miRNAID.index(promi[h])
                        initial[i][n] += miRNA3[m]

    np.savetxt(os.path.join(simpath, "protein-protein similarity.csv"), pro_pro, delimiter=',',fmt='%f')
    np.savetxt(os.path.join(simpath, "lncRNA-lncRNA similarity.csv"), rna_rna, delimiter=',',fmt='%f')
    np.savetxt(os.path.join(simpath, "whole_matrix.csv"), initial, delimiter=',',fmt='%f')

def final_matrix():
    #get the final scoring matrix
    datacontent = pd.read_csv(os.path.join(simpath, "whole_matrix.csv"),header=None,index_col=None)
    data = np.array(datacontent)

    rnarnacontent = pd.read_csv(os.path.join(simpath, "lncRNA-lncRNA similarity.csv"),header=None,index_col=None)
    rnarnadata = initialize(np.array(rnarnacontent))

    proprocontent = pd.read_csv(os.path.join(simpath, "protein-protein similarity.csv"),header=None,index_col=None)
    proprodata = initialize(np.array(proprocontent))

    temparray = np.matmul(rnarnadata,data)
    final_matrix = np.matmul(temparray,proprodata)
    np.savetxt(os.path.join(path, "final scoring.csv"), final_matrix, delimiter=',',fmt='%f')