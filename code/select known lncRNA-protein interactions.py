import pandas as pd
import numpy as np
import csv
import os

path = 'C:\\Users\\zyk\\Desktop\\PAPER\\data\\11_prediction besed on all database\\'
proteincontent = pd.read_csv(os.path.join(path,'protein species.csv'),header=None,index_col=None)
proteindata = np.array(proteincontent)
proteinid = proteindata[:,2]

lncRNAcontent = pd.read_csv(os.path.join(path,'lncRNA species.csv'),header=None,index_col=None)
lncRNAdata = np.array(lncRNAcontent)
lncRNAid = lncRNAdata[:,2]

knowpath = r'C:\Users\zyk\Desktop\PAPER\data\11_prediction besed on all database\all lncRNA-protein interactions.csv'
knowcontent = pd.read_csv(knowpath,header=None,index_col=None)
knowdata = np.array(knowcontent)
knowrna = knowdata[:,3]
knowpro = knowdata[:,7]

rnatemp = []
protemp = []
for i in range(len(knowrna)):
    if knowrna[i] in lncRNAid and knowpro[i] in proteinid:
        rnatemp.append(knowrna[i])
        protemp.append(knowpro[i])

selcetrna = list(set(rnatemp))
selcetpro = list(set(protemp))

proteinkind = open(os.path.join(path,'final protein.csv'),'w',newline='')
proteinkind_csv = csv.writer(proteinkind)
h = 0
for i in range(len(proteinid)):
    if proteinid[i] in selcetpro:
        h += 1
        proteinkind_csv.writerow(proteindata[i])
proteinkind.close()
print('protein species: %d'%(h))

lncRNAkind = open(os.path.join(path,'final lncRNA.csv'),'w',newline='')
lncRNAkind_csv = csv.writer(lncRNAkind)
k = 0
for i in range(len(lncRNAid)):
    if lncRNAid[i] in selcetrna:
        k += 1
        lncRNAkind_csv.writerow(lncRNAdata[i])
lncRNAkind.close()
print('lncRNA species: %d'%(k))

known = open(os.path.join(path,'selected lncRNA-protein.csv'),'w',newline='')
known_csv = csv.writer(known)
m = 0
for i in range(len(knowrna)):
    if knowrna[i] in selcetrna and knowpro[i] in selcetpro:
        m += 1
        known_csv.writerow(knowdata[i])
known.close()
print('selected lncRNA-protein interactions: %d'%(m))