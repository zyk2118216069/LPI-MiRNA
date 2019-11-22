import csv
import os
import pandas as pd
import numpy as np

datapath = 'C:\\Users\\zyk\\Desktop\\PAPER\\data\\2_delete repeat interactions\\'
savepath = 'C:\\Users\\zyk\\Desktop\\PAPER\\data\\3_select based on common miRNA\\'
procontent = pd.read_csv(os.path.join(datapath,'protein-miRNA interactions.csv'),header=None,index_col=None)
prodata = np.array(procontent)
pro_mirna = prodata[:,3]
pro_name = prodata[:,7]

rnacontent = pd.read_csv(os.path.join(datapath,'lncRNA-miRNA interactions.csv'),header=None,index_col=None)
rnadata = np.array(rnacontent)
lncrna_mirna = rnadata[:,3]
lncrna_name = rnadata[:,7]

miRNA = [com for com in pro_mirna if com in lncrna_mirna]
miRNA = list(set(miRNA))
miRNA.sort()
print('common miRNA: %d'%(len(miRNA)))

allprotein = []
for i in range(len(pro_mirna)):
    if pro_mirna[i] in miRNA:
        allprotein.append(pro_name[i])
print('protein-miRNA interactions: %d\nprotein species：%d'%(len(allprotein),len(set(allprotein))))

alllncRNA = []
for i in range(len(lncrna_mirna)):
    if lncrna_mirna[i] in miRNA:
        alllncRNA.append(lncrna_name[i])
print('lncRNA-miRNA interactions: %d\nlncRNA species：%d'%(len(alllncRNA),len(set(alllncRNA))))

protein_interactions = open(os.path.join(savepath,'protein-miRNA interactions.csv'),'w',newline='')
protein_interactions_csv = csv.writer(protein_interactions)
for i in range(len(pro_mirna)):
    if pro_mirna[i] in miRNA: 
        protein_interactions_csv.writerow(prodata[i])
protein_interactions.close()

lncRNA_interactions = open(os.path.join(savepath,'lncRNA-miRNA interactions.csv'),'w',newline='')
lncRNA_interactions_csv = csv.writer(lncRNA_interactions)
for i in range(len(lncrna_mirna)):
    if lncrna_mirna[i] in miRNA:
        lncRNA_interactions_csv.writerow(rnadata[i])
lncRNA_interactions.close()


protein_species = open(os.path.join(savepath,'protein species.csv'),'w',newline='')
protein_species_csv = csv.writer(protein_species)
proteinorder = list(set(allprotein)) 
proteinorder.sort()
for i in proteinorder:
    for j in range(len(pro_name)):
        if i == pro_name[j]:
            protein_species_csv.writerow(prodata[j][5:9])
            break
protein_species.close()

lncRNA_species = open(os.path.join(savepath,'lncRNA species.csv'),'w',newline='')
lncRNA_species_csv = csv.writer(lncRNA_species)
lncRNAorder = list(set(alllncRNA))
lncRNAorder.sort()
for i in lncRNAorder:
    for j in range(len(lncrna_name)):
        if i == lncrna_name[j]:
            lncRNA_species_csv.writerow(rnadata[j][5:9])
            break
lncRNA_species.close()