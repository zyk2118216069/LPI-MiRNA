import pandas as pd
import numpy as np

path = 'C:\\Users\\zyk\\Desktop\\PAPER\\data\\11_prediction besed on all database\\'
def delete_repeat_interactions(path):
    content = pd.read_csv(path,header=None,index_col=None)
    data = np.array(content)
    temp = []
    for i in range(data.shape[0]):
        for j in range(i, data.shape[0]):
            if i != j and data[:,3][i] == data[:,3][j] and data[:,7][i] == data[:,7][j]:
                temp.append(j)
    array = np.delete(data, temp, axis=0)
    savepath = path + 'delete repeat lncRNA-protein interactions.csv'
    np.savetxt(savepath, array, delimiter=',', fmt = '%s')

delete_repeat_interactions('all lncRNA-protein interactions.csv')