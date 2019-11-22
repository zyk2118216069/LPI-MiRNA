import numpy as np
import pandas as pd
import math
import os

path = 'C:\\Users\\zyk\\Desktop\\pseknc\\pseknccsv\\'
savepath = 'C:\\Users\\zyk\\Desktop\\pseknc\\psekncsim\\'
def compute(num):
    datapath = path + str(num) + '.csv'
    content = pd.read_csv(datapath,header=None,index_col=None)
    data = np.array(content)
    result = np.zeros((331,331))
    for i in range(331):
        for j in range(331):
            if i != j:
                temp = 0
                for k in range(len(data[i])):
                    temp += math.pow((data[i][k] - data[j][k]),2)
                if temp == 0:
                    result[i][j] = result[j][i] = 1
                else:
                    temp2 = 1/math.sqrt(temp)
                    result[i][j] = result[j][i] = math.pow(temp2,2)  
    np.savetxt(os.path.join(savepath,str(num)), result, delimiter=',', fmt='%f')

for i in range(1,101):
    compute(i)