import numpy as np
import pandas as pd


lenapth = r'C:\Users\zyk\Desktop\PAPER\data\5_select based on common miRNA\lncRNA species.csv'
content = pd.read_csv(lenapth,header=None,index_col=None,low_memory=False)
data = np.array(content)
lendata = data[:,5]

swpath = r'C:\Users\zyk\Desktop\PAPER\data\9_similarity\Sequence similarity\lncRNA score.csv'
swcontent = pd.read_csv(swpath,header=None,index_col=None,low_memory=False)
swscore = np.array(swcontent)

similarity = np.zeros((331,331))
for i in range(331):
    for j in range(i, 331):
        if i == j:
            similarity[i][j] = similarity[j][i] = 0
        else:
            similarity[i][j] = similarity[j][i] = swscore[i][j]/(lendata[i]+lendata[j])
np.savetxt(r"C:\Users\zyk\Desktop\PAPER\data\9_similarity\Sequence similarity\lncRNA-lncRNA similarity.csv", similarity, delimiter=',',fmt='%f')