import re
import os
import numpy as np
import pandas as pd
import multiprocessing

match    =  3
mismatch = -3
gap      = -2

filename = r'C:\Users\zyk\Desktop\sw\rnaseq.txt'
data = open(filename).read()
pattern = re.compile(r'>.+\n')
allsequence = re.split(pattern,data)
del allsequence[0]

def s_w(seqA, allseq, savepath, num):
    scorelist = [0]*(num+1)
    print('Comparing the %d sequence'%(num+1))
    cols = len(seqA)
    for seqB in allseq:
        rows = len(seqB)
        matrix = [[0 for row in range(rows+1)] for col in range(cols+1)]
        paths = [[0 for row in range(rows+1)] for col in range(cols+1)]
        max_score = 0
        finalscore = 0
        for i in range(cols):
            for j in range(rows):
                if seqA[i] == seqB[j]:
                    diag = matrix[i][j] + match
                else:
                    diag = matrix[i][j] + mismatch
                up    = matrix[i + 1][j] + gap
                left  = matrix[i][j + 1] + gap
                score = max(0,diag, up, left)
                matrix[i+1][j+1] = score
                if score > max_score:
                    max_score = score
                    start_pos = [i+1, j+1]
                if matrix[i+1][j+1] == diag and matrix[i+1][j+1] != 0:
                    paths[i+1][j+1] = 'diag'
                elif matrix[i+1][j+1] == up   and matrix[i+1][j+1] != 0:
                    paths[i+1][j+1] = 'up'
                elif matrix[i+1][j+1] == left and matrix[i+1][j+1] != 0:
                    paths[i+1][j+1] = 'left'
        i, j = start_pos
        start_path = paths[i][j]
        while start_path != 0:
            finalscore += matrix[i][j]
            if start_path == 'diag':
                i, j = i-1, j-1
            elif start_path == 'up':
                j = j-1
            else:
                i = i-1
            start_path = paths[i][j]
        scorelist.append(finalscore)
    np.savetxt(os.path.join(path,num+'.csv'), scorelist, delimiter=',', fmt='%f')

if __name__ == '__main__':
    path = 'C:\\Users\\zyk\\Desktop\\sw\\score\\'
    pool = multiprocessing.Pool(processes = multiprocessing.cpu_count())
    for i in range(len(allsequence)):
        savepath = path + str(i+1) + '.txt'
        sequence1 = allsequence[i]
        sequence2 = allsequence[i+1:]
        pool.apply_async(s_w, (sequence1,sequence2,savepath,i,))
    pool.close()
    pool.join()
    scorematrix = []
    for i in range(len(allsequence)):
        alignpath = path + str(i) + '.csv'
        alignlist = pd.read_csv(alignpath,header=None,index_col=None)
        alignlist = np.array(alignlist)
        alignlist = alignlist.T
        scorematrix.append(alignlist[0])
    finalmatrix = np.array(scorematrix)
    for j in range(finalmatrix.shape[1]):
        for i in range(finalmatrix.shape[0]):
            finalmatrix[i][j] = finalmatrix[j][i]
    np.savetxt(os.path.join(path,r'score matrix.csv'), finalmatrix, delimiter=',', fmt='%f')