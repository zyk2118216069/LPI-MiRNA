Predicting lncRNA-protein interactions with miRNAs as mediators in a heterogeneous network model
本方法选取来自RAID V2.0中经过实验验证的关联作为实验材料，下载地址为http://www.rna-society.org/raid2/download.html，其中lncRNA-miRNA关联为2862条，protein-miRNA关联为2521条，lncRNA-protein关联40668条，见文件1_all interactions
删除其中重复的关联，最后得到lncRNA-miRNA关联2847条，protein-miRNA关联为2521条，lncRNA-protein关联40650条，见文件2_delete repeat interactions，删除代码见delete repeated connection.py
之后选择有共同miRNA的lncRNA-miRNA关联和protein-miRNA关联，最终得到lncRNA-miRNA关联1367条，protein-miRNA关联为1158条，包括341种lncRNA，103种protein和361种共同的miRNA，见文件3_select based on common miRNA，选择代码见select based on common miRNA.py
由于需要lncRNA和protein的序列信息，所以删除其中无法获得序列的lncRNA和protein，有10个lncRNA无法获得序列，最终得到lncRNA-miRNA关联1356条，protein-miRNA关联为1158条，包括331种lncRNA和103种protein，见文件4_select based on sequence
由于删除部分关联，因此再次选择有共同miRNA的lncRNA-miRNA关联和protein-miRNA关联，最终得到lncRNA-miRNA关联1356条，protein-miRNA关联为1156条，包括331种lncRNA，103种protein和361种共同的miRNA，lncRNA和protein的序列和ID见文件5_select based on common miRNA
然后根据lncRNA和protein的种类来选择已知关联来测试最终的结果，最终得到lncRNA-protein关联1925条，包括266种lncRNA和58种protein，见文件6_select lncRNA-protein interactions，选择方法见代码select known lncRNA-protein interactions.py
分别计算miRNA在lncRNA-miRNA异构网络，protein-miRNA异构网络以及lncRNA-miRNA-protein异构网络中的网络贡献度，见文件7_improtance of miRNAs，具体代码见data processe.py
然后处理已知关联便于进行测试，并根据miRNA的网络贡献度来计算初步关联度矩阵，见文件8_known interactions and correlation matrix，具体代码见data processe.py
分别计算lncRNA和protein之间的相似度，本文中计算三种相似度：网络相似度，序列相似度（Smith-Waterman）以及统计特征相似度（PseKNC和PseAAC），见文件9_similarity，Smith-Waterman分数计算方式见代码Smith-Waterman.py，相似度计算见代码Smith-Waterman similarity.py，Pse相似度计算见代码Pse similarity.py，网络相似度计算见data processe.py
最终结果以及图片见文件10_results and curves，ROC代码见ROC curve.py，另外All_ROC curve.py用于画所有的ROC曲线，Single_ROC curve用于画单个相似度的ROC曲线
使用整个RAID V2.0数据库的关联进行实验，具体步骤同上，所有过程见文件11_prediction besed on all database
最后与HeteSim算法进行比较，protein之间的关联以及所有的过程见文件12_cmparison，HeteSim具体计算代码见HeteSim.py
