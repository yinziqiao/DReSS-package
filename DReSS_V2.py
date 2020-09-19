# -*- coding: utf-8 -*-
"""
Created on Thu Jan 4 15:50:00 2020

@author: 殷子樵


import numpy as np             
import itertools               
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import random
import csv
import codecs
import xlwt
import xlrd
import time
import fitter

"""
    # regular graphy
    # generate a regular graph which has 20 nodes & each node has 3 neghbour nodes.
    RG = nx.random_graphs.random_regular_graph(3, 20)
    # the spectral layout
    pos = nx.spectral_layout(RG)
    # draw the regular graphy
    nx.draw(RG, pos, with_labels = False, node_size = 30)
    plt.show()
    
    # erdos renyi graph
    # generate a graph which has n=20 nodes, probablity p = 0.2.
    ER = nx.random_graphs.erdos_renyi_graph(20, 0.2)
    # the shell layout
    pos = nx.shell_layout(ER)
    nx.draw(ER, pos, with_labels = False, node_size = 30)
    plt.show()
    
    # WS network
    # generate a WS network which has 20 nodes,
    # each node has 4 neighbour nodes,
    # random reconnection probability was 0.3.
    WS = nx.random_graphs.watts_strogatz_graph(20, 4, 0.3)
    # circular layout
    pos = nx.circular_layout(WS)
    nx.draw(WS, pos, with_labels = False, node_size = 30)
    plt.show()
    
    # BA scale-free degree network
    # generalize BA network which has 20 nodes, m = 1
    BA = nx.random_graphs.barabasi_albert_graph(20, 1)
    # spring layout
    pos = nx.spring_layout(BA)
    nx.draw(BA, pos, with_labels = False, node_size = 30)
    plt.show()
    
    DF_re = pd.DataFrame(np.zeros([len(G.nodes()),len(G.nodes())]),index = G.nodes(),columns = G.nodes())
    for col_label,row_label in G.edges():
        DF_re.loc[col_label,row_label] = 1 
        DF_re.loc[row_label,col_label] = 1 
"""

def main(A = None, mode = (1,0), threshold = None, print_switch = 'on'):
    """Main Function for calculating DReSS
    There are four parameters
    A is adjacency matrix of origin network
    mode can be chosen from follows:
                     (1,0) --- perturbation that delete existing interactions
                     (0,1) --- perturbation that add activation interactions
                     (0,-1) --- perturbation that add inhibition interactions
    print_switch is a parameter choose to display current progress or not
    """
    start = time.process_time()
    if A == None:
        A = (
                (-1,0,-1,-1,0,0,0,0,0),
                (0,0,-1,-1,0,0,-1,1,0),
                (0,-1,0,0,0,-1,0,0,0),
                (0,-1,0,0,0,-1,0,0,0),
                (0,-1,0,0,-1,-1,0,0,1),
                (0,0,-1,-1,1,0,0,0,0),
                (0,0,0,0,0,-1,0,0,0),
                (0,0,0,0,0,1,0,0,0),
                (0,0,1,1,0,0,1,-1,-1)
            )
    if threshold == None:
        threshold = [0,-1,0,0,0,0,0,0,0]
    """Parameter setting part ends here"""
    A_1 = list(A)
    for i in range(len(A_1)):
        A_1[i] = list(A_1[i])    
    """unzip part ends here"""

    if mode[0] != 0:         #mode determine
        pos = mat_find_nonzero(A_1)    #In deleting mode, find all non-zero position in A
    if mode[0] == 0:         
        pos = mat_find_zero(A_1)       #In adding mode, find all non-zero position in A
    
    basin_ori = Basin_with_plot(A = A,threshold = threshold)   #Calculation of origin attractor basin
    basin_o = basin_ori[1]
        
    pos_haming_simi_basin = {}
    for i in range(len(pos)):
        temp_pos = pos[i]
        A_tem = A_replace(A = A_1,position = temp_pos,interaction = mode[1])
        basin_tem = Basin_with_plot(A = A_tem,threshold = threshold)
        basin_t = basin_tem[1]
        temp = DReSS(A = basin_o, B = basin_t)
        pos_haming_simi_basin[temp_pos] = temp
        A_1 = list(A)
        if print_switch == 'on':
            print((i+1)/len(range(len(pos))))
            
    elapsed = (time.process_time() - start)
    print("总计用时",elapsed)
    return pos_haming_simi_basin

def diag_main(A = None, mode = (1,0), threshold = None, print_switch = 'on'):
    """Main Function for calculating diagDReSS
    There are four parameters
    A is adjacency matrix of origin network
    mode can be chosen from follows:
                     (1,0) --- perturbation that delete existing interactions
                     (0,1) --- perturbation that add activation interactions
                     (0,-1) --- perturbation that add inhibition interactions
    print_switch is a parameter choose to display current progress or not
    """
    start = time.process_time()
    if A == None:
        A = (
                (-1,0,-1,-1,0,0,0,0,0),
                (0,0,-1,-1,0,0,-1,1,0),
                (0,-1,0,0,0,-1,0,0,0),
                (0,-1,0,0,0,-1,0,0,0),
                (0,-1,0,0,-1,-1,0,0,1),
                (0,0,-1,-1,1,0,0,0,0),
                (0,0,0,0,0,-1,0,0,0),
                (0,0,0,0,0,1,0,0,0),
                (0,0,1,1,0,0,1,-1,-1)
            )
    if threshold == None:
        threshold = [0,-1,0,0,0,0,0,0,0]
    """Parameter setting part ends here"""
    A_1 = list(A)
    for i in range(len(A_1)):
        A_1[i] = list(A_1[i])    
    """unzip part ends here"""

    if mode[0] != 0:         #mode determine
        pos = mat_find_nonzero(A_1)    #In deleting mode, find all non-zero position in A
    if mode[0] == 0:         
        pos = mat_find_zero(A_1)       #In adding mode, find all non-zero position in A
    
    basin_ori = Basin_with_plot(A = A,threshold = threshold)   #Calculation of origin attractor basin
    basin_o = basin_ori[1]
    basin_o_diag = diag_array(basin_o)
        
    pos_haming_simi_basin = {}
    for i in range(len(pos)):
        temp_pos = pos[i]
        A_tem = A_replace(A = A_1,position = temp_pos,interaction = mode[1])
        basin_tem = Basin_with_plot(A = A_tem,threshold = threshold)
        basin_t = basin_tem[1]
        basin_t_diag = diag_array(basin_t)
        temp = haming_diag(A = basin_o_diag, B = basin_t_diag)
        pos_haming_simi_basin[temp_pos] = temp
        A_1 = list(A)
        if print_switch == 'on':
            print('当前进度为：', (i+1)/len(range(len(pos))))
            
    elapsed = (time.process_time() - start)
    print("总计用时",elapsed)
    return pos_haming_simi_basin

def Bonferroni(A = None, B = None):
    if A == None:
        A = {(0,1):0.9,(0,4):0.8}
    if B == None:
        B = [0.91,0.81,0.71,0.61]
        
    m = len(A)
    n = len(B)
    Bon_p = {}
    for i in A:
        rig = len([j for j in B if j >= A[i]])
        lef = len([j for j in B if j <= A[i]])
        mp = min(rig,lef)
        p = mp/n
        if p < 0.05 / m:
            Bon_p[i] = (A[i],p,'Ture')
        else:
            Bon_p[i] = (A[i],p,'False')   
            
    position = tuple([i for i in Bon_p])
    index_val = tuple([Bon_p[i][0] for i in Bon_p])
    index_p = tuple([Bon_p[i][1] for i in Bon_p])
    index_state = tuple([Bon_p[i][2] for i in Bon_p])
    dataframe = pd.DataFrame({'position':position,'index_val':index_val,'index_p':index_p,'index_state':index_state})
    dataframe.to_csv("result.csv",index=False,sep=',')
        
    return Bon_p

def random_times(times = None, nettype = 'ER',mode = (1,0)):
    # aa = random_times(times = 3, nettype = 'ER',mode = (1,0))
    start = time.process_time()
    if times == None:
        times = 100

    file_name = str(times) + '_times_' + nettype + '_replace_'+ str(mode[0]) + 'to' + str(mode[1]) + '.xls'    
    random_results = []
    for i in range(times):
        G = graph2matrix(nettype = nettype)
        G = rand_negative(A = G, rate = 0.72)
        temp = main(A = G, mode = mode, threshold = (0,0,0,0,0,0,0,0,0), print_switch = 'off')
        result = tuple(temp.values())        
        random_results.append(result)
        print('当前进度为：', (i+1)/times)
        del G
    
    data_write(file_name, random_results)
    elapsed = (time.process_time() - start)
    print("总计用时",elapsed)
    
def random_times_diag(times = None, nettype = 'ER',mode = (1,0)):
    # aa = random_times(times = 3, nettype = 'ER',mode = (1,0))
    start = time.process_time()
    if times == None:
        times = 100

    file_name = str(times) + '_times_' + nettype + '_replace_'+ str(mode[0]) + 'to' + str(mode[1]) + '_diag.xls'    
    random_results = []
    for i in range(times):
        G = graph2matrix(nettype = nettype)
        G = rand_negative(A = G, rate = 0.72)
        temp = diag_main(A = G, mode = mode, threshold = (0,0,0,0,0,0,0,0,0), print_switch = 'off')
        result = tuple(temp.values())        
        random_results.append(result)
        print('当前进度为：', (i+1)/times)
        del G
    
    data_write(file_name, random_results)
    elapsed = (time.process_time() - start)
    print("总计用时",elapsed)


def data_write(file_path, datas):
    f = xlwt.Workbook()
    sheet1 = f.add_sheet(u'sheet1',cell_overwrite_ok=True) #创建sheet
    
    #将数据写入第 i 行，第 j 列
    i = 0
    for data in datas:
        for j in range(len(data)):
            sheet1.write(j,i,data[j])
        i = i + 1
        
    f.save(file_path) #保存文件

def data_write_csv(file_name, datas):#file_name为写入CSV文件的路径，datas为要写入数据列表
    file_csv = codecs.open(file_name,'a+','utf-8')#追加
    writer = csv.writer(file_csv, delimiter=' ', quotechar=' ', quoting=csv.QUOTE_MINIMAL)
    for data in datas:
        writer.writerow(data)
    print("保存文件成功，处理结束")
    
def result_read(result_mode = 0):
    if type(result_mode) == str:
        data = xlrd.open_workbook('Real_network.xls')
        values=[]
        
        table = data.sheets()[0]
        nrows = table.nrows #行数
        ncols = table.ncols #列数
        if result_mode == 'r0':
            col = 0
        elif result_mode == 'r1':
            col = 1
        elif result_mode == 'r-1':
            col = 2
        
        for x in range(nrows):
            row = table.row_values(x)
            values.append(row[col])
        
        for i in values:
            if '' in values:
                values.remove('')            
        
    if type(result_mode) == int:
        if result_mode == 0:
            data_ER = xlrd.open_workbook('100_times_ER_replace_1to0.xls')
            data_BA = xlrd.open_workbook('100_times_BA_replace_1to0.xls')
            data_WS = xlrd.open_workbook('100_times_WS_replace_1to0.xls')
        elif result_mode == 1:
            data_ER = xlrd.open_workbook('100_times_ER_replace_0to1.xls')
            data_BA = xlrd.open_workbook('100_times_BA_replace_0to1.xls')
            data_WS = xlrd.open_workbook('100_times_WS_replace_0to1.xls')
        elif result_mode == -1:
            data_ER = xlrd.open_workbook('100_times_ER_replace_0to-1.xls')
            data_BA = xlrd.open_workbook('100_times_BA_replace_0to-1.xls')
            data_WS = xlrd.open_workbook('100_times_WS_replace_0to-1.xls')
    
        table_ER = data_ER.sheets()[0]
        nrows = table_ER.nrows #行数
        ncols = table_ER.ncols #列数

        values=[]
        for x in range(nrows):
            row =table_ER.row_values(x)
            for i in range(ncols):
                values.append(row[i])

        table_BA = data_BA.sheets()[0]
        nrows = table_BA.nrows #行数
        ncols = table_BA.ncols #列数
        
        for x in range(nrows):
            row =table_BA.row_values(x)
            for i in range(ncols):
                values.append(row[i])
            
        table_WS = data_WS.sheets()[0]
        nrows = table_WS.nrows #行数
        ncols = table_WS.ncols #列数

        for x in range(nrows):
            row =table_WS.row_values(x)
            for i in range(ncols):
                values.append(row[i])
            
        for i in values:
            if '' in values:
                values.remove('')
        
    return values

def result_read_diag(result_mode = 0):
    if type(result_mode) == str:
        data = xlrd.open_workbook('Real_network.xls')
        values=[]
        
        table = data.sheets()[0]
        nrows = table.nrows #行数
        ncols = table.ncols #列数
        if result_mode == 'r0':
            col = 0
        elif result_mode == 'r1':
            col = 1
        elif result_mode == 'r-1':
            col = 2
        
        for x in range(nrows):
            row = table.row_values(x)
            values.append(row[col])
        
        for i in values:
            if '' in values:
                values.remove('')            
        
    if type(result_mode) == int:
        if result_mode == 0:
            data_ER = xlrd.open_workbook('100_times_ER_replace_1to0_diag.xls')
            data_BA = xlrd.open_workbook('100_times_BA_replace_1to0_diag.xls')
            data_WS = xlrd.open_workbook('100_times_WS_replace_1to0_diag.xls')
        elif result_mode == 1:
            data_ER = xlrd.open_workbook('100_times_ER_replace_0to1_diag.xls')
            data_BA = xlrd.open_workbook('100_times_BA_replace_0to1_diag.xls')
            data_WS = xlrd.open_workbook('100_times_WS_replace_0to1_diag.xls')
        elif result_mode == -1:
            data_ER = xlrd.open_workbook('100_times_ER_replace_0to-1_diag.xls')
            data_BA = xlrd.open_workbook('100_times_BA_replace_0to-1_diag.xls')
            data_WS = xlrd.open_workbook('100_times_WS_replace_0to-1_diag.xls')
    
        table_ER = data_ER.sheets()[0]
        nrows = table_ER.nrows #行数
        ncols = table_ER.ncols #列数

        values=[]
        for x in range(nrows):
            row =table_ER.row_values(x)
            for i in range(ncols):
                values.append(row[i])

        table_BA = data_BA.sheets()[0]
        nrows = table_BA.nrows #行数
        ncols = table_BA.ncols #列数
        
        for x in range(nrows):
            row =table_BA.row_values(x)
            for i in range(ncols):
                values.append(row[i])
            
        table_WS = data_WS.sheets()[0]
        nrows = table_WS.nrows #行数
        ncols = table_WS.ncols #列数

        for x in range(nrows):
            row =table_WS.row_values(x)
            for i in range(ncols):
                values.append(row[i])
            
        for i in values:
            if '' in values:
                values.remove('')
        
    return values

def result_fit(result = 0):
    result_data = result_read(result_mode = result)
    f = fitter.Fitter(result_data)
    f.fit()
    fsummary = f.summary()
    index_0 = fsummary.index[0]
    index = f.fitted_param[index_0]
    print(fsummary)
    print(index)
    

        
def bio_update(A = None,x0 = None,threshold = None):
    """投票更新函数
    调用方法为bio_update(A,x0,threshold)
    A为邻接矩阵，默认值为[[0,1,0],[-1,0,1],[1,-1,0]]
    x0为初始状态,默认值为[1,0,0]
    threshold为投票阈值，默认值为[0,-1,0]
    更新方法：
        网络中存在促进边以及抑制边，当一个节点所有已激活的邻居（状态为1）中：
            ①促进邻居减去抑制邻居高于阈值时，状态变为1；
            ②促进邻居减去抑制邻居低于阈值时，状态变为0；
            ③促进邻居减去抑制邻居等于阈值时，状态不改变；
            
    注意：        
    x0、A、threshold三者大小应匹配，否则报错"""
    if A == None:
        A = (
                (0,1,0),
                (-1,0,1),
                (1,-1,0)
            )    
    if x0 == None:
        x0 = [1,0,0]    
    if threshold == None:
        threshold = [0,-1,0]      
    """以上是默认参数设定环节"""
    A_1 = list(A)
    for i in range(len(A_1)):
        A_1[i] = list(A_1[i])    
    """以上是解压环节""" 
    
    if len(A_1) != len(x0):                                               #A x0 threshold三个的维度必须相同，即网络节点数
        raise AssertionError("邻接矩阵维度与初始状态向量维度不匹配")
    if len(x0) != len(threshold):
        raise AssertionError("初始状态向量维度与阈值向量维度不匹配")
    if len(A_1) != len(threshold):
        raise AssertionError("邻接矩阵维度与阈值向量维度不匹配")           
        
    A_2 = np.matrix(A_1)                       #将输入邻接矩阵列表矩阵化
    x_ini = x0[:]                          #复制保存初始状态
    x0 = np.matrix(x0)                     #将输入的初始状态列表向量化
    
    x_mid = A_2.T * x0.T                     #A第i行第j列代表第i个基因对第j个基因的影响，所以统计i基因收到的影响需将A转置
    x_fin = [0] * len(A_2)                   #初始化最终更新状态向量
    
    for index in range(len(x_mid)):        #对每一个节点进行状态检查，按照上述【更新方法】对状态进行更新
        if int(x_mid[index]) > threshold[index]:
            x_fin[index] = 1
        if int(x_mid[index]) == threshold[index]:
            x_fin[index] = x_ini[index]
        if int(x_mid[index]) < threshold[index]:
            x_fin[index] = 0

    return x_fin            #返回更新后的状态向量

def Basin(A = None,threshold = None):
    """阈值投票更新方式网络的吸引盆求解函数
    A为网络的邻接矩阵
    threshold为网络的更新阈值
    返回值为一字典描述状态空间分布，其中“attactor”键的取值为全部吸引子，“basin”键的取值为吸引子对应的吸引盆
    """
    if A == None:
        A = (
                (-1,0,-1,-1,0,0,0,0,0),
                (0,0,-1,-1,0,0,-1,1,0),
                (0,-1,0,0,0,-1,0,0,0),
                (0,-1,0,0,0,-1,0,0,0),
                (0,-1,0,0,-1,-1,0,0,1),
                (0,0,-1,-1,1,0,0,0,0),
                (0,0,0,0,0,-1,0,0,0),
                (0,0,0,0,0,1,0,0,0),
                (0,0,1,1,0,0,1,-1,-1)
            )  
    if threshold == None:
        threshold = [0,-1,0,0,0,0,0,0,0]       
    """以上是默认参数设定环节"""    
    A_1 = list(A)
    for i in range(len(A_1)):
        A_1[i] = list(A_1[i])    
    """以上是解压环节"""
    
    if len(A_1) != len(threshold):                                       #A threshold的维度必须相同，即网络节点数
        raise AssertionError("邻接矩阵维度与阈值向量维度不匹配")      
    else:
        dim = len(A_1)
        
    state_space = list(itertools.product([0,1],repeat = dim))          #初始化全部状态空间
    state_space_distribution = dict((('attractor',[]),('basin',[])))    #初始化状态空间分布
    for i in range(len(state_space)):                                  #将生成的状态空间中的tuple改为list
        state_space[i] = list(state_space[i])
        
    flag = 0       #flag等于零说明要从现存的状态中随机取一个，等于1说明仍在迭代中
    judge = 0      #judge等于0说明未出现过的吸引子，等于1说明出现过
    x_tem = []     #初始化中间值
    b = []    
    
    while len(state_space) > 0:       #在剩余状态空间非空的情况下，持续迭代
        if flag == 0:                 #如果flag=0，说明已经取到吸引子，则在剩余状态空间中任取一个（默认取第一个）开始新的迭代，并将已经更新的中间状态集合清空
            x0 = state_space[0]
            tem = []
        else:                         #如果flag=1，说明未到吸引子，则继续迭代
            x0 = x_tem
        
        tem.append(x0)                #将x0加入到中间状态集合中
        x_tem = bio_update(A,x0,threshold)   #利用bio_update函数，更新一次x0
        
        if x_tem in tem:                 #如果更新后的x_tem已经出现在过tem集合中，说明取到吸引子，等于tem最后一个说明为不动点，取到tem中间的为极限环
            pos = tem.index(x_tem)       #确定x_tem在tem中出现位置
            if len(state_space_distribution['attractor']) > 0:        #如果已经有吸引子，则需要确认是否为已经出现过的吸引子
                judge = 0                                             #先将judge归为0
                for i in range(len(state_space_distribution['attractor'])):       #检查对比所有已知吸引子
                    if x_tem in state_space_distribution['attractor'][i] or x_tem in state_space_distribution['attractor']:  #如果x_tem出现在了已知极限环，或已知不动点中，更新judge=1，并记录对应吸引子坐标ind
                        judge = 1
                        ind = i
                        
            if judge == 0:   #如果是新的吸引子
                state_space_distribution['attractor'].append(tem[pos:])   #记录吸引子
                state_space_distribution['basin'].append(tem)   #记录吸引盆
                for i in range(len(tem)):      #从剩余状态空间中移除所有中间过程出现过的状态
                    if tem[i] in state_space:
                        state_space.remove(tem[i])
                del tem          #删除中间状态
                flag = 0         #将flag置为0，开始新一轮迭代
            else:    #否则，即为已经出现过得吸引子不记录吸引子，直接将中间状态存储至basin中
                if len(tem) == 1:                 #如果只有一个状态则直接记录       
                    state_space_distribution['basin'][ind].append(tem)   
                else:       #如果有多个状态，则把list拆开一个一个记录
                    for i in range(len(tem)):
                        state_space_distribution['basin'][ind].append(tem[i])
                        
                for i in range(len(tem)):  #从剩余状态空间中移除所有中间过程出现过的状态
                    if tem[i] in state_space:
                        state_space.remove(tem[i])
                del tem #删除中间状态
                del ind #删除吸引子位置坐标
                flag = 0 #重置flag开始新一轮的迭代
        else:
            flag = 1  #如果x_tem不在tem中，说明不是吸引子，直接迭代下一轮
        

    for i in range(len(state_space_distribution['basin'])):   #将basin去重（吸引子会在basin中出现多次）
        b.append(list_duplicate_removal(state_space_distribution['basin'][i]))   #利用list_duplicate_removal函数对list的list去重
                
    state_space_distribution.update(basin=b)   #更新无重复的basin，basin的总长度和应为2^dim
        
    return state_space_distribution  #返回状态空间的分布

def list_duplicate_removal(L = None):        #将列表的列表去重的函数
    """将列表的列表去重，仅限去除一层"""
    if L == None:
        L = [
                [1,0,0],
                [1,0,0],
                [0,1,0]
            ]
    """以上是默认参数设定环节"""    
        
    b_list = []              #建立一个新列表
    for i in L:              
        if i not in b_list:        #扫描A，将A中不存在于B中的元素放在B中
            b_list.append(i)
    return b_list             #返回B

def A_replace(A = None,position = None,interaction = None):
    """
    替换邻接矩阵函数
    可以替换邻接矩阵A中位置为position的连边为interaction
    其中A为邻接矩阵
    position为对应位置，从[0,0]开始计数
    interaction为-1,0,1中的某一值
    """
    if A == None:
        A = (
                (-1,0,-1,-1,0,0,0,0,0),
                (0,0,-1,-1,0,0,-1,1,0),
                (0,-1,0,0,0,-1,0,0,0),
                (0,-1,0,0,0,-1,0,0,0),
                (0,-1,0,0,-1,-1,0,0,1),
                (0,0,-1,-1,1,0,0,0,0),
                (0,0,0,0,0,-1,0,0,0),
                (0,0,0,0,0,1,0,0,0),
                (0,0,1,1,0,0,1,-1,-1)
            )      #使用不可变默认参数，防止多次调用时默认参数改变的问题
    if position == None: #使用不可变默认参数，防止多次调用时默认参数改变的问题
        position = [0,0]        
    if interaction == None: #使用不可变默认参数，防止多次调用时默认参数改变的问题
        interaction = 0
    """以上是默认参数设定环节"""    
    A_1 = list(A)
    for i in range(len(A_1)):
        A_1[i] = list(A_1[i])    
    """以上是解压环节"""
        
    if len(position) != 2:
         raise AssertionError("替换位置坐标position应为一个二维列表")
    if interaction != 0 and interaction != -1 and interaction != 1:
        raise AssertionError("替换interaction应为0或正负1中的一个值")
    
    row = position[0]
    col = position[1]   
    A_1[row][col] = interaction #将A中position对应位置的状态改为interaction对应的状态
    
    for i in range(len(A_1)):
        A_1[i] = tuple(A_1[i])
    """压缩过程"""
    A_2 = tuple(A_1)
    return A_2

def graph_from_A(A = None,draw = 1, Nodes = None):
    """
    通过邻接矩阵A构建网络
    A为邻接矩阵
    draw等于1时画图，draw=0时不画图
    Nodes可以为节点提供名字，维数应与理解矩阵对应
    需要调用networkx库
    """
    
    if A == None:
        A = (
                (-1,0,-1,-1,0,0,0,0,0),
                (0,0,-1,-1,0,0,-1,1,0),
                (0,-1,0,0,0,-1,0,0,0),
                (0,-1,0,0,0,-1,0,0,0),
                (0,-1,0,0,-1,-1,0,0,1),
                (0,0,-1,-1,1,0,0,0,0),
                (0,0,0,0,0,-1,0,0,0),
                (0,0,0,0,0,1,0,0,0),
                (0,0,1,1,0,0,1,-1,-1)
            )      #使用不可变默认参数，防止多次调用时默认参数改变的问题
        
    if Nodes == None:
        Nodes = ['SK','Cdc2/Cdc13','Ste9','Rum1','Slp1','Cdc2/Cdc13*','Wee1/Mik1','Cdc25','PP'] #使用不可变默认参数，防止多次调用时默认参数改变的问题
        
    """以上为默认参数调用环节"""
    A_1 = list(A)
    for i in range(len(A_1)):
        A_1[i] = list(A_1[i])    
    """以上是解压环节"""
    G = nx.DiGraph() #产生一张空的图
    
    for i in Nodes: #添加节点
        G.add_node(i)
     
    for i in range(len(A_1)): #添加连边
        for j in range(len(A_1)):
            if A[i][j] != 0:
                G.add_weighted_edges_from([(Nodes[i], Nodes[j],A_1[i][j])])
    
    if draw == 1: #如果draw=1,就把G画出来
        nx.draw(G,pos = nx.circular_layout(G))  
        plt.show()
    """
    可用layout包括：
    *bipartite_layout(G, nodes[, align, scale, . . . ]), Position nodes in two straight lines.
    circular_layout(G[, scale, center, dim]), Position nodes on a circle.
    kamada_kawai_layout(G[, dist, pos, weight, . . . ]), Position nodes using Kamada-Kawai path-length costfunction.
    planar_layout(G[, scale, center, dim]), Position nodes without edge intersections.
    random_layout(G[, center, dim, seed]), Position nodes uniformly at random in the unit square.
    *rescale_layout(pos[, scale]), Returns scaled position array to (-scale, scale) in all axes.
    spring_layout(G[, k, pos, fixed, . . . ]), Position nodes using Fruchterman-Reingold forcedirected algorithm.
    shell_layout(G[, nlist, scale, center, dim]), Position nodes in concentric circles.
    spectral_layout(G[, weight, scale, center, dim]), Position nodes using the eigenvectors of the graph Laplacian.
    """
    return G

def list_size(A = None):
    """
    计算状态空间划分列表A中的全部状态数
    A的格式应为一个列表，
    列表中划分为同一分支的状态为一个列表，
    每一个状态也为一个列表，
    即使一个分支中仅有一个状态，同样要分为列表的列表，具体示例如下：
    A = [[[0,0],[0,1]],
         [[1,0]],
         [[1,1]]]
    """
    if A == None:
        A = [[[0,0],[0,1]],
             [[1,0]],
             [[1,1]]]
    """以上为默认参数调用环节"""   
    
    length = len(A) #记录A中分支数目
    size_A = 0 #初始化大小
    for i in range(length):  #扫描每一个分支
        tem = len(A[i]) #tem记录分支中状态数
        size_A += tem #将状态数加在大小中
        
    return size_A #返回大小

def state_arrary(A = None):
    """
    状态向量函数，记录状态空间划分列表中的状态分布向量
    假设共有N个状态（其中 N = 2^n,n为节点数量）
    则 dim(state_arrary) = N(N-1)/2
    记录了状态空间中任意两个状态是否处在同一个分支中，1代表在，2代表不在
    列表从 (0,0,...,0) 按照二进制大小顺序依次排列至 (1,1,...,1)
    state_arrary[0]记录(0,0,...,0,0)与(0,0,...,0,1)是否在同一分支中
    state_arrary[1]记录(0,0,...,0,0)与(0,0,...,0,2)是否在同一分支中
    ...
    state_arrary[N-3]记录(1,1,...,0,0)与(1,1,...,1,0)是否在同一分支中
    state_arrary[N-2]记录(1,1,...,0,0)与(1,1,...,1,1)是否在同一分支中
    state_arrary[N-1]记录(1,1,...,1,0)与(1,1,...,1,1)是否在同一分支中
    """
    if A == None:
        A = [[[0,0],[0,1]],[[1,0]],[[1,1]]]
    """以上为默认参数调用环节"""    
    
    dim = len(A[0][0])
    state_space = list(itertools.product([0,1],repeat = dim))          #初始化全部状态空间
    arrary_A = [] #初始记录向量
    
    for i in range(len(state_space)): #扫描全部初始状态——————
        tem_1 = list(state_space[i]) #记录状态1
        pos_1 = None #清空状态1位置
        for k in range(len(A)): #查找状态1在第几个分支里
            if tem_1 in A[k]: 
                pos_1 = k     #记录状态1分支   
        for j in range(len(state_space) - i - 1): #及其之后的状态与其是否为同一分支之中
            tem_2 = list(state_space[i+j+1]) #记录状态2
            pos_2 = None #清空状态2位置
            for k in range(len(A)): #查找状态2在第几个分支里
                if tem_2 in A[k]:
                    pos_2 = k #记录状态2分支   
            if pos_1 == pos_2:
                arrary_A.append(1) #如果状态1与状态2在同一分支中，向arrary中添加1
            else:
                arrary_A.append(0) #如果状态1与状态2不在同一分支中，向arrary中添加0
    return arrary_A #返回向量

def Adjacency_marix(a = None):
    if a == None:
        a = ((0,1,1),(1,0,1))
        
    size_a = max(max(a[0]),max(a[1])) + 1
    A_0 = [0] * size_a
    A_mat = [A_0] * size_a
    A = np.array(A_mat)
    
    for i in range(len(a[0])):
        A[a[0][i]][a[1][i]] = 1

    return A

def Reachability_matrix(A = None):
    
    #start = time.process_time()    
    if A.any() == None:
        A = np.array([[0,1,0,0],[0,0,1,0],[0,0,0,1],[0,0,0,1]])    

        
    size_A = len(A)
    R = A
    
    for i in range(size_A):
        for j in range(size_A):
            for k in range(size_A):
                if R[i][j]*R[j][k] != 0:
                    R[i][k] = 1    
                    
    #elapsed = (time.process_time() - start)
    #print("总计用时",elapsed)
                       
    return R
    
def DReSS(A = np.array([[1,0,0],[0,1,0],[0,0,0]]), B = np.array([[1,0,1],[0,1,0],[1,0,1]])):

    simi = 0 #初始化相似度
    for i in range(len(A)):
        for j in range(len(A)):
            if A[i][j] != B[i][j]:
                simi += 1 #每发现一个相同的位置，相似度加一
    simi_std = simi / (len(A) * len(A)) #标准化相似度，除以向量长度，是的标准相似度最大为1，最小为0
    return simi_std #返回三个值，相似度，向量长度，标准化相似度

def haming_diag(A = (1,0,1),B = (0,1,1)):
    
    simi = 0
    for i in range(len(A)):
        if A[i] != B[i]:
            simi += 1
            
    simi_std = simi / len(A)
    return simi_std

def haming_simi_basin(Basin_A = None, Basin_B = None):
    """
    计算两个状态分布空间对应的状态分布向量的haming相似度
    输入两个状态分布空间
    计算两个状态分布向量
    进而计算相似度
    返回三个值：
    相似维数，状态空间大小，标准化相似度
    标准化相似度为0-1之间取值，1为完全相同，0为完全不同
    """
    if Basin_A == None:
        Basin_A = [[[0,0],[0,1]],[[1,0]],[[1,1]]]
    if Basin_B == None:
        Basin_B = [[[0,0],[0,1],[1,0]],[[1,1]]]
    """以上为默认参数调用环节"""
    
    arrary_A = Reachability_matrix(Adjacency_marix(Basin_A)) #生成状态空间A的分布向量
    arrary_B = Reachability_matrix(Adjacency_marix(Basin_B)) #生成状态空间B的分布向量
    simi = 0 #初始化相似度
    for i in range(len(arrary_A)):
        for j in range(len(arrary_A)):
            if arrary_A[i][j] != arrary_B[i][j]:
                simi += 1 #每发现一个不同的位置，距离加一
    
    simi_std = simi / (len(arrary_A)^2) #标准化相似度，除以向量长度，是的标准相似度最大为1，最小为0
    return simi_std #返回三个值，相似度，向量长度，标准化相似度

def mat_find_zero(A = None):
    """
    查找邻接矩阵中零元素的位置
    输入的必须为方阵，否则报错
    返回值为一个数组tuple，里面每一个元素为对应矩阵中的坐标
    从[0,0]开始
    """
    if A == None:
        A = (
                (-1,0,-1,-1,0,0,0,0,0),
                (0,0,-1,-1,0,0,-1,1,0),
                (0,-1,0,0,0,-1,0,0,0),
                (0,-1,0,0,0,-1,0,0,0),
                (0,-1,0,0,-1,-1,0,0,1),
                (0,0,-1,-1,1,0,0,0,0),
                (0,0,0,0,0,-1,0,0,0),
                (0,0,0,0,0,1,0,0,0),
                (0,0,1,1,0,0,1,-1,-1)
            )      #使用不可变默认参数，防止多次调用时默认参数改变的问题
    """以上为默认参数调用环节"""
    A_1 = list(A)
    for i in range(len(A_1)):
        A_1[i] = list(A_1[i])    
    """以上是解压环节"""
    
    size_1 = len(A_1) #计算A矩阵行数
    for i in range(size_1): #检索每一行元素数是否与总行数相同
        if len(A_1[i]) != size_1:
            raise AssertionError("输入矩阵应为方阵")
    
    pos = [] #初始化位置列表
    for i in range(size_1): #检查每一行
        for j in range(size_1): #检查每一列
            if A_1[i][j] == 0: #如果元素等于0
                pos.append((i,j)) #记录位置
                
    pos = tuple(pos) #将位置列表转化为不可变的tuple
    
    return pos #返回pos

def mat_find_nonzero(A = None):
    """
    查找邻接矩阵中非零元素的位置
    输入的必须为方阵，否则报错
    返回值为一个数组tuple，里面每一个元素为对应矩阵中的坐标
    从[0,0]开始
    """
    if A == None:
        A = (
                (-1,0,-1,-1,0,0,0,0,0),
                (0,0,-1,-1,0,0,-1,1,0),
                (0,-1,0,0,0,-1,0,0,0),
                (0,-1,0,0,0,-1,0,0,0),
                (0,-1,0,0,-1,-1,0,0,1),
                (0,0,-1,-1,1,0,0,0,0),
                (0,0,0,0,0,-1,0,0,0),
                (0,0,0,0,0,1,0,0,0),
                (0,0,1,1,0,0,1,-1,-1)
            )      #使用不可变默认参数，防止多次调用时默认参数改变的问题
    """以上为默认参数调用环节"""
    A_1 = list(A)
    for i in range(len(A_1)):
        A_1[i] = list(A_1[i])    
    """以上是解压环节"""
    
    size_1 = len(A_1) #计算A矩阵行数
    for i in range(size_1): #检索每一行元素数是否与总行数相同
        if len(A_1[i]) != size_1:
            raise AssertionError("输入矩阵应为方阵")
    
    pos = [] #初始化位置列表
    for i in range(size_1):  #检查每一行
        for j in range(size_1): #检查每一列
            if A_1[i][j] != 0: #如果元素不等于0
                pos.append((i,j)) #记录位置
                
    pos = tuple(pos) #将位置列表转化为不可变的tuple
     
    return pos #返回pos

def graph2matrix(G = None,nettype = 'ER'):
    """
    
    """
    if G == None:
        if nettype == 'ER':
            G = nx.random_graphs.erdos_renyi_graph(9, 0.3)
        elif nettype == 'WS':
            G = nx.random_graphs.watts_strogatz_graph(9, 3, 0.3)
        elif nettype == 'BA':
            G = nx.random_graphs.barabasi_albert_graph(9, 1)
        

    G_mid = pd.DataFrame(np.zeros([len(G.nodes()),len(G.nodes())]),index = G.nodes(),columns = G.nodes())
    for col_label,row_label in G.edges():
        G_mid.loc[col_label,row_label] = 1 
        G_mid.loc[row_label,col_label] = 1
    
    G_mat = []
    G_size = len(G_mid)
    
    for i in range(G_size):
        G_mat.append([])
        
    for i in range(G_size):
        for j in range(G_size):
            G_mat[i].append(int(G_mid[i][j]))
    
    for i in range(G_size):
        G_mat[i] = tuple(G_mat[i])
        
    G_mat = tuple(G_mat)
    return G_mat

def rand_negative(A = None, rate = 0.72):
    """将矩阵A中的连边以概率rate变为抑制连边
    原始矩阵中有18个抑制边 7个促进边"""
    if A == None:
        A = (
                (-1,0,-1,-1,0,0,0,0,0),
                (0,0,-1,-1,0,0,-1,1,0),
                (0,-1,0,0,0,-1,0,0,0),
                (0,-1,0,0,0,-1,0,0,0),
                (0,-1,0,0,-1,-1,0,0,1),
                (0,0,-1,-1,1,0,0,0,0),
                (0,0,0,0,0,-1,0,0,0),
                (0,0,0,0,0,1,0,0,0),
                (0,0,1,1,0,0,1,-1,-1)
            )      #使用不可变默认参数，防止多次调用时默认参数改变的问题
    """以上为默认参数调用环节"""
    A_1 = list(A)
    for i in range(len(A_1)):
        A_1[i] = list(A_1[i])    
    """以上是解压环节"""
    pos = mat_find_nonzero(A = A)
    n = len(pos)
    for i in range(n):
        r = random.random()
        if r < rate:
            A_1[pos[i][0]][pos[i][1]] = -1
            
    for i in range(len(A_1)):
        A_1[i] = tuple(A_1[i])
    """压缩过程"""
    A_2 = tuple(A_1)
    return A_2

def Basin_with_plot(A = None,threshold = None):
    """阈值投票更新方式网络的吸引盆求解函数
    A为网络的邻接矩阵
    threshold为网络的更新阈值
    返回值为一字典描述状态空间分布，其中“attactor”键的取值为全部吸引子，“basin”键的取值为吸引子对应的吸引盆
    """
    if A == None:
        A = (
                (-1,0,-1,-1,0,0,0,0,0),
                (0,0,-1,-1,0,0,-1,1,0),
                (0,-1,0,0,0,-1,0,0,0),
                (0,-1,0,0,0,-1,0,0,0),
                (0,-1,0,0,-1,-1,0,0,1),
                (0,0,-1,-1,1,0,0,0,0),
                (0,0,0,0,0,-1,0,0,0),
                (0,0,0,0,0,1,0,0,0),
                (0,0,1,1,0,0,1,-1,-1)
            )  
    if threshold == None:
        threshold = [0,-1,0,0,0,0,0,0,0]       
    """以上是默认参数设定环节"""    
    A_1 = list(A)
    for i in range(len(A_1)):
        A_1[i] = list(A_1[i])    
    """以上是解压环节"""
    
    if len(A_1) != len(threshold):                                       #A threshold的维度必须相同，即网络节点数
        raise AssertionError("邻接矩阵维度与阈值向量维度不匹配")      
    else:
        dim = len(A_1)
        
    state_space = list(itertools.product([0,1],repeat = dim))          #初始化全部状态空间
    state_space_distribution = dict((('attractor',[]),('basin',[])))    #初始化状态空间分布
    for i in range(len(state_space)):                                  #将生成的状态空间中的tuple改为list
        state_space[i] = list(state_space[i])
    Source_index_list = []
    Target_index_list = []
    state_space_tuple = list(itertools.product([0,1],repeat = dim))
        
    flag = 0       #flag等于零说明要从现存的状态中随机取一个，等于1说明仍在迭代中
    judge = 0      #judge等于0说明未出现过的吸引子，等于1说明出现过
    x_tem = []     #初始化中间值
    b = []    
    R_0 = [0] * len(state_space_tuple)
    R_mat = [R_0] * len(state_space_tuple)
    R = np.array(R_mat)
    
    while len(state_space) > 0:       #在剩余状态空间非空的情况下，持续迭代
        if flag == 0:                 #如果flag=0，说明已经取到吸引子，则在剩余状态空间中任取一个（默认取第一个）开始新的迭代，并将已经更新的中间状态集合清空
            x0 = state_space[0]
            tem = []
            tem_index = []
        else:                         #如果flag=1，说明未到吸引子，则继续迭代
            x0 = x_tem
        
        Source_index = state_space_tuple.index(tuple(x0))
        Source_index_list.append(Source_index)
        tem.append(x0)                #将x0加入到中间状态集合中
        tem_index.append(Source_index)
        x_tem = bio_update(A,x0,threshold)   #利用bio_update函数，更新一次x0
        Target_index = state_space_tuple.index(tuple(x_tem))
        Target_index_list.append(Target_index)
        

        
        if x_tem in tem:                 #如果更新后的x_tem已经出现在过tem集合中，说明取到吸引子，等于tem最后一个说明为不动点，取到tem中间的为极限环
            pos = tem.index(x_tem)       #确定x_tem在tem中出现位置
            if len(state_space_distribution['attractor']) > 0:        #如果已经有吸引子，则需要确认是否为已经出现过的吸引子
                judge = 0                                             #先将judge归为0
                for i in range(len(state_space_distribution['attractor'])):       #检查对比所有已知吸引子
                    if x_tem in state_space_distribution['attractor'][i] or x_tem in state_space_distribution['attractor']:  #如果x_tem出现在了已知极限环，或已知不动点中，更新judge=1，并记录对应吸引子坐标ind
                        judge = 1
                        ind = i
                        
            if judge == 0:   #如果是新的吸引子
                state_space_distribution['attractor'].append(tem[pos:])   #记录吸引子
                state_space_distribution['basin'].append(tem)   #记录吸引盆
                for i in range(len(tem)):      #从剩余状态空间中移除所有中间过程出现过的状态
                    if tem[i] in state_space:
                        state_space.remove(tem[i])
                
                tem_index.append(Target_index)
                for r_0 in range(len(tem_index)-1):
                    for r_n in range(len(tem) - r_0):
                        R[tem_index[r_0]][tem_index[r_0 + r_n + 1]] = 1
                
                del tem          #删除中间状态
                del tem_index
                flag = 0         #将flag置为0，开始新一轮迭代
            else:    #否则，即为已经出现过得吸引子不记录吸引子，直接将中间状态存储至basin中
                if len(tem) == 1:                 #如果只有一个状态则直接记录       
                    state_space_distribution['basin'][ind].append(tem)   
                else:       #如果有多个状态，则把list拆开一个一个记录
                    for i in range(len(tem)):
                        state_space_distribution['basin'][ind].append(tem[i])
                        
                for i in range(len(tem)):  #从剩余状态空间中移除所有中间过程出现过的状态
                    if tem[i] in state_space:
                        state_space.remove(tem[i])
                        
                tem_index.append(Target_index)
                for r_0 in range(len(tem_index)-1):
                    for r_n in range(len(tem) - r_0):
                        R[tem_index[r_0]][tem_index[r_0 + r_n + 1]] = 1
                        
                del tem #删除中间状态
                del tem_index
                del ind #删除吸引子位置坐标
                flag = 0 #重置flag开始新一轮的迭代
        else:
            flag = 1  #如果x_tem不在tem中，说明不是吸引子，直接迭代下一轮
        

    for i in range(len(state_space_distribution['basin'])):   #将basin去重（吸引子会在basin中出现多次）
        b.append(list_duplicate_removal(state_space_distribution['basin'][i]))   #利用list_duplicate_removal函数对list的list去重
                
    state_space_distribution.update(basin=b)   #更新无重复的basin，basin的总长度和应为2^dim
    Source_index_list = tuple(Source_index_list)
    Target_index_list = tuple(Target_index_list)
    state_space_network = tuple([Source_index_list,Target_index_list])
    EdgeNumber = len(state_space_network[0])
    Typelabel = ['Directed'] * EdgeNumber
    ID_num = []
    for i in range(EdgeNumber):
        ID_num.append(i+dim)
    dataframe = pd.DataFrame({'Source':Source_index_list,'Target':Target_index_list,'Id':ID_num,'Type':Typelabel})
    dataframe.to_csv("edges.csv",index=False,sep=',')
    Id_2 = list(range(len(state_space_tuple)))
    dataframe = pd.DataFrame({'Id':Id_2,'Label':state_space_tuple})
    dataframe.to_csv("nodes.csv",index=False,sep=',')   
    return state_space_network, R  #返回状态空间的分布


def trace(A = None):
    if A.any == None:
        A = np.array([[1,0,0],[0,1,0],[0,0,0]])
    
    tr = 0
        
    for i in range(len(A)):
        tr += A[i][i]
        
    return tr

def diag_array(A = np.array([[1,0,0],[0,1,0],[0,0,0]])):
    
    diag = []
        
    for i in range(len(A)):
        diag.append(A[i][i])
        
    diag = tuple(diag)
        
    return diag
