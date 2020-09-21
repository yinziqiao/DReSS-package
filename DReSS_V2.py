# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 16:57:23 2019

@author: 殷子樵 Yin Ziqiao


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
            print((i+1)/len(range(len(pos))))
            
    elapsed = (time.process_time() - start)
    print("总计用时",elapsed)
    return pos_haming_simi_basin

def Bonferroni(A = None, B = None):
    """Function for Bonferroni Test with two parameters
    A is the test set with position information
    B is the control set with only values from randomized controlled trial
    result return to a dictionary with Position, p-value, significance
    and write it to an csv file named 'result.csv'
    """
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
    """Function for DReSS randomized controlled trial
    There are 3 parameters
    times control the random times for randomized controlled trial
    nettype can be chosen from followed:
                        'ER' for ER random networks
                        'WS' for WS small world networks
                        'BA' for BA scale free networks
    mode can be chosen from follows:
                     (1,0) --- perturbation that delete existing interactions
                     (0,1) --- perturbation that add activation interactions
                     (0,-1) --- perturbation that add inhibition interactions
    result return to control set with only DReSS values from randomized controlled trial
    and will be written into a xls file named xxx_times_ER/WS/BA_replace_mode[0]tomode[1].xls
    """
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
        print((i+1)/times)
        del G
    
    data_write(file_name, random_results)
    elapsed = (time.process_time() - start)
    print("总计用时",elapsed)
    
def random_times_diag(times = None, nettype = 'ER',mode = (1,0)):
    """Function for diagDReSS randomized controlled trial
    There are 3 parameters
    times control the random times for randomized controlled trial
    nettype can be chosen from followed:
                        'ER' for ER random networks
                        'WS' for WS small world networks
                        'BA' for BA scale free networks
    mode can be chosen from follows:
                     (1,0) --- perturbation that delete existing interactions
                     (0,1) --- perturbation that add activation interactions
                     (0,-1) --- perturbation that add inhibition interactions
    result return to control set with only diagDReSS values from randomized controlled trial
    and will be written into a xls file named xxx_times_ER/WS/BA_replace_mode[0]tomode[1].xls
    """
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
        print((i+1)/times)
        del G
    
    data_write(file_name, random_results)
    elapsed = (time.process_time() - start)
    print(elapsed)


def data_write(file_path, datas):
    """ function for data writing"""
    f = xlwt.Workbook()
    sheet1 = f.add_sheet(u'sheet1',cell_overwrite_ok=True) #create a sheet
    
    #write data into row i column j
    i = 0
    for data in datas:
        for j in range(len(data)):
            sheet1.write(j,i,data[j])
        i = i + 1
        
    f.save(file_path) #save file

def data_write_csv(file_name, datas):
    """function for data csv writing"""
    file_csv = codecs.open(file_name,'a+','utf-8') #add
    writer = csv.writer(file_csv, delimiter=' ', quotechar=' ', quoting=csv.QUOTE_MINIMAL)
    for data in datas:
        writer.writerow(data)
    print("sace successfully")
    
def result_read(result_mode = 0):
    """read function DReSS result and turn it into format that can be read by function Bonferroni
    result_mode can be chosen from follows:
        'r0' for reading real_network deleting interaction results
        'r1' for reading real_network adding a actication interaction results
        'r-1' for reading real_network adding a inhibition interaction results
        '0' for reading randomized trials deleting interaction results
        '1' for reading randomized trials adding a actication interaction results
        '-1' for reading randomized trials adding a inhibition interaction results
    """
    if type(result_mode) == str:
        data = xlrd.open_workbook('Real_network.xls')
        values=[]
        
        table = data.sheets()[0]
        nrows = table.nrows #row
        ncols = table.ncols #column
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
        nrows = table_ER.nrows #row
        ncols = table_ER.ncols #column

        values=[]
        for x in range(nrows):
            row =table_ER.row_values(x)
            for i in range(ncols):
                values.append(row[i])

        table_BA = data_BA.sheets()[0]
        nrows = table_BA.nrows #row
        ncols = table_BA.ncols #column
        
        for x in range(nrows):
            row =table_BA.row_values(x)
            for i in range(ncols):
                values.append(row[i])
            
        table_WS = data_WS.sheets()[0]
        nrows = table_WS.nrows #row
        ncols = table_WS.ncols #column

        for x in range(nrows):
            row =table_WS.row_values(x)
            for i in range(ncols):
                values.append(row[i])
            
        for i in values:
            if '' in values:
                values.remove('')
        
    return values

def result_read_diag(result_mode = 0):
    """read function diagDReSS result and turn it into format that can be read by function Bonferroni
    result_mode can be chosen from follows:
        'r0' for reading real_network deleting interaction results
        'r1' for reading real_network adding a actication interaction results
        'r-1' for reading real_network adding a inhibition interaction results
        '0' for reading randomized trials deleting interaction results
        '1' for reading randomized trials adding a actication interaction results
        '-1' for reading randomized trials adding a inhibition interaction results
    """
    if type(result_mode) == str:
        data = xlrd.open_workbook('Real_network.xls')
        values=[]
        
        table = data.sheets()[0]
        nrows = table.nrows #row
        ncols = table.ncols #column
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
        nrows = table_ER.nrows #row
        ncols = table_ER.ncols #column

        values=[]
        for x in range(nrows):
            row =table_ER.row_values(x)
            for i in range(ncols):
                values.append(row[i])

        table_BA = data_BA.sheets()[0]
        nrows = table_BA.nrows #row
        ncols = table_BA.ncols #column
        
        for x in range(nrows):
            row =table_BA.row_values(x)
            for i in range(ncols):
                values.append(row[i])
            
        table_WS = data_WS.sheets()[0]
        nrows = table_WS.nrows #row
        ncols = table_WS.ncols #column

        for x in range(nrows):
            row =table_WS.row_values(x)
            for i in range(ncols):
                values.append(row[i])
            
        for i in values:
            if '' in values:
                values.remove('')
        
    return values

def result_fit(result = 0):
    """Function for result fitting
    """
    result_data = result_read(result_mode = result)
    f = fitter.Fitter(result_data)
    f.fit()
    fsummary = f.summary()
    index_0 = fsummary.index[0]
    index = f.fitted_param[index_0]
    print(fsummary)
    print(index)
    

        
def bio_update(A = None,x0 = None,threshold = None):
    """state updating function for one step update
    There are 3 parameters
    A is the adjacency matrix for origin system
    x0 is the current state
    threshold is the updating threshold
    function returns to x0's next time step state
    
    Notice: A x0 threshold must with same dimension
        Otherwise function return to error
    """
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
     """Parameter setting part ends here"""
    A_1 = list(A)
    for i in range(len(A_1)):
        A_1[i] = list(A_1[i])    
    """unzip part ends here"""
    
    if len(A_1) != len(x0):                                               #A x0 threshold must with same dimension
        raise AssertionError("A x0 must with same dimension")
    if len(x0) != len(threshold):
        raise AssertionError("x0 threshold must with same dimension")
    if len(A_1) != len(threshold):
        raise AssertionError("A threshold must with same dimension")           
        
    A_2 = np.matrix(A_1)                     
    x_ini = x0[:]                          #save initial state
    x0 = np.matrix(x0)                    
    
    x_mid = A_2.T * x0.T                     #row i column j in A means gene i 's interaction to j, therefore A need to be transpose
    x_fin = [0] * len(A_2)                   
    
    for index in range(len(x_mid)):        #state check
        if int(x_mid[index]) > threshold[index]:
            x_fin[index] = 1
        if int(x_mid[index]) == threshold[index]:
            x_fin[index] = x_ini[index]
        if int(x_mid[index]) < threshold[index]:
            x_fin[index] = 0

    return x_fin            #返回更新后的状态向量

def list_duplicate_removal(L = None):        
    """function for remove list duplicate values"""
    if L == None:
        L = [
                [1,0,0],
                [1,0,0],
                [0,1,0]
            ]
     """Parameter setting part ends here"""   
        
    b_list = []             
    for i in L:              
        if i not in b_list:        
            b_list.append(i)
    return b_list            

def A_replace(A = None,position = None,interaction = None):
    """
    Function for replace one positon's value with interaction value
    A is the adjacency matrix
    position is the replace target positon 
    interaction is the new value to replace the origin value
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
    if position == None: 
        position = [0,0]        
    if interaction == None: 
        interaction = 0
     """Parameter setting part ends here"""
    A_1 = list(A)
    for i in range(len(A_1)):
        A_1[i] = list(A_1[i])    
    """unzip part ends here"""
        
    if len(position) != 2:
         raise AssertionError("position should be 2 dimension")
    if interaction != 0 and interaction != -1 and interaction != 1:
        raise AssertionError("interaction should be either 0 or -1 or 1")
    
    row = position[0]
    col = position[1]   
    A_1[row][col] = interaction
    
    for i in range(len(A_1)):
        A_1[i] = tuple(A_1[i])
    """zip part"""
    A_2 = tuple(A_1)
    return A_2

def graph_from_A(A = None,draw = 1, Nodes = None):
    """
    Function to turn agjacency matrix A into a graph which can be read by networkx package
    draw can be chosen from follows:
        1 plot network
        0 do not plot network
    Nodes can provide names for each node
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
        
    if Nodes == None:
        Nodes = ['SK','Cdc2/Cdc13','Ste9','Rum1','Slp1','Cdc2/Cdc13*','Wee1/Mik1','Cdc25','PP'] 
        
     """Parameter setting part ends here"""
    A_1 = list(A)
    for i in range(len(A_1)):
        A_1[i] = list(A_1[i])    
    """unzip part ends here"""
    G = nx.DiGraph() #generalized a new graph
    
    for i in Nodes: #add nodes
        G.add_node(i)
     
    for i in range(len(A_1)): #add edges
        for j in range(len(A_1)):
            if A[i][j] != 0:
                G.add_weighted_edges_from([(Nodes[i], Nodes[j],A_1[i][j])])
    
    if draw == 1: #if draw = 1, plot
        nx.draw(G,pos = nx.circular_layout(G))  
        plt.show()
    """
    layout including:
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

def Adjacency_marix(a = None):
    """Function to get adjacency matrix from a adjacency list
    input a is a 2 row tuple
    first row is source
    second row is target
    function return to a's corresponding adjacency matrix
    """
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
    """
    Function to get modified Reachability matrix R
    input A is the adjacency matrix of state space
    function returns to the modified reachability matrix of state space with adjacency matrix A
    """
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
    #print(elapsed)
                       
    return R
    
def DReSS(A = np.array([[1,0,0],[0,1,0],[0,0,0]]), B = np.array([[1,0,1],[0,1,0],[1,0,1]])):
    """
    Function for calculating DReSS value
    A and B are two modified reachablity matrix
    """
    simi = 0
    for i in range(len(A)):
        for j in range(len(A)):
            if A[i][j] != B[i][j]:
                simi += 1 
    simi_std = simi / (len(A) * len(A)) you 
    return simi_std 

def haming_diag(A = (1,0,1),B = (0,1,1)):
    """
    function for calculating diagDReSS
    A and B are diag array of two modified reachablity matrix
    """
    simi = 0
    for i in range(len(A)):
        if A[i] != B[i]:
            simi += 1
            
    simi_std = simi / len(A)
    return simi_std

def mat_find_zero(A = None):
    """
    Function for find all position of zero values in an adjacency matrix A 
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
     """Parameter setting part ends here"""
    A_1 = list(A)
    for i in range(len(A_1)):
        A_1[i] = list(A_1[i])    
    """unzip part ends here"""
    
    size_1 = len(A_1) 
    for i in range(size_1): 
        if len(A_1[i]) != size_1:
            raise AssertionError("A should be N*N matrix")
    
    pos = [] 
    for i in range(size_1): 
        for j in range(size_1): 
            if A_1[i][j] == 0: 
                pos.append((i,j)) 
                
    pos = tuple(pos) 
    
    return pos 

def mat_find_nonzero(A = None):
    """
    Function for find all position of non-zero values in an adjacency matrix A 
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
     """Parameter setting part ends here"""
    A_1 = list(A)
    for i in range(len(A_1)):
        A_1[i] = list(A_1[i])    
    """unzip part ends here"""
    
    size_1 = len(A_1) 
    for i in range(size_1): 
        if len(A_1[i]) != size_1:
            raise AssertionError("A should be N*N matrix")
    
    pos = [] 
    for i in range(size_1):  
        for j in range(size_1): 
            if A_1[i][j] != 0: 
                pos.append((i,j)) 
                
    pos = tuple(pos)
     
    return pos 

def graph2matrix(G = None,nettype = 'ER'):
    """
    function to turn network generalized by networkx into an adjacency matrix
    
    """
    if G == None:
        if nettype == 'ER':
            G = nx.random_graphs.erdos_renyi_graph(9, 0.3)
        elif nettype == 'WS':
            G = nx.random_graphs.watts_strogatz_graph(9, 3, 0.3)
        elif nettype == 'BA':
            G = nx.random_graphs.barabasi_albert_graph(9, 1)
   """network parameters should be changed here"""     

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
    """function to randomly turn interaction in a random network into an inhibition interaction
    rate is the probability of turning a normal interaction into an inhibition interaction
    the defult value is set to 0.72 which is the rate of inhibition interaction in yeast cell cycle network
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
     """Parameter setting part ends here"""
    A_1 = list(A)
    for i in range(len(A_1)):
        A_1[i] = list(A_1[i])    
    """unzip part ends here"""
    pos = mat_find_nonzero(A = A)
    n = len(pos)
    for i in range(n):
        r = random.random()
        if r < rate:
            A_1[pos[i][0]][pos[i][1]] = -1
            
    for i in range(len(A_1)):
        A_1[i] = tuple(A_1[i])
    """zip part"""
    A_2 = tuple(A_1)
    return A_2

def Basin_with_plot(A = None,threshold = None):
    """Function for calculating the attractor basion of a system with adjacency matrix A
    there are two parameters
    A is the adjacency matrix of origin system
    threshold is the updating threshold
    function return to a dictionary
    keys are atrractors
    and key's value is the basin of corresponding key 
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
     """Parameter setting part ends here"""    
    A_1 = list(A)
    for i in range(len(A_1)):
        A_1[i] = list(A_1[i])    
    """unzip part ends here"""
    
    if len(A_1) != len(threshold):                                       #A threshold must be with same dimension
        raise AssertionError("A threshold must be with same dimension")      
    else:
        dim = len(A_1)
        
    state_space = list(itertools.product([0,1],repeat = dim))          #initialize state space
    state_space_distribution = dict((('attractor',[]),('basin',[])))    #initialize state space distribution
    for i in range(len(state_space)):                                  
        state_space[i] = list(state_space[i])
    Source_index_list = []
    Target_index_list = []
    state_space_tuple = list(itertools.product([0,1],repeat = dim))
        
    flag = 0       #flag = 0 means need to get a new remaining state, flag = 1 means are continuously updating
    judge = 0      #judge=0 means there is a new attractor, judge = 1 means attractor has been collected
    x_tem = []     
    b = []    
    R_0 = [0] * len(state_space_tuple)
    R_mat = [R_0] * len(state_space_tuple)
    R = np.array(R_mat)
    
    while len(state_space) > 0:       #when there are still states not have been updated, continously update
        if flag == 0:                 #flag = 0 means updated to an attractor, need to get a new remaining state
            x0 = state_space[0]
            tem = []
            tem_index = []
        else:                         #flag = 1 means current state is not attractor, need to continuously updating
            x0 = x_tem
        
        Source_index = state_space_tuple.index(tuple(x0))
        Source_index_list.append(Source_index)
        tem.append(x0)                
        tem_index.append(Source_index)
        x_tem = bio_update(A,x0,threshold)   
        Target_index = state_space_tuple.index(tuple(x_tem))
        Target_index_list.append(Target_index)
        

        
        if x_tem in tem:                 #if x_tem in tem, means x_tem is an attractor
            pos = tem.index(x_tem)       #find x_tem's position
            if len(state_space_distribution['attractor']) > 0:        #if there are already other attractor, need to know is x_tem already there
                judge = 0                                             
                for i in range(len(state_space_distribution['attractor'])):       
                    if x_tem in state_space_distribution['attractor'][i] or x_tem in state_space_distribution['attractor']:
                        judge = 1
                        ind = i
                        
            if judge == 0:   #if x_tem is a new attractor
                state_space_distribution['attractor'].append(tem[pos:])   #attractor recording
                state_space_distribution['basin'].append(tem)   #attractor basin recording
                for i in range(len(tem)):      #remove all states in basin from remaining states
                    if tem[i] in state_space:
                        state_space.remove(tem[i])
                
                tem_index.append(Target_index)
                for r_0 in range(len(tem_index)-1):
                    for r_n in range(len(tem) - r_0):
                        R[tem_index[r_0]][tem_index[r_0 + r_n + 1]] = 1
                
                del tem          #clear tem
                del tem_index
                flag = 0         #clear flag
            else:    #else, x_tem is already there, just add new basin states to existing attractor basin
                if len(tem) == 1:                 
                    state_space_distribution['basin'][ind].append(tem)   
                else:       
                    for i in range(len(tem)):
                        state_space_distribution['basin'][ind].append(tem[i])
                        
                for i in range(len(tem)):  #remove all states in basin from remaining states
                    if tem[i] in state_space:
                        state_space.remove(tem[i])
                        
                tem_index.append(Target_index)
                for r_0 in range(len(tem_index)-1):
                    for r_n in range(len(tem) - r_0):
                        R[tem_index[r_0]][tem_index[r_0 + r_n + 1]] = 1
                        
                del tem #clear tem
                del tem_index
                del ind 
                flag = 0 #clear flag
        else:
            flag = 1 #x_tem is not attractor, keep on updating
        

    for i in range(len(state_space_distribution['basin'])):   #remove basin duplicated value
        b.append(list_duplicate_removal(state_space_distribution['basin'][i]))   
                
    state_space_distribution.update(basin=b)  
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
    return state_space_network, R


def trace(A = None):
    """
    function for get adjacency matrix A's trace
    """
    if A.any == None:
        A = np.array([[1,0,0],[0,1,0],[0,0,0]])
    
    tr = 0
        
    for i in range(len(A)):
        tr += A[i][i]
        
    return tr

def diag_array(A = np.array([[1,0,0],[0,1,0],[0,0,0]])):
    """
    function for get adjacency matrix A's diag array
    """
    diag = []
        
    for i in range(len(A)):
        diag.append(A[i][i])
        
    diag = tuple(diag)
        
    return diag
