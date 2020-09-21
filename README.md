# DReSS-package
Defult Network is yeast cell cycle regulatory network whose adjacency matrix is 
A = ((-1,0,-1,-1,0,0,0,0,0),
     (0,0,-1,-1,0,0,-1,1,0),
     (0,-1,0,0,0,-1,0,0,0),
     (0,-1,0,0,0,-1,0,0,0),
     (0,-1,0,0,-1,-1,0,0,1),
     (0,0,-1,-1,1,0,0,0,0),
     (0,0,0,0,0,-1,0,0,0),
     (0,0,0,0,0,1,0,0,0),
     (0,0,1,1,0,0,1,-1,-1))
     
## main(A = None, mode = (1,0), threshold = None, print_switch = 'on')
Main Function for calculating DReSS
    There are four parameters
    A is adjacency matrix of origin network
    mode can be chosen from follows:
                     (1,0) --- perturbation that delete existing interactions
                     (0,1) --- perturbation that add activation interactions
                     (0,-1) --- perturbation that add inhibition interactions
    print_switch is a parameter choose to display current progress or not

## diag_main(A = None, mode = (1,0), threshold = None, print_switch = 'on')
Main Function for calculating diagDReSS
    There are four parameters
    A is adjacency matrix of origin network
    mode can be chosen from follows:
                     (1,0) --- perturbation that delete existing interactions
                     (0,1) --- perturbation that add activation interactions
                     (0,-1) --- perturbation that add inhibition interactions
    print_switch is a parameter choose to display current progress or not
    
## Bonferroni(A = None, B = None):
Function for Bonferroni Test with two parameters
    A is the test set with position information
    B is the control set with only values from randomized controlled trial
    result return to a dictionary with Position, p-value, significance
    and write it to an csv file named 'result.csv'
    
## random_times(times = None, nettype = 'ER',mode = (1,0))
Function for DReSS randomized controlled trial
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
    
## random_times_diag(times = None, nettype = 'ER',mode = (1,0)):
Function for diagDReSS randomized controlled trial
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

## result_read(result_mode = 0):
read function DReSS result and turn it into format that can be read by function Bonferroni
    result_mode can be chosen from follows:
        'r0' for reading real_network deleting interaction results
        'r1' for reading real_network adding a actication interaction results
        'r-1' for reading real_network adding a inhibition interaction results
        '0' for reading randomized trials deleting interaction results
        '1' for reading randomized trials adding a actication interaction results
        '-1' for reading randomized trials adding a inhibition interaction results
        
## result_read_diag(result_mode = 0)
read function diagDReSS result and turn it into format that can be read by function Bonferroni
    result_mode can be chosen from follows:
        'r0' for reading real_network deleting interaction results
        'r1' for reading real_network adding a actication interaction results
        'r-1' for reading real_network adding a inhibition interaction results
        '0' for reading randomized trials deleting interaction results
        '1' for reading randomized trials adding a actication interaction results
        '-1' for reading randomized trials adding a inhibition interaction results
        
## bio_update(A = None,x0 = None,threshold = None)
state updating function for one step update
    There are 3 parameters
    A is the adjacency matrix for origin system
    x0 is the current state
    threshold is the updating threshold
    function returns to x0's next time step state
    
    Notice: A x0 threshold must with same dimension
        Otherwise function return to error
        
## A_replace(A = None,position = None,interaction = None)
Function for replace one positon's value with interaction value
    A is the adjacency matrix
    position is the replace target positon 
    interaction is the new value to replace the origin value
    
## graph_from_A(A = None,draw = 1, Nodes = None)
    Function to turn agjacency matrix A into a graph which can be read by networkx package
    draw can be chosen from follows:
        1 plot network
        0 do not plot network
    Nodes can provide names for each node
    
## Adjacency_marix(a = None)
Function to get adjacency matrix from a adjacency list
    input a is a 2 row tuple
    first row is source
    second row is target
    function return to a's corresponding adjacency matrix
    
## Reachability_matrix(A = None)
Function to get modified Reachability matrix R
    input A is the adjacency matrix of state space
    function returns to the modified reachability matrix of state space with adjacency matrix A
    
## DReSS(A = np.array([[1,0,0],[0,1,0],[0,0,0]]), B = np.array([[1,0,1],[0,1,0],[1,0,1]]))
Function for calculating DReSS value
    A and B are two modified reachablity matrix
    
## haming_diag(A = (1,0,1),B = (0,1,1))
    function for calculating diagDReSS
    A and B are diag array of two modified reachablity matrix
    
## mat_find_zero(A = None)
Function for find all position of zero values in an adjacency matrix A 

## mat_find_nonzero(A = None)
Function for find all position of non-zero values in an adjacency matrix A 

## graph2matrix(G = None,nettype = 'ER')
function to turn network generalized by networkx into an adjacency matrix

## rand_negative(A = None, rate = 0.72)
function to randomly turn interaction in a random network into an inhibition interaction
    rate is the probability of turning a normal interaction into an inhibition interaction
    the defult value is set to 0.72 which is the rate of inhibition interaction in yeast cell cycle network
    
## Basin_with_plot(A = None,threshold = None)
Function for calculating the attractor basion of a system with adjacency matrix A
    there are two parameters
    A is the adjacency matrix of origin system
    threshold is the updating threshold
    function return to a dictionary
    keys are atrractors
    and key's value is the basin of corresponding key 
    
## trace(A = None):
function for get adjacency matrix A's trace

## diag_array(A = np.array([[1,0,0],[0,1,0],[0,0,0]]))
function for get adjacency matrix A's diag array
