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
function main is to purterbate all possible position of one edge in the network, and return to each purterbation's DReSS value.

There is three mode, (1,0) represents to delete an existing interaction, (0,1) represents to add an activation interaction, (0,-1) represents to add an inhibition interaction.
