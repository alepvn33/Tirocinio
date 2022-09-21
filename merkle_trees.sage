import hashlib;
load('merkle_tree_utils.sage');

debug = True
data = "Ciao sono Alessio e sto studiando i Merkle Tree";
num_blocks = 10; #Number of data blocks or number of tree leaves
#Leaves in perfect binary trees must be a power of 2 number, so num_blocks is approximated upper to the next power of 2
i = 0;
while 2^i < num_blocks: 
    i+=1;
num_blocks = 2^i;


#Debug printing
if debug:   print("\n*****DEBUG MODE*****\n\nData: "+str(data)+'\n')

tree = merkle_tree_generator(data,num_blocks,debug)
