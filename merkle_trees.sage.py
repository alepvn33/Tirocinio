

# This file was *autogenerated* from the file merkle_trees.sage
from sage.all_cmdline import *   # import sage library

_sage_const_10 = Integer(10); _sage_const_0 = Integer(0); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1)
import hashlib;
load('merkle_tree_utils.sage');

data = "Ciao sono Alessio e sto studiando i Merkle Tree";

num_blocks = _sage_const_10 ; #Number of data blocks or number of tree leaves
#Leaves in perfect binary trees must be a power of 2 number, so num_blocks is approximated upper to the next power of 2
i = _sage_const_0 ;
while _sage_const_2 **i < num_blocks: 
    i+=_sage_const_1 ;
num_blocks = _sage_const_2 **i;

tree = merkle_tree_generator(data,num_blocks)


print(tree)
tree[len(tree)-_sage_const_1 ][_sage_const_0 ].plot()

