#Hashes a couple of the last level nodes into a next level node
def joiner(tree):
    last_level = len(tree)-1; #Last filled level
    tree.append([]);    #Beginning of next level

    for i in range(0,len(tree[last_level]),2):
        label_0 = tree[last_level][i].sort_key()[0][1];
        label_1 = tree[last_level][i+1].sort_key()[0][1];
        hash_input = str(label_0)+str(label_1);
        hash_output = hashlib.sha256();
        hash_output.update(hash_input.encode('utf-8'));
        hash_output = hash_output.hexdigest();

        tree[last_level+1].append(LabelledOrderedTree([tree[last_level][i],tree[last_level][i+1]], label = hash_output));
    
    return tree

########################################################

def merkle_tree_generator(data,num_blocks,debug):
    data_bin = ''.join(format(ord(x), 'b') for x in str(data)); #Input data in binary representation
    blocks = []; leaves = [];
    block_len = int(len(data_bin)/num_blocks)+1; #Bit lenght of each block

    for i in range(0,num_blocks):
        #Debug printing
        if debug:
            print("Data block no."+str(i).zfill(2)+": "+str(data_bin[i*block_len:(i+1)*block_len]))   
        blocks.append(data_bin[i*block_len:(i+1)*block_len]);
        leaves.append(hashlib.sha256());
        leaves[i].update(blocks[i].encode('utf-8'));
        leaves[i] = leaves[i].hexdigest();

    #Building the lower level of the tree (leaves)
    tree = [[]];
    for x in range(0,num_blocks):
        tree[0].append(LabelledOrderedTree([], label = leaves[x]));

    #Building the upper levels until the root is reached
    while len(tree[len(tree)-1]) != 1:
        tree = joiner(tree);

    #Debug printing
    if debug: 
        print("\n****************\nBlock length: "+str(block_len)+"\nBinary data lenght: "+str(len(data_bin))+"\nBlocks number: "+str(num_blocks)+"\n****************\n\nTree:")
        print(tree)
    
    return tree