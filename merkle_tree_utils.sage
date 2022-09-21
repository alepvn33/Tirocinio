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

def merkle_tree_generator(data,num_blocks):
    data_bin = ''.join(format(ord(x), 'b') for x in data);
    blocks = []; leaves = [];
    block_len = int(len(data_bin)/num_blocks)+1; #Bit lenght of each block

    #Padding the data string
    pad_bits_num = len(data_bin)-block_len*num_blocks;
    ''.join('0' for x in range(0,pad_bits_num));

    for i in range(0,num_blocks):
        blocks.append(data_bin[i:i*block_len]);
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
    
    return tree