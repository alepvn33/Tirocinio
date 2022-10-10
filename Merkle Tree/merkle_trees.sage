#Docs available at https://github.com/Tierion/pymerkletools
import merkletools;
import hashlib;
load('merkletools.sage')

mt = MerkleTools()

data = "Ciao sono Alessio e sto studiando i Merkle Tree"
num_blocks = 10 #Number of data blocks or number of tree leaves
#Leaves in perfect binary trees must be a power of 2 number, so num_blocks is approximated upper to the next power of 2
i = 0
while 2^i < num_blocks: 
    i+=1
num_blocks = 2^i

data_bin = ''.join(format(ord(x), 'b') for x in str(data)) #Input data in binary representation
blocks = []
block_len = int(len(data_bin)/num_blocks)+1 #Bit lenght of each block

for i in range(0,num_blocks): 
    blocks.append(data_bin[i*block_len:(i+1)*block_len])
    mt.add_leaf(blocks[i], True)

mt.make_tree()

if mt.is_ready:
    merkleRoot =  mt.get_merkle_root()
    for i in range(0,num_blocks):
        targetHash = hashlib.sha256()
        targetHash.update(blocks[i].encode('utf-8'))
        targetHash = targetHash.hexdigest()
        proof = mt.get_proof(i)
        is_valid = mt.validate_proof(proof, targetHash, merkleRoot)
        if is_valid:    print("Block no."+str(i+1).zfill(2)+": VALID")
        else:   print("Block no."+str(i+1).zfill(2)+": NOT VALID")
else: print("Tree not created")



