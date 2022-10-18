#Docs available at https://github.com/Tierion/pymerkletools
import merkletools;
import hashlib;
load('merkletools.sage')
load('merkle_trees_utils.sage')

mt = MerkleTools()

data = ["1","2","3","4","5"]

#Create the Merkle Tree
MerkleTree(data)

if mt.is_ready:

    #Get tree root
    merkleRoot =  mt.get_merkle_root()

    #Get authentication path
    auth_path = mt.get_proof(0)

    #Validate 
    is_valid = mt.validate_proof(proof, targetHash, merkleRoot)

else:
    print("Tree not created")



