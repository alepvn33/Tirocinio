load("merkletools.sage")

mt = MerkleTools()

#Hash input
def hash(input):

    dig = hashlib.sha256()
    dig.update(input.encode('utf-8'))
    dig = dig.hexdigest()

    return dig

#Function that takes a list a of elements as input and generates a Merkle Tree from it
def MerkleTree(a):

    while not log(len(a),2)==ceil(log(len(a),2)):
        set_random_seed()
        ### CONTROLLARE (numero a caso + hash perchè i commitment sono hash, quindi boh mi è sembrato logico)
        a.append(hash(str(randint(0,10^40)))) 

    d = log(len(a),2)

    for i in range(0,2^(d)-1): 
        mt.add_leaf(a[i], True) #Automatically hashed ("True")

    mt.make_tree()