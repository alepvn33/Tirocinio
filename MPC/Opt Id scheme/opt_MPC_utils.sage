#MPC functions

from random import seed


class Bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

#Class that represents a node of the seed tree
class Node:

    def __init__(self,val='',lev=-1,ind=-1,dir=''):
        self.val = val  #Value of the node
        self.lev = lev  #Level within the tree of the node
        self.ind = ind  #Index of the node in that level
        self.dir = dir  #Direction of the node (left or right)

    def level(self):
        return self.lev

    def value(self):
        return self.val
    
    def index(self):
        return self.ind

    def direction(self):
        return self.dir
    

#Key generation
def key_gen():

    set_random_seed()
    Htr_unsys = random_matrix(Fq,n-r,r) #public matrix (only non-systematic portion)
    e = FWV()

    #compute syndrome
    s = e[0:r] + e[r:n]*Htr_unsys

    return e,Htr_unsys,s

##################################################################

#Generate a random vector with weight w and lenght n
def FWV():

    a = vector(Fq,w)

    for i in range(0,w):
        a[i] = Fq_star.random_element()

    b = zero_vector(n) #List of zeros of lenght n
    P = Permutations(range(0,n))
    rnd_supp = P.random_element()

    for i in range(0,w):
        b[rnd_supp[i]] = a[i]

    return b

##################################################################

#Find an isometry function tau such that e = tau(e_tilde)
def find_isometry(e,e_tilde):

    supp_e = []
    supp_e_tilde = []

    for i in range(0,n):
        if e[i]!=0:
            supp_e.append(i)
        if e_tilde[i]!=0:
            supp_e_tilde.append(i)

    P_supp = Permutations(supp_e_tilde)
    tau_perm_supp = P_supp.random_element()

    zeros_e_tilde = [x for x in list(range(0,n)) if x not in supp_e_tilde] #Complement of e_tilde support
    P_zeros = Permutations(zeros_e_tilde)
    tau_perm_zeros = P_zeros.random_element()

    #Joining the permutations of e_tilde support and e_tilde support complementary, according to e support
    tau_perm = [] #Total permutation
    j = 0
    k = 0

    for i in range(0,n):
        if i in supp_e:
            tau_perm.append(tau_perm_supp[j])
            j = j + 1
        else:
            tau_perm.append(tau_perm_zeros[k])
            k = k + 1

    tau_values = vector(Fq,n)

    for i in range(0,n):
        #Find value vector through inverse calculation
        if e_tilde[tau_perm[i]] == 0:
            tau_values[tau_perm[i]] = Fq_star.random_element()
        else:
            tau_values[tau_perm[i]] = Fq(e[i]/e_tilde[tau_perm[i]])
    
    return tau_perm, tau_values

##################################################################   

#Apply monomial transformation (isometry)
def apply_rest_monomial(tau_perm,tau_values,a):

    b = vector(Fq,n)

    for i in range(0,n):
        p = tau_perm[i]
        b[i] = tau_values[p]*a[p]

    return b

##################################################################   

#Hash input
def hash(input):

    dig = hashlib.sha256()
    dig.update(input.encode('utf-8'))
    dig = dig.hexdigest()

    return dig

##################################################################   

#Compute the hash of an input node and separates the result in two parts
def hash_forking(input):

    dig = hash(input) 
    sx = dig[:len(dig)//2]
    dx = dig[len(dig)//2:]

    return sx, dx

#Generates the N seeds as leaves of a binary tree obtained by a starting seed
def SeedTree(starting_seed):

    last_level = [str(starting_seed)]

    while len(last_level) < N:

        next_level = []

        for i in last_level:
            sx, dx = hash_forking(i)
            next_level.append(sx) # Sibling node (SX)
            next_level.append(dx) # Sibling node (DX)
        
        last_level = next_level

    return last_level[0:N] # Return the N required seeds


################################################################## 

#Returns the seed path of I-th seed which is the I-th leaf of the binary tree obtained from starting_seed
def get_seed_path(S,starting_seed):
    
    path = []

    #Get number of leaves of the seeds binary tree
    n_leaves = 1
    while n_leaves < N:
        n_leaves = 2*n_leaves
    indexes = list(range(0,n_leaves)) #List of leaves indexes
    
    #First, we get the path of just the first seed whose index is in S
    last_level = str(starting_seed)   #Root node
    
    levels = log(n_leaves,2)

    lv_counter = levels
    last_index = 0
    while not len(indexes) == 1:
        lv_counter -= 1
        sx, dx = hash_forking(last_level)   
        if S[0] in indexes[:len(indexes)//2]: #S[0] is on the left part
            path.append(Node(dx,lv_counter,2*last_index+1,'right'))
            last_index = 2*last_index
            last_level = sx     #New parent node
            del indexes[len(indexes)//2:] #Eliminate indexes we don't need anymore
        else: #S[0] is on the right part
            path.append(Node(sx,lv_counter,2*last_index,'left'))
            last_index = 2*last_index + 1
            last_level = dx     #New parent node
            del indexes[:len(indexes)//2] #Eliminate indexes we don't need anymore

    if not indexes[0] == S[0]:
        print("ERRORE NEL CALCOLO DEL SEED PATH")
        return
  
    for i in S: 
        if not i == S[0]:
            indexes = [list(range(0,n_leaves/(2^x))) for x in range(0,levels)] #List of node index lists for every level
            for j in range(levels-1,-1,-1): #j is the level index (from levels-1 to 0 (leaves level) included)
                #i is on the left part
                if i in indexes[0][:len(indexes[0])//2]:  
                    #Cut the indexes in half at every level
                    for ind in indexes: #Eliminate indexes we don't need anymore
                        if len(ind) == 1:
                            del ind[0]
                        else:
                            del ind[len(ind)//2:] 
                    for x in path:
                        if x.level()==j and x.index()==indexes[j][0]:
                            #Remove the node and add its siblings (only if not included in S)
                            path.remove(x)
                            #Calculate sibling nodes
                            sx, dx = hash_forking(x.value())
                            if j-1 == 0 and not indexes[j-1][0] in S or j-1 > 0:
                                #Left sibling node
                                path.append(Node(sx,j-1,indexes[j-1][0],'left'))
                            if j-1 == 0 and not indexes[j-1][1] in S or j-1 > 0:
                                #Rigth sibling node
                                path.append(Node(dx,j-1,indexes[j-1][1],'right'))
                #i is on the right part
                else: 
                    #Cut the indexes in half at every level
                    for ind in indexes: #Eliminate indexes we don't need anymore
                        if len(ind) == 1:
                            del ind[0]
                        else:
                            del ind[:len(ind)//2] 
                    for x in path:
                        if x.level()==j and x.index()==indexes[j][0]:
                            #Remove the node and add its siblings (only if not included in S)
                            path.remove(x)
                            #Calculate sibling nodes
                            sx, dx = hash_forking(x.value())
                            if j-1 == 0 and not indexes[j-1][0] in S or j-1 > 0:
                                #Left sibling node
                                path.append(Node(sx,j-1,indexes[j-1][0],'left'))
                            if j-1 == 0 and not indexes[j-1][1] in S or j-1 > 0:
                                #Rigth sibling node
                                path.append(Node(dx,j-1,indexes[j-1][1],'right'))

    return path        
            
        
################################################################## 

#Function that takes a list a of elements as input, generates a Merkle Tree from it and returns the call to the constructor
def MerkleTree(input):

    a = input

    #https://stackoverflow.com/questions/18754180/create-multiple-instances-of-a-class
    mt = MerkleTools()

    counter = 0
    while not log(len(a),2)==ceil(log(len(a),2)):
        counter += 1
        #We add 'custom' leaves to get the right number of leaves (power of 2). This is done by hashing the last leaf (we can't add random leaves because this would
        #cause the reconstructed tree (created by Verifier) to be different to the first one (created by Prover) in the custom leaves if the set seed isn't the same)
        a.append(hash(a[-1])) 

    d = log(len(a),2)

    for i in range(0,2^(d)): 
        mt.add_leaf(a[i]) #Not automatically hashed, the input a is intended to be already hashed

    mt.make_tree()

    #Removing custom leaves now that we no longer need them
    for i in range(0,counter):
        a.pop(-1)

    return mt

################################################################## 

#Outputs the list of reconstructed seeds from the seed path obtained from get_seed_path (obviously, seeds with j in S are not included in the list)
def reconstruct_seeds(S,path_seed):

    #Get number of leaves of the seeds binary tree
    n_leaves = 1
    while n_leaves < N:
        n_leaves = 2*n_leaves

    seed = ['']*n_leaves #Final list of seeds we are looking for

    for node in path_seed:

        #First, we have to reconstruct the leaves associated to each node in path_seed
        last_level = [node]

        for i in range(0,node.level()): #As said before, i is the level of the tree where the digest is, so is also the number of hash forkings that we have to compute

            next_level = []
            next_level_indexes = []

            for x in last_level:
                sx, dx = hash_forking(x.value()) 
                next_level += [Node(sx,x.index()-1,2*x.index(),'left'),Node(dx,x.index()-1,2*x.index()+1,'right')]
                next_level_indexes += [2*x.index(), 2*x.index()+1]

            last_level = next_level
        
        if node.level() == 0:
            next_level_indexes = [node.index()]

        #After recontructing the leaves, we have to place them in the final seed list
        counter = 0
        for i in next_level_indexes:        
            seed[i] = last_level[counter].value()
            counter += 1
    
    test = 1
    for i in S:
        if not seed[i] == '':
            test = 0
    if test and len(seed)==n_leaves:
        return seed
    else:
        print("ERRORE NEL RICALCOLO DEI SEED")

    
##################################################################  
#Here the actor (Helper, Prover, Verifier) classes description begins. 
#If the comments above the method definitions are in:
# - lower case letters: the method is one of the steps of the first MPC-in-the-head scheme (the one with the Helper)
# - upper case letters: the method is one of the steps of the second MPC-in-the-head scheme (the one without the Helper), which use the "lower case" methods

#HELPER
class Helper:

    #Setup (H)
    def H(seed):

        set_random_seed(seed)
        u_helper = random_vector(Fq,n)
        e_tilde_helper = FWV()

        r_v = matrix(GF(2),q,_lambda)
        c_v = [""]*len(Fq)

        set_random_seed(seed)
        for v in Fq:
            r_v[v] = random_vector(GF(2),_lambda)
            c_v[v] = hash(str(r_v[v])+str(u_helper+v*e_tilde_helper))

        aux = c_v

        return aux     

##################################################################  

#PROVER
class Prover:

    def __init__(self,hash_leaf=False):

        self.hash_leaf = hash_leaf # True to hash the input leaves, False otherwise

        self.u = []
        self.e_tilde = []
        self.tau_perm = []
        self.tau_values = []
        self.r_vec = []
        
        self.rsp_S = []                             # list of s responses
        self.path_aux_S = []                        # list of s aux for the selected rounds
        self.path_c_S = []                          # list of s paths of the commitments (c) for the selected rounds
        self.aux = []                               # list of N aux lists of strings
        self.aux_hash = [[] for x in range(N)]      # hashed version of aux (leaves must be obtained by hashing the aux values)
        self.c = []                                 # list of N commitments
        self.c_hash = []                            # hashed version of c
        self.root_aux = [""]*N                      # list of the N roots of each aux tree
        self.T_aux = []                             # list of calls to the T_aux merkle trees constructors

    

    #Commitment (P1)
    def P1(self,Htr_unsys,e,seed):

        set_random_seed(seed)
        self.u.append(random_vector(Fq,n))
        self.e_tilde.append(FWV())
        tau_p, tau_v = find_isometry(e,self.e_tilde[-1])
        self.tau_perm.append(tau_p)
        self.tau_values.append(tau_v)
        self.r_vec.append(random_vector(GF(2),_lambda))
        tau_u = apply_rest_monomial(self.tau_perm[-1],self.tau_values[-1],self.u[-1])
        
        c = hash(str(self.r_vec[-1])+str(self.tau_perm[-1])+str(self.tau_values[-1])+str(tau_u[0:r]+tau_u[r:n]*Htr_unsys))
        
        return c

    #Response (P2)
    def P2(self,ch,seed,I):

        set_random_seed(seed)
        r_v_prover = matrix(GF(2),q,_lambda)

        for v in Fq:
            r_v_prover[v] = random_vector(GF(2),_lambda)

        r_z = r_v_prover[ch]
        y = self.u[I] + ch*self.e_tilde[I]
        rsp = [self.r_vec[I],r_z,self.tau_perm[I],self.tau_values[I],y]

        return rsp

    #I. COMMITMENT
    def Commitment(self,Htr_unsys,e):

        set_random_seed()
        self.starting_seed = initial_seed()         # Random starting seed, root of the seed tree
        self.seed = SeedTree(self.starting_seed)    # list of N seeds obtained from the seed tree

        for i in range(0,N):
            self.aux.append(Helper.H(int(self.seed[i],16)))
            if self.hash_leaf:
                for j in range(0,len(self.aux[i])):
                    self.aux_hash[i].append(hash(self.aux[i][j]))
            else:
                self.aux_hash[i] = self.aux[i] #We don't do the hash of the leaves
            mt = MerkleTree(self.aux_hash[i])            
            self.T_aux.append(mt)
            self.root_aux[i] = mt.get_merkle_root()
            self.c.append(self.P1(Htr_unsys,e,int(self.seed[i],16)))

        self.h = hash(''.join(self.root_aux)) # Total hash of aux trees roots
        if self.hash_leaf:
            for i in range(0,len(self.c)):
                self.c_hash.append(hash(self.c[i]))
        else:
            self.c_hash = self.c

        self.T_c = MerkleTree(self.c_hash)
        mt = self.T_c
        self.root_c = mt.get_merkle_root()
    
        return self.h, self.root_c

    #III RESPONSE
    def Response(self,z_j,S):
        
        for i in S:
            idx = S.index(i) #Index of i in the S list
            self.rsp_S.append(self.P2(z_j[idx],int(self.seed[i],16),i))

            mt = self.T_aux[i]  
            self.path_aux_S.append(mt.get_proof(int(z_j[idx])))
            mt = self.T_c
            self.path_c_S.append(mt.get_proof(int(i)))
     
        path_seed = get_seed_path(S,self.starting_seed)

        return self.rsp_S, self.path_aux_S, self.path_c_S, path_seed

##################################################################  

#VERIFIER
class Verifier:

    def __init__(self,hash_leaf,s_num):

        self.hash_leaf = hash_leaf  # True to hash the input leaves, False otherwise
        self.s_num = s_num                  # Number of rounds to select

        self.prover = Prover()

        self.S = []                                             # list of s_num selected rounds (from 0 to N-1)
        self.z = []                                             # list of s_num selected aux (from 0 to q-1)
        self.T_aux_bar = []                                     # list of calls to the T_aux_bar merkle trees constructors
        self.aux_bar_not_S = ['']*N                             # list of the N aux_bar reconstructed by every seed not included in S (each entry is a list of length q)
        self.aux_bar_not_S_hash = [[] for x in range(N)]        # hashed version of aux_bar_not_S (leaves must be obtained by hashing the aux values)
        self.root_aux_bar = ['']*N                              # list of the roots of all the T_aux_bar merkle trees

    #Challenge (V1)
    def V1(self):

        set_random_seed()
        z = Fq.random_element()

        return z
    
    #Verification (V2)
    def V2(self,Htr_unsys,s,aux,c,rsp):

        tau_y = apply_rest_monomial(rsp[2],rsp[3],rsp[4])
        t = tau_y[0:r]+tau_y[r:n]*Htr_unsys - self.z*s

        c_verifier = hash(str(rsp[0])+str(rsp[2])+str(rsp[3])+str(t))
        
        #The 2nd and 3rd conditions check if tau is an isometry
        if c_verifier == c and numpy.unique(rsp[2]).size == n and 0 not in rsp[3]:

            c_z_verifier = hash(str(rsp[1])+str(rsp[4]))
            c_z = aux[self.z]

            if c_z_verifier == c_z:
                return 1
            else:
                return 0
        else:
            return 0
    
    #II. CHALLENGE
    def Challenge(self):
        
        set = Set(Integers(N))

        for i in range(0,self.s_num):
            set_random_seed()
            self.S.append(int(set.random_element()))
            set = set.difference([self.S[i]])
            self.z.append(self.V1())

        return self.S, self.z

    #II. VERIFICATION
    def Verification(self,Htr_unsys,s,h,root_c,rsp_S,path_aux_S,path_c_S,path_seed):

        seed_bar = reconstruct_seeds(self.S,path_seed) #According to the paper, seeds should be taken one by one for every j in range(0,N), but is easier to get the whole seeds with j not in S at the same time
        for j in range(0,N):
            if j in self.S:
                idx = self.S.index(j) #Index of j in the self.S list
                tau_y_S = apply_rest_monomial(rsp_S[idx][2],rsp_S[idx][3],rsp_S[idx][4])
                t_S = tau_y_S[0:r]+tau_y_S[r:n]*Htr_unsys - self.z[idx]*s
                c_S = hash(str(rsp_S[idx][0])+str(rsp_S[idx][2])+str(rsp_S[idx][3])+str(t_S))

                if self.hash_leaf:
                    c_S_hash = hash(c_S)
                else:
                    c_S_hash = c_S

                mt = MerkleTools()
                b = int(mt.validate_proof(path_c_S[idx], c_S_hash, root_c)) #validate_proof() is equal to get_merkle_root() (the function that they call "ReconstructRoot" in the paper) + check with target root

                c_z_S = hash(str(rsp_S[idx][1])+str(rsp_S[idx][4]))
                if self.hash_leaf: 
                    c_z_S_hash = hash(c_z_S)
                else:
                    c_z_S_hash = c_z_S

                self.root_aux_bar[j] = mt.get_root_from_path(path_aux_S[idx], c_z_S_hash)
                self.T_aux_bar.append('') #T_aux_bar must have j entries, even if this one won't be used
            else:
                self.aux_bar_not_S[j] = Helper.H(int(seed_bar[j],16))

                if self.hash_leaf:
                    for i in range(0,len(self.aux_bar_not_S[j])):
                        self.aux_bar_not_S_hash[j].append(hash(self.aux_bar_not_S[j][i]))
                else:
                    self.aux_bar_not_S_hash[j] = self.aux_bar_not_S[j]

                self.T_aux_bar.append(MerkleTree(self.aux_bar_not_S_hash[j]))
                mt = self.T_aux_bar[j]
                self.root_aux_bar[j] = mt.get_merkle_root()

        h_bar = hash(''.join(self.root_aux_bar))
        if h_bar == h:
            b_prime = 1
            print(Bcolors.OKCYAN+"\nb_prime = 1"+Bcolors.ENDC)
        else:
            b_prime = 0
            print(Bcolors.WARNING+"\nb_prime = 0"+Bcolors.ENDC)

        if b == 1:
            print(Bcolors.OKCYAN+"\nb = 1"+Bcolors.ENDC)
        else:
            print(Bcolors.WARNING+"\nb = 0"+Bcolors.ENDC)
        
        return b and b_prime

##################################################################  

def MPC_id(e,Htr_unsys,s,s_num):
    
    prover = Prover(True)       #True for abilitating hashing of input leaves
    verifier = Verifier(True,s_num)   #True for abilitating hashing of input leaves

    #I. COMMITMENT (Prover)
    h, root_c = Prover.Commitment(prover,Htr_unsys,e)

    #II. CHALLENGE (Verifier)
    S, z_j = Verifier.Challenge(verifier)

    #III. RESPONSE (Prover)
    rsp_S, path_aux_S, path_c_S, path_seed = Prover.Response(prover,z_j,S)

    #IV. VERIFICATION (Verifier)
    return Verifier.Verification(verifier,Htr_unsys,s,h,root_c,rsp_S,path_aux_S,path_c_S,path_seed)