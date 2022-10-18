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
def get_seed_path(I,starting_seed):
    
    path = []

    #Get number of leaves of the seeds binary tree
    n_leaves = 1
    while n_leaves < N:
        n_leaves = 2*n_leaves
    indexes = list(range(0,n_leaves)) #List of leaves indexes

    last_level = str(starting_seed)   #Root node
    
    while not len(indexes) == 1:

        sx, dx = hash_forking(last_level)   

        if I in indexes[:len(indexes)//2]: #I is on the left part
            path.append(dx)
            last_level = sx     #New parent node
            del indexes[len(indexes)//2:] #Eliminate indexes we don't need anymore
        else: #I is on the right part
            path.append(sx)
            last_level = dx     #New parent node
            del indexes[:len(indexes)//2] #Eliminate indexes we don't need anymore

    if indexes[0] == I:
        return path
    else:
        print("ERRORE NEL CALCOLO DEL SEED PATH")
        
################################################################## 

#Function that takes a list a of elements as input, generates a Merkle Tree from it and returns the call to the constructor
def MerkleTree(a):

    #https://stackoverflow.com/questions/18754180/create-multiple-instances-of-a-class
    mt = MerkleTools()

    while not log(len(a),2)==ceil(log(len(a),2)):
        set_random_seed()
        ### CONTROLLARE (numero a caso + hash perchè i commitment sono hash, quindi boh mi è sembrato logico)
        a.append(hash(str(randint(0,10^40)))) 

    d = log(len(a),2)

    for i in range(0,2^(d)-1): 
        mt.add_leaf(a[i]) #Not automatically hashed, the input a is intended to be already hashed

    mt.make_tree()

    return mt

################################################################## 

#Outputs the list of reconstructed seeds from the seed path obtained from get_seed_path (obviously, seeds with j in S are not included in the list)
def reconstruct_seeds(I,path_seed):

    #Get number of leaves of the seeds binary tree
    n_leaves = 1
    while n_leaves < N:
        n_leaves = 2*n_leaves
    indexes = list(range(0,n_leaves)) #List of leaves indexes

    seed = ['']*n_leaves #Final list of seeds we are looking for

    #Given that the seed path list has been filled from higher levels (near root) to lower levels (near leaves), i is basically the level of the tree that digest belongs to
    for i in range(len(path_seed)-1,-1,-1): #This loop goes from i=len(path_seed)-1 to i=0 (included)

        last_level = [path_seed[len(path_seed)-1-i]]

        for j in range(0,i): #As said before, i is the level of the tree where the digest is, so is also the number of hash forkings that we have to compute

            next_level = []

            for node in last_level:
                sx, dx = hash_forking(node) 
                next_level += [sx,dx]

            last_level = next_level

        #At this point, the list last_level contains all the leaves associated with the i-th element of path_seed
        #Now we have to decide where to put them on the final seed list

        if I in indexes[:len(indexes)//2]: #I is on the left part, so the seed list we have just computed is on the right
            considered_idx = indexes[len(indexes)//2:]
            seed[slice(min(considered_idx),max(considered_idx)+1)] = last_level
            del indexes[len(indexes)//2:] #Eliminate indexes we don't need anymore
        else: #I is on the right part, so the seed list we have just computed is on the left
            considered_idx = indexes[:len(indexes)//2]
            seed[slice(min(considered_idx),max(considered_idx)+1)] = last_level
            del indexes[:len(indexes)//2] #Eliminate indexes we don't need anymore

    if seed[I] == '':
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
        
        
        self.aux = []                               # list of N aux lists of strings
        self.aux_hash = [[] for x in range(N)]      # hashed version of aux (leaves must be obtained by hashing the aux values)
        self.c = []                               # list of N commitments
        self.c_hash = []                          # hashed version of c
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
    def Response(self,ch,I):

        print("I: "+str(I))#########################
        rsp_I = self.P2(ch,int(self.seed[I],16),I)

        mt = self.T_aux[I]  
        path_aux_I = mt.get_proof(int(ch))
        mt = self.T_c
        path_c_I = mt.get_proof(int(I))
        path_seed = get_seed_path(I,self.starting_seed)

        #To remove the I-th seed, we insert a 0 in the I-th position of the list. The verifier knows he doesn't have to check that position
        seed_not_I = self.seed
        seed_not_I[I] = 0

        return rsp_I, path_aux_I, path_c_I, path_seed

##################################################################  

#VERIFIER
class Verifier:

    def __init__(self,hash_leaf):

        self.hash_leaf = hash_leaf # True to hash the input leaves, False otherwise

        self.prover = Prover()

        self.T_aux_bar = []                                     # list of calls to the T_aux_bar merkle trees constructors
        self.aux_bar_not_I = ['']*N                             # list of the N aux_bar reconstructed by every seed (each entry is a list of length q)
        self.aux_bar_not_I_hash = [[] for x in range(N)]        # hashed version of aux (leaves must be obtained by hashing the aux values)
        self.root_aux_bar = ['']*N                              # list of the roots of all the T_aux_bar merkle trees

    #Challenge (V1)
    def V1(self):

        set_random_seed()
        self.z = Fq.random_element()
        ch = self.z

        return ch
    
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

        set_random_seed()
        self.I = int(Integers(N).random_element())
        ch = self.V1()

        return self.I, ch

    #II. VERIFICATION
    def Verification(self,Htr_unsys,s,h,root_c,rsp_I,path_aux_I,path_c_I,path_seed):

        seed_bar = reconstruct_seeds(int(self.I),path_seed) #According to the paper, seeds should be taken one by one for every j in range(0,N), but is easier to get the whole seeds with j not in S at the same time
        for j in range(0,N):
            if j == self.I:
                tau_y_I = apply_rest_monomial(rsp_I[2],rsp_I[3],rsp_I[4])
                t_I = tau_y_I[0:r]+tau_y_I[r:n]*Htr_unsys - self.z*s
                c_I = hash(str(rsp_I[0])+str(rsp_I[2])+str(rsp_I[3])+str(t_I))

                if self.hash_leaf:
                    c_I_hash = hash(c_I)
                else:
                    c_I_hash = c_I

                mt = MerkleTools() 
                b = int(mt.validate_proof(path_c_I, c_I_hash, root_c)) #validate_proof() is equal to get_merkle_root() (the function that they call "ReconstructRoot" in the paper) + check with target root

                c_z_I = hash(str(rsp_I[1])+str(rsp_I[4]))
                if self.hash_leaf: 
                    c_z_I_hash = hash(c_z_I)
                else:
                    c_z_I_hash = c_z_I

                self.root_aux_bar[j] = mt.get_root_from_path(path_aux_I, c_z_I_hash)
                self.T_aux_bar.append('') #T_aux_bar must have j entries, even if this one won't be used
            else:
                self.aux_bar_not_I[j] = Helper.H(int(seed_bar[j],16))

                if self.hash_leaf:
                    for i in range(0,len(self.aux_bar_not_I[j])):
                        self.aux_bar_not_I_hash[j].append(hash(self.aux_bar_not_I[j][i]))
                else:
                    self.aux_bar_not_I_hash[j] = self.aux_bar_not_I[j]

                self.T_aux_bar.append(MerkleTree(self.aux_bar_not_I_hash[j]))
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

def MPC_id(e,Htr_unsys,s):
    
    prover = Prover(True)       #True for abilitating hashing of input leaves
    verifier = Verifier(True)   #True for abilitating hashing of input leaves

    #I. COMMITMENT (Prover)
    h, root_c = Prover.Commitment(prover,Htr_unsys,e)

    #II. CHALLENGE (Verifier)
    I, ch = Verifier.Challenge(verifier)

    #III. RESPONSE (Prover)
    rsp_I, path_aux_I, path_c_I, path_seed = Prover.Response(prover,ch,I)

    #IV. VERIFICATION (Verifier)
    return Verifier.Verification(verifier,Htr_unsys,s,h,root_c,rsp_I,path_aux_I,path_c_I,path_seed)