#MPC functions

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

#SIGNER
class Signer:

    def __init__(self):
        self.u = []
        self.e_tilde = []
        self.tau_perm = []
        self.tau_values = []
        self.r_vec = []

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

        self.seed = vector(ZZ,N)
        aux = [""]*N
        c = [""]*N
        for i in range(0,N):
            set_random_seed()
            self.seed[i] = initial_seed()
            aux[i] = Helper.H(self.seed[i])
            c[i] = self.P1(Htr_unsys,e,self.seed[i])
    
        return aux, c

    #III RESPONSE
    def Response(self,ch,I):

        rsp_I = self.P2(ch,self.seed[I],I)

        #To remove the I-th seed, we insert a 0 in the I-th position of the list. The verifier knows he doesn't have to check that position
        seed_not_I = self.seed
        seed_not_I[I] = 0

        return rsp_I, seed_not_I

    #II. CHALLENGE
    def Challenge(self,mex,aux,c):

        set_random_seed(Integer(hash(str(mex)+str(aux)+str(c)),16))
        self.I = Integers(N).random_element()
        self.z = Fq.random_element()
        ch = self.z

        return self.I, ch

    #Main execution of the Signer
    def main(self,Htr_unsys,e,mex):

        #I. COMMITMENT
        aux, c = self.Commitment(Htr_unsys,e)

        #II. CHALLENGE
        I, ch = self.Challenge(mex,aux,c)

        #III. RESPONSE
        rsp_I, seed_not_I = self.Response(ch,I)

        sigma = [aux,c,rsp_I,seed_not_I]
        return sigma

##################################################################  

#VERIFIER
class Verifier:

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
    
    #II. VERIFICATION
    def Verification(self,Htr_unsys,s,aux,c,rsp_I,seed_not_I):

        #Checking the aux from seeds for i!=I
        aux_bar = [""]*N
        b = 1

        for i in range(0,N):
            if i != self.I_verifier:
                aux_bar[i] = Helper.H(seed_not_I[i])
                if aux_bar[i] != aux[i]:
                    b = 0
                    break
        
        #If b = 1, Checking the rsp for i=I
        if b == 1:
            print(Bcolors.OKCYAN+"\nSimulated Helper is reliable!"+Bcolors.ENDC)
            b_prime = self.V2(Htr_unsys,s,aux[self.I_verifier],c[self.I_verifier],rsp_I)
            if b_prime == 1:
                print(Bcolors.OKCYAN+"\nProver authentication successful!"+Bcolors.ENDC)
            else:
                print(Bcolors.WARNING+"\nProver authentication failed..."+Bcolors.ENDC)
        else:
            print(Bcolors.WARNING+"\nInvalid simulated Helper..."+Bcolors.ENDC)

        return b and b_prime

    #Function the Verifier uses to reconstruct the challenge (is the same as Signer.Challenge())
    def Challenge_verifier(self,mex,aux,c):

        set_random_seed(Integer(hash(str(mex)+str(aux)+str(c)),16))
        self.I_verifier = Integers(N).random_element()
        self.z = Fq.random_element()

        return None

    #Main execution of the Verifier
    def main(self,Htr_unsys,s,sigma,mex):

        aux_verifier = sigma[0]
        c_verifier = sigma[1]
        rsp_I_verifier = sigma[2]
        seed_not_I_verifier = sigma[3]

        self.Challenge_verifier(mex,aux_verifier,c_verifier)

        return self.Verification(Htr_unsys,s,aux_verifier,c_verifier,rsp_I_verifier,seed_not_I_verifier)

##################################################################  

def MPC_id(e,Htr_unsys,s,mex):
    
    signer = Signer()
    verifier = Verifier()

    #SIGNATURE
    sigma = Signer.main(signer,Htr_unsys,e,mex)

    #VERIFICATION
    return Verifier.main(verifier,Htr_unsys,s,sigma,mex)