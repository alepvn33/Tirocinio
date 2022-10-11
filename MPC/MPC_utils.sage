#MPC functions

#Key generation
def key_gen(Fq,n,r,w):

    set_random_seed()
    Htr_unsys = random_matrix(Fq,n-r,r) #public matrix (only non-systematic portion)
    e = FWV(n,w)

    #compute syndrome
    s = e[0:r] + e[r:n]*Htr_unsys

    return e,Htr_unsys,s

##################################################################

#Generate a random vector with weight w
def FWV(n,w):

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
def find_isometry(e,e_tilde,n):

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
def apply_rest_monomial(Fq,tau_perm,tau_values,a,n):

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

#SETUP (H)
def H(seed):

    u_helper = random_vector(Fq,n)
    e_tilde_helper = FWV(n,w)

    set_random_seed(seed)
    r_v = matrix(GF(2),q,_lambda)
    c_v = [""]*len(Fq)

    for v in Fq:
        r_v[v] = random_vector(GF(2),_lambda)
        c_v[v] = hash(str(r_v[v])+str(u_helper+v*e_tilde_helper))

    aux = c_v

    return aux     

##################################################################  

#PROVER
class Prover:

    #COMMITMENT (P1)
    def P1(self,Htr_unsys,e,seed):

        set_random_seed(seed)
        self.u = random_vector(Fq,n)
        self.e_tilde = FWV(n,w)
        self.tau_perm, self.tau_values = find_isometry(e,self.e_tilde,n)
        self.r_vec = random_vector(GF(2),_lambda)
        tau_u = apply_rest_monomial(Fq,self.tau_perm,self.tau_values,self.u,n)
        c = hash(str(self.r_vec)+str(self.tau_perm)+str(self.tau_values)+str(tau_u[0:r]+tau_u[r:n]*Htr_unsys))

        return c

    #RESPONSE (P2)
    def P2(self,ch,seed):

        set_random_seed(seed)
        r_v_prover = matrix(GF(2),q,_lambda)

        for v in Fq:
            r_v_prover[v] = random_vector(GF(2),_lambda)

        r_z = r_v_prover[ch]
        y = self.u + ch*self.e_tilde
        rsp = [self.r_vec,r_z,self.tau_perm,self.tau_values,y]

        return rsp

##################################################################  

#VERIFIER
class Verifier:

    #CHALLENGE (V1)
    def V1(self):

        set_random_seed()
        self.z = Fq.random_element()
        ch = self.z

        return ch
    
    #VERIFICATION (V2)
    def V2(self,Htr_unsys,s,aux,c,rsp):

        tau_y = apply_rest_monomial(Fq,rsp[2],rsp[3],rsp[4],n)
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

##################################################################  


def MPC_id(e,Htr_unsys,s):
    
    set_random_seed()
    seed = initial_seed()

    #SETUP (H)
    aux = H(seed)
    
    #COMMITMENT (P1)
    prover = Prover()
    c = Prover.P1(prover,Htr_unsys,e,seed)

    #CHALLENGE (V1)
    verifier = Verifier()
    ch = Verifier.V1(verifier)

    #RESPONSE (P2)
    rsp = Prover.P2(prover,ch,seed)

    #VERIFICATION (V2)
    return Verifier.V2(verifier,Htr_unsys,s,aux,c,rsp)