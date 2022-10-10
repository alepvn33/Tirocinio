#Key generation
from hmac import digest


def key_gen(Fq,n,r,w):
    set_random_seed()
    Htr_unsys = random_matrix(Fq,n-r,r) #public matrix (only non-systematic portion)
    e = FWV(Fq,n,w)
    
    #compute syndrome
    s = e[0:r] + e[r:n]*Htr_unsys
    return e,Htr_unsys,s

##################################################################

#Generate a random vector with weight w
def FWV(Fq,n,w):
    a = vector(Fq,w)
    Fq_star = Set(Fq).difference([0])
    for i in range(0,w):
        a[i] = choice(Fq_star);
    b = zero_vector(n) #List of zeros of lenght n
    P = Permutations(range(0,n));
    rnd_supp = P.random_element()
    for i in range(0,w):
        b[rnd_supp[i]] = a[i]
    return b

##################################################################

#Find an isometry function tau such that e = tau(e_tilde)
def find_isometry(e,e_tilde,n):
    #print("wth(e) = {}, wth(e_tilde) = {}".format(numpy.count_nonzero(e),numpy.count_nonzero(e_tilde)))
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
    j = 0;
    k = 0;
    for i in range(0,n):
        if i in supp_e:
            #print('j = {}, len(tau_perm_supp) = {}'.format(j,len(tau_perm_supp)))
            tau_perm.append(tau_perm_supp[j])
            j = j + 1
        else:
            #print('k = {}, len(tau_perm_zeros) = {}'.format(k,len(tau_perm_zeros)))
            tau_perm.append(tau_perm_zeros[k])
            k = k + 1

    tau_values = vector(Fq,n);
    for i in range(0,n):
        #Find value vector through inverse calculation
        if e_tilde[tau_perm[i]] == 0:
            tau_values[tau_perm[i]] = Fq.random_element()
        else:
            tau_values[tau_perm[i]] = Fq(e[i]/e_tilde[tau_perm[i]])
    
    return tau_perm, tau_values

##################################################################   

#Apply monomial transformation (isometry)
def apply_rest_monomial(Fq,tau_perm,tau_values,a,n):
    b = vector(Fq,n);
    for i in range(0,n):
        p = tau_perm[i];
        b[i] = tau_values[p]*a[p];
    return b;

##################################################################   

#Hash input
def hash(input):
    digest = hashlib.sha256()
    digest.update(input.encode('utf-8'))
    digest = digest.hexdigest()
    return digest;

##################################################################         

def MPC_id(Fq,n,r,w,e,Htr_unsys,s):
    
    #SETUP (H)
    set_random_seed()
    seed = initial_seed()
    _lambda = 128 #Seed bit-lenght
    u_helper = random_vector(Fq,n)
    e_tilde_helper = FWV(Fq,n,w)

    set_random_seed(seed)
    r_v = matrix(GF(2),q,_lambda)
    c_v = [""]*len(Fq)
    for v in Fq:
        r_v[v] = random_vector(GF(2),_lambda)
        c_v[v] = hash(str(r_v[v])+str(u_helper+v*e_tilde_helper))
    aux = c_v

    #COMMITMENT (P1)
    set_random_seed(seed)
    u = random_vector(Fq,n)
    e_tilde = FWV(Fq,n,w)
    tau_perm, tau_values = find_isometry(e,e_tilde,n)
    r_vec = random_vector(GF(2),_lambda)
    tau_u = apply_rest_monomial(Fq,tau_perm,tau_values,u,n)
    c = hash(str(r_vec)+str(tau_perm)+str(tau_values)+str(tau_u[0:r]+tau_u[r:n]*Htr_unsys))

    #CHALLENGE (V1)
    set_random_seed()
    z = Fq.random_element()
    ch = z

    #RESPONSE (P2)
    set_random_seed(seed)
    r_v_prover = matrix(GF(2),q,_lambda)
    for v in Fq:
        r_v_prover[v] = random_vector(GF(2),_lambda)
    r_z = r_v_prover[ch]
    y = u + ch*e_tilde
    rsp = [r_vec,r_z,tau_perm,tau_values,y]

    #VERIFICATION (V2)
    tau_y = apply_rest_monomial(Fq,rsp[2],rsp[3],rsp[4],n)
    t = tau_y[0:r]+tau_y[r:n]*Htr_unsys - z*s

    c_verifier = hash(str(rsp[0])+str(rsp[2])+str(rsp[3])+str(t))
    

    wth_e_tilde = numpy.count_nonzero(e_tilde)

    if c_verifier == c and wth_e_tilde == w:

        c_z_verifier = hash(str(rsp[1])+str(rsp[4]))
        c_z = aux[z]

        if c_z_verifier == c_z:
            return 1
        else:
            return 0
    else:
        return 0