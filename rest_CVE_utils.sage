#Key generation
def key_gen(Fq,m,n,r):
    set_random_seed();
    Htr_unsys = random_matrix(Fq,n-r,r); #public matrix (only non-systematic portion)
    e = rnd_restricted_vector(Fq,m,n);
    
    #compute syndrome
    s = e[0,0:r] + e[0,r:n]*Htr_unsys;
    return e,Htr_unsys,s;    

##################################################################

#Generate a random vector made of elements in the restricted set (zero=True (default) if 0 can be included in the vector, False if not)
def rnd_restricted_vector(Fq,m,n,zero=True):

    #Define estricted set
    E = [Fq(m)^i for i in range(multiplicative_order(Fq(m)))];
    if zero:    E.insert(0,0); #Add the 0 element to the set
    #Generate random restricted vector
    a = matrix(Fq,1,n);
    for i in range(0,n):
        a[0,i] = choice(E);
    
    return a;

##################################################################
 

#Apply restricted monomial transformation
def apply_rest_monomial(Fq,tau_perm, tau_values, a, n):
    b = matrix(Fq,1,n);
    for i in range(0,n):
        p = tau_perm[i];
        b[0,i] = tau_values[0,p]*a[0,p];
    return b;

##################################################################

#Apply inverse of restricted monomial transformation
def apply_inv_rest_monomial(Fq,tau_perm, tau_values, a, n):
    b = matrix(Fq,1,n);
    for i in range(0, n):
        p = tau_perm[i];
        b[0,p] = a[0,i]/tau_values[0,p];
    return b;    

##################################################################

#Convert a hexadecimal string into a binary list (from LSB to MSB)
def hex_to_bin(hex):
    dec = Integer(hex, 16); #From HEX to DEC
    bin = dec.digits(2); #From HEX to BIN
    return bin

##################################################################

##Simulate one round of the protocol
def one_round_sim(Fq,m,n,r,e,Htr_unsys,s):
    
    ok=0;
    
    ##Generating committments
    u = random_matrix(Fq,1,n);
    set_random_seed();
    tau_seed = initial_seed();
    tau_perm = P.random_element();
    tau_values  = rnd_restricted_vector(Fq,m,n,False); #0 cannot be included in the value vector

    tau_e = apply_rest_monomial(Fq,tau_perm,tau_values,e,n);
    tau_u = apply_rest_monomial(Fq,tau_perm,tau_values,u,n);

    u_Htr = u[0,0:r]+u[0,r:n]*Htr_unsys;


    #Hashing and sending c0 and c1
    c0_before_hash = str(tau_perm)+str(tau_values)+str(u_Htr);
    c0 = hashlib.sha256();
    c0.update(c0_before_hash.encode('utf-8'));
    c0 = c0.hexdigest();

    c1_before_hash = str(tau_u)+str(tau_e);
    c1 = hashlib.sha256();
    c1.update(c1_before_hash.encode('utf-8'));
    c1 = c1.hexdigest();

    ##Verifier chooses z
    z = Fq_star.random_element();

    #Prover computes y
    y = tau_u + z*tau_e;

    ##Verifier chooses b
    b = GF(2).random_element();
    b=0;

    #Creating response
    if b==0:
        ##Calculating tau from seed
        set_random_seed(tau_seed);
        tau_perm_verifier = P.random_element();
        tau_values_verifier  = rnd_restricted_vector(Fq,m,n,False); #0 cannot be included in the value vector
        
        tau_inv_y = apply_inv_rest_monomial(Fq,tau_perm_verifier,tau_values_verifier,y,n);
        final_val = tau_inv_y[0,0:r]+tau_inv_y[0,r:n]*Htr_unsys-z*s;
        prover_c0_before_hash = str(tau_perm_verifier)+str(tau_values_verifier)+str(final_val);
        prover_c0 = hashlib.sha256();
        prover_c0.update(prover_c0_before_hash.encode('utf-8'));
        prover_c0 = prover_c0.hexdigest();
        if prover_c0 == c0:
            ok=1;
    else:
        final_val = y-z*tau_e;
        prover_c1_before_hash = str(final_val)+str(tau_e);
        prover_c1 = hashlib.sha256();
        prover_c1.update(prover_c1_before_hash.encode('utf-8'));
        prover_c1 = prover_c1.hexdigest();
        if prover_c1 == c1:
            ok=1;
    return ok;

########################################################################

def multiple_rounds_sim(Fq,mex,m,n,r,e,Htr_unsys,s,N):


    ##Generating N committments
    big_c_before_hash = str(mex); #it is the overall commitment, initialized with the message to be signed
    tau_u_matrix = matrix(Fq,N,n); #it contains the N vectors u, for all rounds
    tau_perm_matrix = matrix(ZZ,N,n); #it contains the N permutations, for all rounds
    tau_values_matrix = matrix(Fq,N,n); #it contains the N scaling vectors, for all rounds
    tau_e_matrix = matrix(Fq,N,n);
    tau_seed_matrix = matrix(ZZ,N,1); #it contains the N tau seeds, for all rounds

    comm_hashes=[];    
    #Proceeding with remaining rounds
    for i in range(0,N):

        #Generating a random u and a random restricted monomial
        u = random_matrix(Fq,1,n);
        set_random_seed();
        tau_seed = initial_seed();
        tau_perm = P.random_element();
        tau_values = rnd_restricted_vector(Fq,m,n,False); #0 cannot be included in the value vector

        #Apply monomial to tau_e and tau_u
        tau_e = apply_rest_monomial(Fq,tau_perm,tau_values,e,n);
        tau_u = apply_rest_monomial(Fq,tau_perm,tau_values,u,n);

        u_Htr = u[0,0:r]+u[0,r:n]*Htr_unsys;

        #Hashing, appending commitments to c and storing c0 and c1
        c0_before_hash = str(matrix(tau_perm))+str(tau_values)+str(u_Htr);
        c0 = hashlib.sha256();
        c0.update(c0_before_hash.encode('utf-8'));
        c0 = c0.hexdigest();

        c1_before_hash = str(tau_u)+str(tau_e);
        c1 = hashlib.sha256();
        c1.update(c1_before_hash.encode('utf-8'));
        c1 = c1.hexdigest();

        big_c_before_hash+=c0+c1;
        comm_hashes.append(c0); 
        comm_hashes.append(c1);

        #Storing the random generated values
        tau_u_matrix[i,:] = tau_u;
        tau_e_matrix[i,:] = tau_e;
        tau_perm_matrix[i,:] = matrix(tau_perm);
        tau_values_matrix[i,:] = tau_values;
        tau_seed_matrix[i,:] =  tau_seed;

    #Final hashed commitment
    big_c = hashlib.sha256();
    big_c.update(big_c_before_hash.encode('utf-8'));
    big_c = big_c.hexdigest();
    big_c_bin = hex_to_bin(big_c);

    #Re-Hashing of final commitment (if bits frome one digest are not enough)
    #2 bits needed per round, every digest has 256 bits: considering the first round already done (so -1) and considering that int() approximates down and we need to approximate up (so +1),
    #the still needed number of hashes is int(2*N/256) - 1 + 1 = int(2*N/256), which is the upper bound of range() if the lower is 0.
    big_c_last = big_c;
    for i in range(0,int(2*N/256)): 
        big_c_next = hashlib.sha256();
        big_c_next.update(big_c_last.encode('utf-8'));
        big_c_next = big_c_next.hexdigest();
        big_c_next_bin = hex_to_bin(big_c_next);
        big_c_bin = Integer(''.join(str(e) for e in big_c_bin)+''.join(str(e) for e in big_c_next_bin)).digits(2); #Hashed commitment is converted from HEX to BIN ('list' object from LSB to MSB)
        big_c_last = big_c_next;
    
    
    #verification over N rounds starts    
    verifier_c_before_hash= str(mex); #The overall verifier commitment, initialized with the signed message

    for j in range(0,N):

        ##Verifier chooses z
        z = big_c_bin[2*j];

        #Prover computes y
        y = tau_u_matrix[j,:]+z*tau_e_matrix[j,:];

        ##Challenge bit b is chosen from hashed commitment (from LSB to MSB)
        b = big_c_bin[2*j+1];

        #Creating response: for each value of b, we also send the opposite hash commitment
        if b==0:
            #Calculating tau from seeds
            tau_seed = tau_seed_matrix[j,0];
            set_random_seed(tau_seed);
            tau_perm_verifier = P.random_element();
            tau_values_verifier  = rnd_restricted_vector(Fq,m,n,False); #0 cannot be included in the value vector
            
            tau_inv_y = apply_inv_rest_monomial(Fq,vector(tau_perm_verifier),tau_values_verifier,y,n);
            final_val = tau_inv_y[0,0:r]+tau_inv_y[0,r:n]*Htr_unsys-z*s;
            prover_c0_before_hash = str(matrix(tau_perm_verifier))+str(tau_values_verifier)+str(final_val);
            prover_c0 = hashlib.sha256();
            prover_c0.update(prover_c0_before_hash.encode('utf-8'));
            prover_c0 = prover_c0.hexdigest();
            verifier_c_before_hash+=prover_c0+comm_hashes[2*j+1];
        else:
            final_val = y-z*tau_e_matrix[j,:];
            prover_c1_before_hash = str(final_val)+str(tau_e_matrix[j,:]);
            prover_c1 = hashlib.sha256();
            prover_c1.update(prover_c1_before_hash.encode('utf-8'));
            prover_c1 = prover_c1.hexdigest();
            verifier_c_before_hash+=comm_hashes[2*j]+prover_c1;
            
    #Final hashed verified commitment
    verifier_c = hashlib.sha256();
    verifier_c.update(verifier_c_before_hash.encode('utf-8'));
    verifier_c = verifier_c.hexdigest();

    #Final test
    if big_c == verifier_c:
        ok = 1;
    else:
        ok = 0;
    return(ok);    