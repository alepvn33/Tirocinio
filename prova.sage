#prova

q = 5

Fq=GF(q); #finite field with q elements

e       = vector(Fq,[0,1,3,0,0,1,4,2,0,0])
e_tilde = vector(Fq,[0,0,0,2,3,4,1,0,2,0])

supp_e = []
supp_e_tilde = []
for i in range(0,10):
    if e[i]!=0:
        supp_e.append(i)
    if e_tilde[i]!=0:
        supp_e_tilde.append(i)

P_supp = Permutations(supp_e_tilde)
tau_perm_supp = P_supp.random_element()

zeros_e_tilde = [x for x in list(range(0,10)) if x not in supp_e_tilde] #Complement of e_tilde support
P_zeros = Permutations(zeros_e_tilde)
tau_perm_zeros = P_zeros.random_element()

#Joining the permutations of e_tilde support and e_tilde support complementary, according to e support
tau_perm = [] #Total permutation
j = 0;
k = 0;
for i in range(0,10):
    if i in supp_e:
        tau_perm.append(tau_perm_supp[j])
        j = j + 1
    else:
        tau_perm.append(tau_perm_zeros[k])
        k = k + 1

tau_values = matrix(Fq,1,10);
tau_e_tilde = matrix(Fq,1,10);
for i in range(0,10):
    #Find value vector through inverse calculation
    if e_tilde[tau_perm[i]] == 0:
        tau_values[0,tau_perm[i]] = Fq.random_element()
    else:
        tau_values[0,tau_perm[i]] = Fq(e[i]/e_tilde[tau_perm[i]])
    #Calculate tau(e_tilde) just for checking it's equal to e
    tau_e_tilde[0,i] = tau_values[0,tau_perm[i]]*e_tilde[tau_perm[i]]


if e == tau_e_tilde[0]:
    print("DAJEEEE")
    print(tau_perm,tau_values)
else:
    print("NONEE")
    print(e,tau_e_tilde)