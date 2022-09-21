#one round simulation    
reset();

import hashlib;
from sage.misc.prandom import choice
load('rest_CVE_utils.sage')

#scheme parameters
q = 31; # Select prime integer q (size of the field)
m = 2; # Select power base m which generates the restricted set
n = 256;
r = 204;
N = 135;
mex = "Ciao core"; # Message to be signed

#Preparations to execute the protocol
Fq=GF(q); #finite field with q elements
Fq_set = Set(Fq);
Fq_star = Fq_set.difference([0]); #multiplicative group of Fq
P = Permutations(range(0,n)); #set of permutations of {0, ... , n-1}

#Key generation
e, Htr_unsys, s = key_gen(Fq,m,n,r);

#Single round verification
print("Checking one single round...");
ok = one_round_sim(Fq,m,n,r,e,Htr_unsys,s);
if ok==1:
    print("ACCEPTED!");
else:
    print("REJECTED!");
    
print("----------------------------");


#Multiple rounds (with compression) verification    
print("Signed message: "+str(mex)+"\nRunning N rounds with compression...");
ok = multiple_rounds_sim(Fq,mex,m,n,r,e,Htr_unsys,s,N);
if ok==1:
    print("Signature ACCEPTED!");
else:
    print("Signature REJECTED!");
