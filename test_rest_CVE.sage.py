

# This file was *autogenerated* from the file test_rest_CVE.sage
from sage.all_cmdline import *   # import sage library

_sage_const_31 = Integer(31); _sage_const_2 = Integer(2); _sage_const_256 = Integer(256); _sage_const_204 = Integer(204); _sage_const_135 = Integer(135); _sage_const_0 = Integer(0); _sage_const_1 = Integer(1)#one round simulation    
reset();

import hashlib;
from sage.misc.prandom import choice
load('rest_CVE_utils.sage')

#scheme parameters
# Select prime integer q (size of the field)
q = _sage_const_31 ;
# Select power base m which generates the restricted set
m = _sage_const_2 ;
n = _sage_const_256 ;
r = _sage_const_204 ;
N = _sage_const_135 ;

#Preparations to execute the protocol
Fq=GF(q); #finite field with q elements
Fq_set = Set(Fq);
Fq_star = Fq_set.difference([_sage_const_0 ]); #multiplicative group of Fq
P = Permutations(range(_sage_const_0 ,n)); #set of permutations of {0, ... , n-1}

#Key generation
e, Htr_unsys, s = key_gen(Fq,m,n,r);

#Single round verification
print("Checking one single round...");
ok = one_round_sim(Fq,m,n,r,e,Htr_unsys,s);
if ok==_sage_const_1 :
    print("ACCEPTED!");
else:
    print("REJECTED!");
    
print("----------------------------");


#Multiple rounds (with compression) verification    
print("Running N rounds with compression...");
ok = multiple_rounds_sim(Fq,m,n,r,e,Htr_unsys,s,N);
if ok==_sage_const_1 :
    print("ACCEPTED!");
else:
    print("REJECTED!");

