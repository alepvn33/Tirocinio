reset();

import hashlib
import numpy
load('MPC_utils.sage')

#Scheme parameters
q = 31; #Select prime integer q (size of the field)
n = 256;
r = 204;
w = 10; #Private key weight
_lambda = 128 #Seed bit lenght
N = 100 #Number of instances

#Preparations to execute the protocol
Fq=GF(q); #finite field with q elements
Fq_set = Set(Fq);
Fq_star = Fq_set.difference([0]); #multiplicative group of Fq

#Key generation
e, Htr_unsys, s = key_gen(Fq,n,r,w);

#MPC-in-the-head identification
ok = MPC_id(e,Htr_unsys,s);
if ok==1:
    print("Id ACCEPTED!");
else:
    print("Id REJECTED!");
