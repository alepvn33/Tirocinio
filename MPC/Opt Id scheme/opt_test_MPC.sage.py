

# This file was *autogenerated* from the file opt_test_MPC.sage
from sage.all_cmdline import *   # import sage library

_sage_const_31 = Integer(31); _sage_const_256 = Integer(256); _sage_const_204 = Integer(204); _sage_const_10 = Integer(10); _sage_const_128 = Integer(128); _sage_const_400 = Integer(400); _sage_const_20 = Integer(20); _sage_const_0 = Integer(0); _sage_const_1 = Integer(1)
reset();

import hashlib
import numpy
from time import perf_counter

load('/home/sage/OneDrive/Desktop/Cose per tirocinio/Tirocinio/MPC/Opt Id scheme/opt_MPC_utils.sage')
load('/home/sage/OneDrive/Desktop/Cose per tirocinio/Tirocinio/MPC/Opt Id scheme/merkletools.sage')

#Scheme parameters
q = _sage_const_31 ;         #Select prime integer q (size of the field)
n = _sage_const_256 ;
r = _sage_const_204 ;
w = _sage_const_10 ;          #Private key weight
_lambda = _sage_const_128     #Seed bit lenght
N = _sage_const_400           #Number of instances
s_num = _sage_const_20        #Number of chosen rounds (|S|)

#Preparations to execute the protocol
Fq=GF(q); #finite field with q elements
Fq_set = Set(Fq);
Fq_star = Fq_set.difference([_sage_const_0 ]); #multiplicative group of Fq

start_time = perf_counter()

#Key generation
e, Htr_unsys, s = key_gen();

#MPC-in-the-head identification
ok = MPC_id(e,Htr_unsys,s,s_num);
if ok==_sage_const_1 :
    print(Bcolors.BOLD+Bcolors.UNDERLINE+Bcolors.OKGREEN+"\nId ACCEPTED!\n"+Bcolors.ENDC+Bcolors.ENDC);
else:
    print(Bcolors.BOLD+Bcolors.UNDERLINE+Bcolors.FAIL+"\nId REJECTED!\n"+Bcolors.ENDC+Bcolors.ENDC);

stop_time = perf_counter()

print(stop_time - start_time, "seconds")

