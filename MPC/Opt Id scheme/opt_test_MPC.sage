reset();

import hashlib
import numpy
from time import perf_counter

load('/home/sage/OneDrive/Desktop/Cose per tirocinio/Tirocinio/MPC/Opt Id scheme/opt_MPC_utils.sage')
load('/home/sage/OneDrive/Desktop/Cose per tirocinio/Tirocinio/MPC/Opt Id scheme/merkletools.sage')

#Scheme parameters
q = 31;         #Select prime integer q (size of the field)
n = 256;
r = 204;
w = 10;          #Private key weight
_lambda = 128    #Seed bit lenght
N = 400          #Number of instances
s_num = 20       #Number of chosen rounds (|S|)

#Preparations to execute the protocol
Fq=GF(q); #finite field with q elements
Fq_set = Set(Fq);
Fq_star = Fq_set.difference([0]); #multiplicative group of Fq

start_time = perf_counter()

#Key generation
e, Htr_unsys, s = key_gen();

#MPC-in-the-head identification
ok = MPC_id(e,Htr_unsys,s,s_num);
if ok==1:
    print(Bcolors.BOLD+Bcolors.UNDERLINE+Bcolors.OKGREEN+"\nId ACCEPTED!\n"+Bcolors.ENDC+Bcolors.ENDC);
else:
    print(Bcolors.BOLD+Bcolors.UNDERLINE+Bcolors.FAIL+"\nId REJECTED!\n"+Bcolors.ENDC+Bcolors.ENDC);

stop_time = perf_counter()

print(stop_time - start_time, "seconds")