reset()

import hashlib
import numpy
import sys
from time import perf_counter

load('/home/sage/OneDrive/Desktop/Cose per tirocinio/Tirocinio/MPC/Opt Signature scheme/opt_MPC_utils.sage')
load('/home/sage/OneDrive/Desktop/Cose per tirocinio/Tirocinio/MPC/Opt Signature scheme/merkletools.sage')

#Scheme parameters
q = 31;         #Select prime integer q (size of the field)
n = 256
r = 204
w = 10;          #Private key weight
_lambda = 128    #Seed bit lenght
N = 400          #Number of instances
s_num = 20       #Number of chosen rounds (|S|)
mex = "Ciao!"
avg = True       #True if average measurements over 5 executions are needed, False otherwise

#Preparations to execute the protocol
Fq=GF(q) #finite field with q elements
Fq_set = Set(Fq)
Fq_star = Fq_set.difference([0]) #multiplicative group of Fq

#Key generation
e, Htr_unsys, s, pk_seed = key_gen()

#MPC-in-the-head identification
ok, sig_size, tot_time = MPC_id(e,Htr_unsys,s,mex,s_num)
if avg:
    for i in range(4):
        ok_tmp, sig_size_tmp, tot_time_tmp = MPC_id(e,Htr_unsys,s,mex,s_num)
        sig_size += sig_size_tmp
        tot_time += tot_time_tmp
    sig_size = sig_size/5
    tot_time = tot_time/5




print("\nSigned message: "+mex)
print("\nSignature size: "+str(sig_size/1024.)+" KB")
print("\nPublic Key size (seed): "+str(sys.getsizeof(pk_seed))+" B")
if ok==1:
    print(Bcolors.BOLD+Bcolors.UNDERLINE+Bcolors.OKGREEN+"\nSignature ACCEPTED!\n"+Bcolors.ENDC+Bcolors.ENDC);
else:
    print(Bcolors.BOLD+Bcolors.UNDERLINE+Bcolors.FAIL+"\nSignature REJECTED!\n"+Bcolors.ENDC+Bcolors.ENDC);



print("Total execution time: "+str(tot_time)+" seconds\n")