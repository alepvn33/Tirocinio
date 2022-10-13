reset()

import hashlib
import numpy
load('MPC_utils.sage')


#Scheme parameters
q = 31; #Select prime integer q (size of the field)
n = 256
r = 204
w = 10; #Private key weight
_lambda = 128 #Seed bit lenght
N = 100 #Number of instances
mex = "Ciao a tutti, sono il Signer" #Message to be signed

#Preparations to execute the protocol
Fq=GF(q); #finite field with q elements
Fq_set = Set(Fq)
Fq_star = Fq_set.difference([0]) #multiplicative group of Fq

#Key generation
e, Htr_unsys, s = key_gen()

#MPC-in-the-head identification
ok = MPC_id(e,Htr_unsys,s,mex)
if ok==1:
    print(Bcolors.BOLD+Bcolors.UNDERLINE+Bcolors.OKGREEN+"\nSignature ACCEPTED!\n"+Bcolors.ENDC+Bcolors.ENDC)
else:
    print(Bcolors.BOLD+Bcolors.UNDERLINE+Bcolors.FAIL+"\nSignature REJECTED!\n"+Bcolors.ENDC+Bcolors.ENDC)
