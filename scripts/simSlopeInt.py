import sys
from anthSmut import simAndInfer
from random import random


for k in range(0,10000):
    betaF = 6.#20*random()
    betaVJ = 0.#2.3 #20*random()
    a = simAndInfer.simsInf(readQ=True,betaF_fac=betaF,betaVJ_fac=betaVJ)
    a.simSlopeInt()
    #a.simTransect()
