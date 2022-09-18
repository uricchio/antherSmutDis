import sys
from anthSmut import simAndInfer
from random import random

betaF = 20.*random()
betaVJ = 20.*random()

a = simAndInfer.simsInf(readQ=True,betaF_fac=betaF,betaVJ_fac=betaVJ,migSpores=0.5,migSeeds=0.5)
print ("#", betaF, betaVJ)
a.simTransectRepGen()
