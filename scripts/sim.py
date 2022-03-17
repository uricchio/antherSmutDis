import sys
from anthSmut import simAndInfer


betaF = float(sys.argv[1])
betaVJ = float(sys.argv[2])

for i in range(0,1000):
    a = simAndInfer.simsInf(readQ=True,betaF_fac=betaF,betaVJ_fac=betaVJ)
    a.simTransect()
