import sys
from anthSmut import simAndInfer


betaF = float(sys.argv[1])
betaVJ = float(sys.argv[2])

for i in range(0,1):
    a = simAndInfer.simsInf(readQ=False,betaF_fac=betaF,betaVJ_fac=betaVJ)
    a.idealQual()
    a.simTransect(printLinTrans=True,nInf=50,maxGens=100)


