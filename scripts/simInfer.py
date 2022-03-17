from anthSmut import simAndInfer

for i in range(0,1):
    a = simAndInfer.simsInf(readQ=True,betaF_fac=10.,betaVJ_fac=10.)
    a.simTransect()
