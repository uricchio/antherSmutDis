import numpy as np
import sys
import math
from scipy.stats import poisson
from scipy.optimize import minimize
from scipy.optimize import brute

class simsInf:

    def __init__(self,readQ = True, determ = False,betaF_fac=1.,betaVJ_fac=1.):

        self.readQ = readQ
        self.determ = determ

        # containers for spatial grids
        self.Vhgrid = np.zeros((15,100)) # veg healthy
        self.Fhgrid = np.zeros((15,100)) # flow healthy
        self.Fdgrid = np.zeros((15,100)) # flow disease
        self.Vdgrid = np.zeros((15,100)) # veg disease
        self.jgrid = np.zeros((15,100))  # juvenile 
        self.qGrid = np.ones((15,100))  # quality function
        self.qGridOrig = np.ones((15,100))  # quality function original
        self.N = np.zeros((15,100))        
        
        # containers for spatial grids prev gen
        self.Vhgrid0 = np.zeros((15,100)) # veg healthy
        self.Fhgrid0 = np.zeros((15,100)) # flow healthy
        self.Fdgrid0 = np.zeros((15,100)) # flow disease
        self.Vdgrid0 = np.zeros((15,100)) # veg disease
        self.jgrid0 = np.zeros((15,100))  # juvenile 

        # containers for real data
        self.Fh15 = np.zeros((15,100))
        self.Fd15 = np.zeros((15,100))
        self.Fh14 = np.zeros((15,100))
        self.Fd14 = np.zeros((15,100))
        self.Fh05 = np.zeros((15,100))
        self.Fd05 = np.zeros((15,100))
        self.Fh07 = np.zeros((15,100))
        self.Fd07 = np.zeros((15,100))
        self.Fh18 = np.zeros((15,100))
        self.Fd18 = np.zeros((15,100))
        self.Fh19 = np.zeros((15,100))
        self.Fd19 = np.zeros((15,100))
  
        #data files
        self.AllData = 'empDensData.txt'
        self.dataDens15 = 'data/dens.15.txt' 
        self.dataDis15 = 'data/dis.15.txt' 
        self.dataDens14 = 'data/dens.14.txt' 
        self.dataDis14 = 'data/dis.14.txt' 
        self.dataDens05 = 'data/dens.05.txt' 
        self.dataDis05 = 'data/dis.05.txt' 
        self.dataDens07 = 'data/dens.07.txt' 
        self.dataDis07 = 'data/dis.07.txt' 
        self.dataDens18 = 'data/dens.18.txt' 
        self.dataDis18 = 'data/dis.18.txt' 
        self.dataDens19 = 'data/dens.19.txt' 
        self.dataDis19 = 'data/dis.19.txt' 

        # "known" parameters
        self.betaF  = betaF_fac*0.069   # freq dependent transmission rate
        self.betaV = betaVJ_fac*0.0000338*250 # density dependent transmission rate in adult plants 
        self.betaJ = betaVJ_fac*0.000338*250 # density dependent transmission rate in juv plants

        self.rec  = 0.038          # recoverrate
        self.muF = 0.116           # mortality flowering
        self.muV =  0.258          # mortality vegetative         
        self.phiFh = 0.664         # probabiity of flowering
        self.phiFd = 0.704
        self.phiV  = 0.356
        self.phiJ  = 0.17          # calculated from implant -proportion that flowered in first year
        
        # UNKNOWN PARAMETERS
        self.b= 1.3                  # 1 note - this is establishment that takes into account fixed probability of seed/seedling mortality
        self.migSeeds  = 0.1          # migration rate of seeds = fraction of seeds added to neighbors
        self.migSpores = 0.1
        self.k = 0.00015*250

        self.VhProb = []
        self.VdProb = []
        self.FdProb = []
        self.FhProb = []

    def get_N(self):
        
        self.N = np.add(self.jgrid,self.Vdgrid)
        self.N = np.add(self.N,self.Fhgrid)
        self.N = np.add(self.N,self.Fdgrid)
        self.N = np.add(self.N,self.Vhgrid)
    
    def get_N_old(self):
        
        self.N = np.add(self.jgrid0,self.Vdgrid0)
        self.N = np.add(self.N,self.Fhgrid0)
        self.N = np.add(self.N,self.Fdgrid0)
        self.N = np.add(self.N,self.Vhgrid0)
    
        #self.N = np.multiply(2,self.Fhgrid0)
        #self.N = np.add(self.N,np.multiply(2,self.Fdgrid0))

    def get_N_old_Real(self):
        
        #self.N = self.jgrid0[:]
        #self.N = np.add(np.multiply(self.Fhgrid0,np.multiply(1-self.phiJ,Fj)),np.multiply((1-self.phiFd*(1-self.muF))/(self.phiV*(1-self.muV)),self.Fdgrid0))
        self.N = np.multiply(2,self.Fhgrid0)
        self.N = np.add(self.N,np.multiply(2,self.Fdgrid0))

        #self.N = np.add(self.N,np.multiply(np.add(np.multiply(self.Fhgrid0,np.multiply(1-self.phiJ,np.subtract(1,Fj))),np.multiply((1-(1-self.muF)*self.phiFh)/(self.phiV*(1-self.muV)),self.Fhgrid0))))

    def idealQual(self):
        for i in range(0,15):
            for j in range(0,100):
                if j < 10 or j > 90:
                    self.qGrid[i][j] = 0.1
                elif j >= 10 and j < 50:
                    self.qGrid[i][j] = 0.1 + (j-10)*(0.9/40)
                else:
                    self.qGrid[i][j] = self.qGrid[i][100-j]

        #for i in range(0,15):
        #    for j in range(0,100):
        #        print self.qGrid[i][j],
        #    print

    def readQual(self,qFile='data/krumQual.txt',printq=False,trim=True):

        # method to read quality scores in

        # first initialize qGrid
        fh = open(qFile, 'r')
        xMax = 0
        yMax = 0
        for line in fh:
            data = line.strip().split()
            if yMax < float(data[1]):
                yMax = int(round(float(data[1])))
            if xMax < float(data[0]):
                xMax = int(round(float(data[0])))
    
        # next fill in qGrid
        self.qGrid = np.empty((yMax,xMax))
        self.qGrid[:] = np.nan
        fh.close()   

        fh = open(qFile, 'r')
        for line in fh:
            data = line.strip().split()
            if data[2] == 'NA':
                data[2] = np.nan
            self.qGridOrig[int(round(float(data[1])))-1,int(round(float(data[0])))-1] = float(data[2])
            if math.isnan(float(data[2])):
                self.qGrid[int(round(float(data[1])))-1,int(round(float(data[0])))-1] = 0.
            else:
                self.qGrid[int(round(float(data[1])))-1,int(round(float(data[0])))-1] = float(data[2])
        fh.close()
 
        if trim:
            holder = np.empty((15,100))
            for i in range(0,15):
                for j in range(0,100):
                    holder[i,j] = self.qGrid[i,j]
            self.qGrid = holder[:]
        
        if printq:
            for i in range(0,15):
                for j in range(0,100):
                    sys.stdout.write(str(self.qGrid[i,j])+' ')
                sys.stdout.write('\n')    
 
    def readRealData(self):
        fh = open(self.dataDens15, 'r')        
        for line in fh:
            data = line.strip().split()
            self.Fh15[int(data[1])-1][int(data[0])-1]  = int(data[2])
        fh.close()        
        
        fh = open(self.dataDis15, 'r')        
        for line in fh:
            data = line.strip().split()
            self.Fd15[int(data[1])-1][int(data[0])-1]  = int(data[2])
        fh.close()        
        
        fh = open(self.dataDis14, 'r')        
        for line in fh:
            data = line.strip().split()
            self.Fd14[int(data[1])-1][int(data[0])-1]  = int(data[2])
        fh.close()        
        
        fh = open(self.dataDens14, 'r')        
        for line in fh:
            data = line.strip().split()
            self.Fh14[int(data[1])-1][int(data[0])-1]  = int(data[2])
        fh.close()        

        fh = open(self.dataDis07, 'r')
        for line in fh:
            data = line.strip().split()
            if data[2] == 'NA': 
                data[2] = np.nan
                self.Fd07[int(data[1])-1][int(data[0])-1]  = data[2]
            else:  
                self.Fd07[int(data[1])-1][int(data[0])-1]  = int(data[2])
        fh.close()

        fh = open(self.dataDens07, 'r')
        for line in fh:
            data = line.strip().split()
            if data[2] == 'NA': 
                data[2] = np.nan
                self.Fh07[int(data[1])-1][int(data[0])-1]  = data[2]
            else:
                self.Fh07[int(data[1])-1][int(data[0])-1]  = int(data[2])
        fh.close()

        fh = open(self.dataDis05, 'r')
        for line in fh:
            data = line.strip().split()
            if data[2] == 'NA': 
                data[2] = np.nan
                self.Fd05[int(data[1])-1][int(data[0])-1]  = data[2]
            else:
                self.Fd05[int(data[1])-1][int(data[0])-1]  = int(round(float(data[2])))
        fh.close()

        fh = open(self.dataDens05, 'r')
        for line in fh:
            data = line.strip().split()
            if data[2] == 'NA': 
                data[2] = np.nan
                self.Fh05[int(data[1])-1][int(data[0])-1]  = data[2]
            else:
                self.Fh05[int(data[1])-1][int(data[0])-1]  = int(round(float(data[2])))
        fh.close()

        fh = open(self.dataDis18, 'r')
        for line in fh:
            data = line.strip().split()
            if data[2] == 'NA': 
                data[2] = np.nan
                self.Fd18[int(data[1])-1][int(data[0])-1]  = data[2]
            else:
                self.Fd18[int(data[1])-1][int(data[0])-1]  = int(data[2])
        fh.close()

        fh = open(self.dataDens18, 'r')
        for line in fh:
            data = line.strip().split()
            if data[2] == 'NA': 
                data[2] = np.nan
                self.Fh18[int(data[1])-1][int(data[0])-1]  = data[2]
            else:
                self.Fh18[int(data[1])-1][int(data[0])-1]  = int(data[2])
        fh.close()
        
        fh = open(self.dataDis19, 'r')
        for line in fh:
            data = line.strip().split()
            if data[2] == 'NA':
                data[2] = np.nan
                self.Fd19[int(data[1])-1][int(data[0])-1]  = data[2]
            else:
                self.Fd19[int(data[1])-1][int(data[0])-1]  = int(data[2])
        fh.close()

        fh = open(self.dataDens19, 'r')
        for line in fh:
            data = line.strip().split()
            if data[2] == 'NA':
                data[2] = np.nan
                self.Fh19[int(data[1])-1][int(data[0])-1]  = data[2]
            else:
                self.Fh19[int(data[1])-1][int(data[0])-1]  = int(data[2])
        fh.close()



        for i in range(0,len(self.Fh14)):
            for j in range(0,len(self.Fh14[i])):
                self.Fh14[i][j] = self.Fh14[i][j] - self.Fd14[i][j]
                self.Fh15[i][j] = self.Fh15[i][j] - self.Fd15[i][j]
                self.Fh05[i][j] = self.Fh05[i][j] - self.Fd05[i][j]
                self.Fh07[i][j] = self.Fh07[i][j] - self.Fd07[i][j]
                self.Fh18[i][j] = self.Fh18[i][j] - self.Fd18[i][j]
                #self.Fh19[i][j] = self.Fh19[i][j] - self.Fd19[i][j]
                # this is a hack I added because of the issues with the ML inference for these cells
                if self.Fh19[i][j] > 2 and self.Fh18[i][j] == 0:
                    self.Fh19[i][j] = np.nan
                    self.Fh18[i][j] = np.nan

        return

    def initPop(self):
 
        # set each cell to an arbitrary level

        initArb = 50
        self.Vhgrid = np.empty((15,100))
        self.Vdgrid = np.empty((15,100))
        self.Fhgrid = np.empty((15,100))
        self.Fdgrid = np.empty((15,100))
        self.jgrid = np.empty((15,100))

        self.Vhgrid[:] = initArb
        self.Vdgrid[:] = 0.
        self.Fhgrid[:] = initArb
        self.Fdgrid[:] = 0.
        self.jgrid[:] = initArb
    
    def compute_expectations_Real(self,stay0,stay1,betaF,betaV,betaJ):

        # first, get force of infection
        Fj = 0 # for junveniles
        Ff = 0 # for flowering
        Fv = 0 # for vegetative

        # get nearest neighbor counts for juveniles, which are weighted by self.migSeeds
        # contributions, others do not
        nn_seeds = np.zeros((15,100))
        # nn_seeds is effective number of juveniles from neighboring cells
        # get component of juveniles from within the current cell
        stay = stay0
        self.get_N_old_Real()
        nn_seeds = np.multiply(self.Fhgrid0,np.multiply(stay,np.divide(np.multiply(self.b,self.qGrid),(np.add(1.,np.multiply(self.k,self.N))))))
        # get components from neighboring cells -- use focal cell in neighbor does not exist
        for i in range(0,15):
            for j in range(0,100):
                if i == 0 or math.isnan(self.qGridOrig[i-1,j]):
                    nn_seeds[i,j] += np.multiply(self.Fhgrid0[i,j],np.multiply((1-stay)*0.25,np.divide(np.multiply(self.b,self.qGrid[i,j]),np.add(1.,np.multiply(self.k,self.N[i,j])))))
                else:
                    nn_seeds[i,j] +=  np.multiply(self.Fhgrid0[i-1,j],np.multiply((1-stay)*0.25,np.divide(np.multiply(self.b,self.qGrid[i-1,j]),np.add(1.,np.multiply(self.k,self.N[i-1,j])))))
                if i == 14 or math.isnan(self.qGridOrig[i+1,j]):    
                    nn_seeds[i,j] +=  np.multiply(self.Fhgrid0[i,j],np.multiply((1-stay)*0.25,np.divide(np.multiply(self.b,self.qGrid[i,j]),np.add(1.,np.multiply(self.k,self.N[i,j])))))
                else:
                    nn_seeds[i,j] +=  np.multiply(self.Fhgrid0[i+1,j],np.multiply((1-stay)*0.25,np.divide(np.multiply(self.b,self.qGrid[i+1,j]),np.add(1.,np.multiply(self.k,self.N[i+1,j])))))
                if j == 0 or math.isnan(self.qGridOrig[i,j-1]):
                    nn_seeds[i,j] +=  np.multiply(self.Fhgrid0[i,j],np.multiply((1-stay)*0.25,np.divide(np.multiply(self.b,self.qGrid[i,j]),np.add(1.,np.multiply(self.k,self.N[i,j])))))
                else:
                    nn_seeds[i,j] +=  np.multiply(self.Fhgrid0[i,j-1],np.multiply((1-stay)*0.25,np.divide(np.multiply(self.b,self.qGrid[i,j-1]),np.add(1.,np.multiply(self.k,self.N[i,j-1])))))
                if j == 99 or math.isnan(self.qGridOrig[i,j+1]):
                    nn_seeds[i,j] +=  np.multiply(self.Fhgrid0[i,j],np.multiply((1-stay)*0.25,np.divide(np.multiply(self.b,self.qGrid[i,j]),np.add(1.,np.multiply(self.k,self.N[i,j])))))
                else:
                    nn_seeds[i,j] +=  np.multiply(self.Fhgrid0[i,j+1],np.multiply((1-stay)*0.25,np.divide(np.multiply(self.b,self.qGrid[i,j+1]),np.add(1.,np.multiply(self.k,self.N[i,j+1])))))

        # Force of infection and gets weighted by self.migSpores
        stay = stay1
        nn_spores = np.multiply(stay,self.Fdgrid0)
        
        # get components from neighboring cells -- use focal cell in neighbor does not exist
        for i in range(0,15):
            for j in range(0,100):
                if i == 0 or math.isnan(self.qGridOrig[i-1,j]):
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid0[i,j])
                else:
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid0[i-1,j])
                if i == 14 or math.isnan(self.qGridOrig[i+1,j]):    
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid0[i,j])
                else:
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid0[i+1,j])
                if j == 0 or math.isnan(self.qGridOrig[i,j-1]):
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid0[i,j])
                else:
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid0[i,j-1])
                if j == 99 or math.isnan(self.qGridOrig[i,j+1]):
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid0[i,j])
                else:
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid0[i,j+1])

        # calc Ff, which is vector-borne and depends only on whole grid frequency of flowering diseased
        Ff = np.minimum((betaF*np.nansum(self.Fdgrid0))/(np.nansum(self.Fdgrid0)+np.nansum(self.Fhgrid0)),1)

        # now Fj, which is diff for each cell in grid
        Fj = np.minimum(np.multiply(betaJ,nn_spores),np.ones((15,100)))

        # Fv, which is just Fj scaled
        Fv = np.minimum(np.multiply(betaV,nn_spores),np.ones((15,100)))       

        jgridInf = nn_seeds[:]
        VhInf = np.add(np.multiply(jgridInf,np.multiply(1-self.phiJ,np.subtract(1,Fj))),np.multiply((1-(1-self.muF)*self.phiFh)/(self.phiV*(1-self.muV)),self.Fhgrid0))
        VdInf = np.add(np.multiply(jgridInf,np.multiply(1-self.phiJ,Fj)),np.multiply((1-self.phiFd*(1-self.muF))/(self.phiV*(1-self.muV)),self.Fdgrid0))
    
        Fhtemp = np.multiply(jgridInf,np.multiply(self.phiJ,np.subtract(1,Fj)))
        Fhtemp = np.add(Fhtemp,np.multiply(self.Fhgrid0,np.multiply((1-self.muF)*self.phiFh*(1-Ff),np.subtract(1,Fv))))
        Fhtemp = np.add(Fhtemp,np.multiply(self.Fdgrid0,(1-self.muF)*self.phiFd*self.rec))
        Fhtemp = np.add(Fhtemp,np.multiply(VhInf,np.multiply((1-self.muV)*self.phiV,np.subtract(1,Fv))))
        Fhtemp = np.add(Fhtemp,np.multiply(VdInf,(1-self.muV)*self.phiV*self.rec))

        Fdtemp = np.multiply(jgridInf,np.multiply(self.phiJ,Fj))
        Fdtemp = np.add(Fdtemp,np.multiply(self.Fhgrid0,np.multiply((1-self.muF)*self.phiFh,np.add(Ff,Fv))))
        Fdtemp = np.add(Fdtemp,np.multiply(self.Fdgrid0,(1-self.muF)*self.phiFd*(1-self.rec)))
        Fdtemp = np.add(Fdtemp,np.multiply(VhInf,np.multiply((1-self.muV)*self.phiV,Fv)))
        Fdtemp = np.add(Fdtemp,np.multiply(VdInf,(1-self.muV)*self.phiV*(1-self.rec)))

        self.FdProb = np.add(Fdtemp[:],1e-3)
        self.FhProb = np.add(Fhtemp[:],1e-3)

    def compute_expectations(self,stay0,stay1,betaF,betaV,betaJ):

        # first, get force of infection
        Fj = 0 # for junveniles
        Ff = 0 # for flowering
        Fv = 0 # for vegetative

        # get nearest neighbor counts for juveniles, which are weighted by self.migSeeds
        # contributions, others do not
        nn_seeds = np.zeros((15,100))
        # nn_seeds is effective number of juveniles from neighboring cells
        # get component of juveniles from within the current cell

        stay = stay0
        self.get_N_old()
        nn_seeds = np.multiply(self.Fhgrid0,np.multiply(stay,np.divide(np.multiply(self.b,self.qGrid),(np.add(1.,np.multiply(self.k,self.N))))))
        # get components from neighboring cells -- use focal cell in neighbor does not exist
        for i in range(0,15):
            for j in range(0,100):
                if i == 0 or math.isnan(self.qGridOrig[i-1,j]):
                    nn_seeds[i,j] += np.multiply(self.Fhgrid0[i,j],np.multiply(0.25*(1-stay),np.divide(np.multiply(self.b,self.qGrid[i,j]),(np.add(1.,np.multiply(self.k,self.N[i,j]))))))
                else:
                    nn_seeds[i,j] += np.multiply(self.Fhgrid0[i-1,j],np.multiply(0.25*(1-stay),np.divide(np.multiply(self.b,self.qGrid[i-1,j]),(np.add(1.,np.multiply(self.k,self.N[i-1,j]))))))
                if i == 14 or math.isnan(self.qGridOrig[i+1,j]):    
                    nn_seeds[i,j] += np.multiply(self.Fhgrid0[i,j],np.multiply(0.25*(1-stay),np.divide(np.multiply(self.b,self.qGrid[i,j]),(np.add(1.,np.multiply(self.k,self.N[i,j]))))))
                else:
                    nn_seeds[i,j] += np.multiply(self.Fhgrid0[i+1,j],np.multiply(0.25*(1-stay),np.divide(np.multiply(self.b,self.qGrid[i+1,j]),(np.add(1.,np.multiply(self.k,self.N[i+1,j]))))))
                if j == 0 or math.isnan(self.qGridOrig[i,j-1]):
                    nn_seeds[i,j] += np.multiply(self.Fhgrid0[i,j],np.multiply(0.25*(1-stay),np.divide(np.multiply(self.b,self.qGrid[i,j]),(np.add(1.,np.multiply(self.k,self.N[i,j]))))))
                else:
                    nn_seeds[i,j] += np.multiply(self.Fhgrid0[i,j-1],np.multiply(0.25*(1-stay),np.divide(np.multiply(self.b,self.qGrid[i,j-1]),(np.add(1.,np.multiply(self.k,self.N[i,j-1]))))))
                if j == 99 or math.isnan(self.qGridOrig[i,j+1]):
                    nn_seeds[i,j] += np.multiply(self.Fhgrid0[i,j],np.multiply(0.25*(1-stay),np.divide(np.multiply(self.b,self.qGrid[i,j]),(np.add(1.,np.multiply(self.k,self.N[i,j]))))))
                else:
                    nn_seeds[i,j] += np.multiply(self.Fhgrid0[i,j+1],np.multiply(0.25*(1-stay),np.divide(np.multiply(self.b,self.qGrid[i,j+1]),(np.add(1.,np.multiply(self.k,self.N[i,j+1]))))))
        
        # Force of infection and gets weighted by self.migSpores
        stay = stay1
        nn_spores = np.multiply(stay,self.Fdgrid0)
        
        # get components from neighboring cells -- use focal cell in neighbor does not exist
        for i in range(0,15):
            for j in range(0,100):
                if i == 0 or math.isnan(self.qGridOrig[i-1,j]):
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid0[i,j])
                else:
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid0[i-1,j])
                if i == 14 or math.isnan(self.qGridOrig[i+1,j]):    
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid0[i,j])
                else:
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid0[i+1,j])
                if j == 0 or math.isnan(self.qGridOrig[i,j-1]):
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid0[i,j])
                else:
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid0[i,j-1])
                if j == 99 or math.isnan(self.qGridOrig[i,j+1]):
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid0[i,j])
                else:
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid0[i,j+1])

        # calc Ff, which is vector-borne and depends only on whole grid frequency of flowering diseased
        Ff = np.minimum((betaF*np.nansum(self.Fdgrid0))/(np.nansum(self.Fdgrid0)+np.nansum(self.Fhgrid0)),1)

        # now Fj, which is diff for each cell in grid
        Fj = np.minimum(np.multiply(betaJ,nn_spores),np.ones((15,100)))

        # Fv, which is just Fj scaled
        Fv = np.minimum(np.multiply(betaV,nn_spores),np.ones((15,100)))       

        #print np.nansum(Ff),np.nansum(Fv)/1500,np.nansum(Fj) /1500

        # next get expectation of count in each life/infection stage
        # these eacFhtemp = np.multiply(self.jgrid0,np.multiply(self.phiJ,np.subtract(1,Fj)))
       
        jgridInf = nn_seeds[:]
        VhInf = np.add(np.multiply(jgridInf,np.multiply(1-self.phiJ,np.subtract(1,Fj))),np.multiply((1-(1-self.muF)*self.phiFh)/(self.phiV*(1-self.muV)),self.Fhgrid0))
        VdInf = np.add(np.multiply(jgridInf,np.multiply(1-self.phiJ,Fj)),np.multiply((1-self.phiFd*(1-self.muF))/(self.phiV*(1-self.muV)),self.Fdgrid0))

        #these are the values if you don't actually know Fh
        Fhtemp = np.multiply(self.jgrid0,np.multiply(self.phiJ,np.subtract(1,Fj)))
        Fhtemp = np.add(Fhtemp,np.multiply(self.Fhgrid0,np.multiply((1-self.muF)*self.phiFh*(1-Ff),np.subtract(1,Fv))))
        Fhtemp = np.add(Fhtemp,np.multiply(self.Fdgrid0,(1-self.muF)*self.phiFd*self.rec))
        Fhtemp = np.add(Fhtemp,np.multiply(self.Vhgrid0,np.multiply((1-self.muV)*self.phiV,np.subtract(1,Fv))))
        Fhtemp = np.add(Fhtemp,np.multiply(self.Vdgrid0,(1-self.muV)*self.phiV*self.rec))

        Fdtemp = np.multiply(self.jgrid0,np.multiply(self.phiJ,Fj))
        Fdtemp = np.add(Fdtemp,np.multiply(self.Fhgrid0,np.multiply((1-self.muF)*self.phiFh,np.add(Ff,Fv))))
        Fdtemp = np.add(Fdtemp,np.multiply(self.Fdgrid0,(1-self.muF)*self.phiFd*(1-self.rec)))
        Fdtemp = np.add(Fdtemp,np.multiply(self.Vhgrid0,np.multiply((1-self.muV)*self.phiV,Fv)))
        Fdtemp = np.add(Fdtemp,np.multiply(self.Vdgrid0,(1-self.muV)*self.phiV*(1-self.rec)))
 
        #out = np.multiply((1-self.phiFd*(1-self.muF))/(self.phiV*(1-self.muV)),self.Fdgrid0)
        #out1 =np.multiply((1-(1-self.muF)*self.phiFh)/(self.phiV*(1-self.muV)),self.Fhgrid0)
        #for i in range(0,15):
        #    for j in range(0,100):
        #        print out[i][j], self.Vdgrid0[i][j], out1[i][j], self.Vhgrid0[i][j]
        #exit()

        #Fhtemp = np.multiply(self.jgrid0,np.multiply(self.phiJ,np.subtract(1,Fj)))     
        #Fhtemp = np.add(Fhtemp,np.multiply(self.Fhgrid0,np.multiply((1-self.muF)*self.phiFh*(1-Ff),np.subtract(1,Fv))))
        #Fhtemp = np.add(Fhtemp,np.multiply(self.Fdgrid0,(1-self.muF)*self.phiFd*self.rec))       
        #Fhtemp = np.add(Fhtemp,np.multiply(self.Vhgrid0,np.multiply((1-self.muV)*self.phiV,np.subtract(1,Fv))))
        #Fhtemp = np.add(Fhtemp,np.multiply(self.Vdgrid0,(1-self.muV)*self.phiV*self.rec))

        #Fdtemp = np.multiply(self.jgrid0,np.multiply(self.phiJ,Fj))
        #Fdtemp = np.add(Fdtemp,np.multiply(self.Fhgrid0,np.multiply((1-self.muF)*self.phiFh,np.add(Ff,Fv))))
        #Fdtemp = np.add(Fdtemp,np.multiply(self.Fdgrid0,(1-self.muF)*self.phiFd*(1-self.rec)))       
        #Fdtemp = np.add(Fdtemp,np.multiply(self.Vhgrid0,np.multiply((1-self.muV)*self.phiV,Fv)))
        #Fdtemp = np.add(Fdtemp,np.multiply(self.Vdgrid0,(1-self.muV)*self.phiV*(1-self.rec)))
        
        Vhtemp = np.multiply(self.jgrid0,np.multiply(1-self.phiJ,np.subtract(1,Fj)))
        Vhtemp = np.add(Vhtemp,np.multiply(self.Fhgrid0,np.multiply((1-self.muF)*(1-self.phiFh)*(1-Ff),np.subtract(1,Fv))))
        Vhtemp = np.add(Vhtemp,np.multiply(self.Fdgrid0,(1-self.muF)*(1-self.phiFd)*self.rec))       
        Vhtemp = np.add(Vhtemp,np.multiply(self.Vhgrid0,np.multiply((1-self.muV)*(1-self.phiV),np.subtract(1,Fv))))
        Vhtemp = np.add(Vhtemp,np.multiply(self.Vdgrid0,(1-self.muV)*(1-self.phiV)*self.rec))
        
        Vdtemp = np.multiply(self.jgrid0,np.multiply(1-self.phiJ,Fj))
        Vdtemp = np.add(Vdtemp,np.multiply(self.Fhgrid0,np.multiply((1-self.muF)*self.phiFh,np.add(Ff,Fv))))
        Vdtemp = np.add(Vdtemp,np.multiply(self.Fdgrid0,(1-self.muF)*(1-self.phiFd)*(1-self.rec)))       
        Vdtemp = np.add(Vdtemp,np.multiply(self.Vhgrid0,np.multiply((1-self.muV)*(1-self.phiV),Fv)))
        Vdtemp = np.add(Vdtemp,np.multiply(self.Vdgrid0,(1-self.muV)*(1-self.phiV)*(1-self.rec)))

        self.VhProb = Vhtemp[:]
        self.VdProb = Vdtemp[:]
        self.FdProb = Fdtemp[:]
        self.FhProb = Fhtemp[:]

    def nextGen(self,printOut=False,opt=False,printFI=False):

        # first, get force of infection
        Fj = 0 # for junveniles
        Ff = 0 # for flowering
        Fv = 0 # for vegetative

        # get nearest neighbor counts for juveniles, which are weighted by self.migSeeds
        # contributions, others do not
        nn_seeds = np.zeros((15,100))
        stay = 1-self.migSeeds
        # nn_seeds is effective number of juveniles from neighboring cells
        # get component of juveniles from within the current cell
        self.get_N()
        nn_seeds = np.multiply(self.Fhgrid,np.multiply(stay,np.divide(np.multiply(self.b,self.qGrid),(np.add(1.,np.multiply(self.k,self.N))))))

        # get components from neighboring cells -- use focal cell if neighbor does not exist
        for i in range(0,15):
            for j in range(0,100):
                if math.isnan(self.qGridOrig[i,j]):
                    continue
                if i == 0 or math.isnan(self.qGridOrig[i-1,j]):
                    nn_seeds[i,j] += np.multiply(self.Fhgrid[i,j],np.multiply((1-stay)*0.25,np.divide(np.multiply(self.b,self.qGrid[i,j]),np.add(1.,np.multiply(self.k,self.N[i,j])))))
                else:
                    nn_seeds[i,j] += np.multiply(self.Fhgrid[i-1,j],np.multiply((1-stay)*0.25,np.divide(np.multiply(self.b,self.qGrid[i-1,j]),np.add(1.,np.multiply(self.k,self.N[i-1,j])))))
                if i == 14 or math.isnan(self.qGridOrig[i+1,j]):    
                    nn_seeds[i,j] += np.multiply(self.Fhgrid[i,j],np.multiply((1-stay)*0.25,np.divide(np.multiply(self.b,self.qGrid[i,j]),np.add(1.,np.multiply(self.k,self.N[i,j])))))
                else:
                    nn_seeds[i,j] +=  np.multiply(self.Fhgrid[i+1,j],np.multiply((1-stay)*0.25,np.divide(np.multiply(self.b,self.qGrid[i+1,j]),np.add(1.,np.multiply(self.k,self.N[i+1,j])))))
                if j == 0 or math.isnan(self.qGridOrig[i,j-1]):
                    nn_seeds[i,j] +=  np.multiply(self.Fhgrid[i,j],np.multiply((1-stay)*0.25,np.divide(np.multiply(self.b,self.qGrid[i,j]),np.add(1.,np.multiply(self.k,self.N[i,j])))))
                else:
                    nn_seeds[i,j] +=  np.multiply(self.Fhgrid[i,j-1],np.multiply((1-stay)*0.25,np.divide(np.multiply(self.b,self.qGrid[i,j-1]),np.add(1.,np.multiply(self.k,self.N[i,j-1])))))
                if j == 99 or math.isnan(self.qGridOrig[i,j+1]):
                    nn_seeds[i,j] +=  np.multiply(self.Fhgrid[i,j],np.multiply((1-stay)*0.25,np.divide(np.multiply(self.b,self.qGrid[i,j]),np.add(1.,np.multiply(self.k,self.N[i,j])))))
                else:
                    nn_seeds[i,j] +=  np.multiply(self.Fhgrid[i,j+1],np.multiply((1-stay)*0.25,np.divide(np.multiply(self.b,self.qGrid[i,j+1]),np.add(1.,np.multiply(self.k,self.N[i,j+1])))))

        # Force of infection gets weighted by self.migSpores
        stay = 1.-self.migSpores
        nn_spores = np.multiply(stay,self.Fdgrid)
      
        # get components from neighboring cells -- use focal cell in neighbor does not exist
        for i in range(0,15):
            for j in range(0,100):
                if i == 0 or math.isnan(self.qGridOrig[i-1,j]):
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid[i,j])
                else:
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid[i-1,j])
                if i == 14 or math.isnan(self.qGridOrig[i+1,j]):    
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid[i,j])
                else:
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid[i+1,j])
                if j == 0 or math.isnan(self.qGridOrig[i,j-1]):
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid[i,j])
                else:
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid[i,j-1])
                if j == 99 or math.isnan(self.qGridOrig[i,j+1]):
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid[i,j])
                else:
                    nn_spores[i,j] += np.multiply((1-stay)*0.25,self.Fdgrid[i,j+1])

        # calc Ff, which is vector-borne and depends only on whole grid frequency of flowering diseased
        Ff = 0.
        if np.nansum(self.Fdgrid)+np.nansum(self.Fhgrid) != 0.:
            Ff = np.minimum((self.betaF*np.nansum(self.Fdgrid))/(np.nansum(self.Fdgrid)+np.nansum(self.Fhgrid)),1)

        # now Fj, which is diff for each cell in grid
        Fj = np.minimum(np.multiply(self.betaJ,nn_spores),np.ones((15,100)))

        # Fv, which is just Fj scaled
        Fv = np.minimum(np.multiply(self.betaV,nn_spores),np.ones((15,100)))       
        
        if printFI:
            print ("#", np.nansum(Ff), np.nansum(Fj), np.nansum(Fv))

        if printOut:
            print( (np.sum(self.Fdgrid)+0.)/np.sum(np.add(self.Fhgrid,self.Fdgrid)), np.sum(self.Fhgrid), np.sum(self.Fdgrid))

        # next get expectation of count in each life/infection stage
        # these each just correspond to the equations in Bruns 2017 (J. Ecology)
        Fhtemp = np.multiply(self.jgrid,np.multiply(self.phiJ,np.subtract(1,Fj)))     
        #Fhmean= np.nanmean(Fhtemp,axis=0)

        #for thing in Fhmean:
        #    print thing,
        #print
        
        Fhtemp = np.add(Fhtemp,np.multiply(self.Fhgrid,np.multiply((1-self.muF)*self.phiFh*(1-Ff),np.subtract(1,Fv))))
        Fhtemp = np.add(Fhtemp,np.multiply(self.Fdgrid,(1-self.muF)*self.phiFd*self.rec))       
        Fhtemp = np.add(Fhtemp,np.multiply(self.Vhgrid,np.multiply((1-self.muV)*self.phiV,np.subtract(1,Fv))))
        Fhtemp = np.add(Fhtemp,np.multiply(self.Vdgrid,(1-self.muV)*self.phiV*self.rec))


        Fdtemp = np.multiply(self.jgrid,np.multiply(self.phiJ,Fj))
        Fdtemp = np.add(Fdtemp,np.multiply(self.Fhgrid,np.multiply((1-self.muF)*self.phiFh,np.add(Ff,Fv))))
        Fdtemp = np.add(Fdtemp,np.multiply(self.Fdgrid,(1-self.muF)*self.phiFd*(1-self.rec)))       
        Fdtemp = np.add(Fdtemp,np.multiply(self.Vhgrid,np.multiply((1-self.muV)*self.phiV,Fv)))
        Fdtemp = np.add(Fdtemp,np.multiply(self.Vdgrid,(1-self.muV)*self.phiV*(1-self.rec)))
        
        Vhtemp = np.multiply(self.jgrid,np.multiply(1-self.phiJ,np.subtract(1,Fj)))
        Vhtemp = np.add(Vhtemp,np.multiply(self.Fhgrid,np.multiply((1-self.muF)*(1-self.phiFh)*(1-Ff),np.subtract(1,Fv))))
        Vhtemp = np.add(Vhtemp,np.multiply(self.Fdgrid,(1-self.muF)*(1-self.phiFd)*self.rec))       
        Vhtemp = np.add(Vhtemp,np.multiply(self.Vhgrid,np.multiply((1-self.muV)*(1-self.phiV),np.subtract(1,Fv))))
        Vhtemp = np.add(Vhtemp,np.multiply(self.Vdgrid,(1-self.muV)*(1-self.phiV)*self.rec))
        
        Vdtemp = np.multiply(self.jgrid,np.multiply(1-self.phiJ,Fj))
        Vdtemp = np.add(Vdtemp,np.multiply(self.Fhgrid,np.multiply((1-self.muF)*self.phiFh,np.add(Ff,Fv))))
        Vdtemp = np.add(Vdtemp,np.multiply(self.Fdgrid,(1-self.muF)*(1-self.phiFd)*(1-self.rec)))       
        Vdtemp = np.add(Vdtemp,np.multiply(self.Vhgrid,np.multiply((1-self.muV)*(1-self.phiV),Fv)))
        Vdtemp = np.add(Vdtemp,np.multiply(self.Vdgrid,(1-self.muV)*(1-self.phiV)*(1-self.rec)))

        #self.compute_expectations(1.-self.migSeeds,1.-self.migSpores,self.betaF,0.0001,0.0001)
        def mLL(x):

            stay0 = x[0]
            stay1 = x[1]
            
            bF = x[2]
            bV = x[3]
            bJ = x[4]
            
            if stay0 > 1 or stay1 > 1 or stay0 < 0 or stay1 < 0 or bF < 0 or bV < 0 or bJ < 0:
                return 1000000           
            
            self.compute_expectations(stay0,stay1,bF,bV,bJ)


            Fh = np.nansum(poisson.logpmf(self.Fhgrid,self.FhProb))          
            Fd = np.nansum(poisson.logpmf(self.Fdgrid,self.FdProb))          
            Vh = np.nansum(poisson.logpmf(self.Vhgrid,self.VhProb))          
            Vd = np.nansum(poisson.logpmf(self.Vdgrid,self.VdProb))          
            
            return -(Fh+Fd+Vh+Vd)
            #return -(Fh+Fd)
        

        if not self.determ:
            Fhgrid = poisson.rvs(Fhtemp)
            Vhgrid = poisson.rvs(Vhtemp)
            Fdgrid = poisson.rvs(Fdtemp)
            Vdgrid = poisson.rvs(Vdtemp)

            #nn_seeds_mean = np.nanmean(nn_seeds,axis=0)
            #for thing in nn_seeds_mean:
            #    print thing,
            #print
            jgrid = poisson.rvs(np.where(np.isnan(nn_seeds),0, nn_seeds))
            
            self.Fhgrid0 = self.Fhgrid[:]    
            self.Fdgrid0 = self.Fdgrid[:]    
            self.Vhgrid0 = self.Vhgrid[:]
            self.Vdgrid0 = self.Vdgrid[:]
            self.jgrid0 = self.jgrid[:]

            #print np.mean(Fhtemp),

            self.Fhgrid= Fhgrid[:]
            self.Fdgrid= Fdgrid[:]
            self.Vhgrid= Vhgrid[:]
            self.Vdgrid= Vdgrid[:]
            self.jgrid= jgrid[:]

            self.jgrid = self.jgrid
            self.Fhgrid = self.Fhgrid
            self.Fdgrid =  self.Fdgrid
            self.Vhgrid = self.Vhgrid
            self.Vdgrid = self.Vdgrid
         
            #for i in range(0,15):
            #    for j in range(0,100):
            #        if math.isnan(self.qGridOrig[i][j]):
            #            print self.qGrid[i][j], self.Fhgrid[i][j]
            #            self.Fhgrid[i][j] = 0
            #            self.Vhgrid[i][j] = 0
            #            self.Vdgrid[i][j] = 0
            #            self.Fdgrid[i][j] = 0
            #            self.jgrid[i][j] = 0
   
            if opt:
                #for i in range(0,15):
                #    for j in range(0,100):
                #        if self.Fhgrid[i][j] > 0 and self.Fhgrid0[i][j] == 0:
                #            print self.Fhgrid[i][j]-self.Fhgrid0[i][j], i, j
                #exit()
                xT = [1.-self.migSeeds,1.-self.migSpores,self.betaF,self.betaV,self.betaJ]
                #x0 = xT[:]
                #print x0, xT
                b0 = 0.7
                b1 = 0.7
                b2 = 0.5
                b3 = 0.1
                b4 = 0.5
 
                x0 = [np.random.random()*(1-b0)+b0,np.random.random()*(1-b1)+b1,np.random.random()*b2,np.random.random()*b3,np.random.random()*b4]
                myOpt = minimize(mLL, x0, method='Nelder-Mead', tol=1e-6)# ( bounds=((b0,1),(b1,1),(b2,0.5),(b3,0.1),(b4,0.5)))

                for thing in myOpt.x:
                    sys.stdout.write(thing + ' ')
                for thing in xT:
                    sys.stdout.write(thing+' ')
                print (np.sum(self.Fhgrid),np.sum(self.Fdgrid),np.sum(self.Vhgrid),np.sum(self.Vdgrid), mLL(xT), mLL(myOpt.x)) #,end='')

 
        if self.determ:
            self.Fhgrid= Fhtemp[:]
            self.Fdgrid= Fdtemp[:]
            self.Vhgrid= Vhtemp[:]
            self.Vdgrid= Vdtemp[:]
            self.jgrid= np.multiply(nn_seeds,self.Fhgrid)
        
        #self.Fhgrid = np.maximum( self.Fhgrid , np.zeros((15,100)))
        #self.Vhgrid = np.maximum( self.Vhgrid , np.zeros((15,100)))
        #self.Fdgrid = np.maximum( self.Fdgrid , np.zeros((15,100)))
        #self.Vdgrid = np.maximum( self.Vdgrid , np.zeros((15,100)))

        #Fhprob =  np.sum(poisson.logpmf(self.Fhgrid,Fhtemp))
        #Fdprob =  np.sum(poisson.logpmf(self.Fdgrid,Fdtemp)) 
        #Vhprob =  np.sum(poisson.logpmf(self.Vhgrid,Vhtemp)) 
        #Vdprob =  np.sum(poisson.logpmf(self.Vdgrid,Vdtemp)) 
        
        #FhprobH =  np.sum(poisson.logpmf(self.Fhgrid,self.FhProb))
        #FdprobH =  np.sum(poisson.logpmf(self.Fdgrid,self.FdProb)) 
        #VhprobH =  np.sum(poisson.logpmf(self.Vhgrid,self.VhProb)) 
        #VdprobH =  np.sum(poisson.logpmf(self.Vdgrid,self.VdProb)) 
     

        #mLLloc = -( Fhprob+ Fdprob+ Vhprob+ Vdprob)
        #mLLH = -( FhprobH+ FdprobH+ VhprobH+ VdprobH)

        #print np.sum(self.Fhgrid), np.sum(self.Fdgrid), np.sum(self.Vhgrid), np.sum(self.Vdgrid), mLL, mLLH, mLL-mLLH

    def simTransect(self, burnGens = 50, maxGens = 100, printData = False, printQual = False, printLinTrans = False, nInf= 400):
        if self.readQ:
            self.readQual(printq=False,trim=False)
        self.initPop()
   
        # burn in
        for i in range(0,burnGens):
            self.nextGen(printOut=False)

        # introduce disease
        num = 0
        while num < nInf:
            irand = np.random.randint(5)
            jrand = np.random.randint(100)
            self.Fdgrid[irand,jrand] += 1
            num += 1 
        
        if printData:
            self.printSimData()
            return
        
        if printQual:
            self.printQualData()
            return

        #self.nextGen(printOut=False,opt=False,printFI=True)
        #for i in range(0,20):
        #    self.nextGen(printOut=False,opt=False)
        #if np.nansum(self.Vdgrid) == 0:
        #    print >> sys.stderr, "No output"
        #    exit()
        prev = []
        abun_Fh = [] 
        abun_Fd = [] 
        for i in range(0,maxGens):
            self.nextGen(printOut=False,opt=False)
            Fha = np.nansum(self.Fhgrid)
            Fda = np.nansum(self.Fdgrid)
            Vha = np.nansum(self.Vhgrid)
            Vda = np.nansum(self.Vdgrid)
            ja = np.nansum(self.jgrid)
        
            #print (Fha, Fda, Vha, Vda, ja, Fda/(Fha+Fda+0.))
            
            abun_Fh.append(Fha)
            abun_Fd.append(Fda)
            if Fha + Fda > 0:
                prev.append(Fda/(Fha+Fda+0.))
            else:
                
                prev.append(0.)
            
            if printLinTrans and i == (maxGens - 1):
                Fha = np.nanmean(self.Fhgrid,axis=0)
                Fda = np.nanmean(self.Fdgrid,axis=0)
                Vha = np.nanmean(self.Vhgrid,axis=0)
                Vda = np.nanmean(self.Vdgrid,axis=0)
                ja = np.nanmean(self.jgrid,axis=0)
                for l in range(0,len(Fha)):
                    sys.stdout.write(str(Fha[l])+' ') 
                for l in range(0,len(Fda)):
                    sys.stdout.write(str(Fda[l])+' ') 
                print()
            if not printLinTrans and (np.nansum(self.Fdgrid) + np.nansum(self.Vdgrid) == 0):
                print (np.mean(abun_Fh), np.max(abun_Fh), np.min(abun_Fh), np.mean(abun_Fd), np.max(abun_Fd), np.min(abun_Fd), np.mean(prev), np.max(prev), np.min(prev), i,end=' ')
                if np.nansum(self.Fdgrid) + np.nansum(self.Vdgrid) + np.nansum(self.Fhgrid) + np.nansum(self.Vhgrid) + np.nansum(self.jgrid) > 0:
                    print (1)
                else:
                    print (0) 
                return
        if not printLinTrans:    
            print (np.mean(abun_Fh), np.max(abun_Fh), np.min(abun_Fh), np.mean(abun_Fd), np.max(abun_Fd), np.min(abun_Fd), np.mean(prev), np.max(prev), np.min(prev), i,end=' ')
            if np.nansum(self.Fdgrid) + np.nansum(self.Vdgrid) + np.nansum(self.Fhgrid) + np.nansum(self.Vhgrid) + np.nansum(self.jgrid) > 0:
                print (1)
            else:
                print (0) 

    def simTransectRepGen(self, burnGens = 50, maxGens = 1000):
        if self.readQ:
            self.readQual(printq=False,trim=False)
 
        self.initPop()
   
        # burn in
        for i in range(0,burnGens):
            self.nextGen(printOut=False)

        # introduce disease
        num = 0
        while num < 1:
            irand = np.random.randint(5)
            jrand = np.random.randint(100)
            self.Fdgrid[irand,jrand] += 1
            num += 1 
        
        #self.nextGen(printOut=False,opt=False,printFI=True)
        #for i in range(0,20):
        #    self.nextGen(printOut=False,opt=False)
        #if np.nansum(self.Vdgrid) == 0:
        #    print >> sys.stderr, "No output"
        #    exit()
        prev = []
        abun_Fh = [] 
        abun_Fd = [] 
        for i in range(0,maxGens):
            self.nextGen(printOut=False,opt=False)
            Fha = np.nansum(self.Fhgrid)
            Fda = np.nansum(self.Fdgrid)
            Vha = np.nansum(self.Vhgrid)
            Vda = np.nansum(self.Vdgrid)
            ja = np.nansum(self.jgrid)
        
            if Fha+Fda+0. > 0: 
                #print (Fha+Fda, Fha, Fda, Vha, Vda, ja, Fda/(Fha+Fda+0.),end=' ')
                print ( Fha, Fda, Vha, Vda, ja, Fda/(Fha+Fda+0.),end=' ')
            else:            
                #print (Fha+Fda, Fha, Fda, Vha, Vda, ja, 'NaN',end=' ')
                print ( Fha, Fda, Vha, Vda, ja, 'NaN',end=' ')
            abun_Fh.append(Fha)
            abun_Fd.append(Fda)
            if Fha + Fda > 0:
                prev.append(Fda/(Fha+Fda+0.))
            else:
                prev.append(0.)

            Fha = np.nanmean(self.Fhgrid,axis=0)
            Fda = np.nanmean(self.Fdgrid,axis=0)
            Vha = np.nanmean(self.Vhgrid,axis=0)
            Vda = np.nanmean(self.Vdgrid,axis=0)
            ja = np.nanmean(self.jgrid,axis=0)
            for l in range(0,len(Fha)):
                print (Fha[l],end=' ') 
            for l in range(0,len(Fda)):
                print (Fda[l],end=' ') 
            print()
            if np.nansum(self.Fdgrid) + np.nansum(self.Vdgrid) == 0:
                break
            #    
            #    print np.mean(abun_Fh), np.max(abun_Fh), np.min(abun_Fh), np.mean(abun_Fd), np.max(abun_Fd), np.min(abun_Fd), np.mean(prev), np.max(prev), np.min(prev), i,
            #    if np.nansum(self.Fdgrid) + np.nansum(self.Vdgrid) + np.nansum(self.Fhgrid) + np.nansum(self.Vhgrid) + np.nansum(self.jgrid) > 0:
            #        print 1
            #    else:
            #        print 0 
            #    return
            
        #print np.mean(abun_Fh), np.max(abun_Fh), np.min(abun_Fh), np.mean(abun_Fd), np.max(abun_Fd), np.min(abun_Fd), np.mean(prev), np.max(prev), np.min(prev), i,
        #if np.nansum(self.Fdgrid) + np.nansum(self.Vdgrid) + np.nansum(self.Fhgrid) + np.nansum(self.Vhgrid) + np.nansum(self.jgrid) > 0:
        #    print 1
        #else:
        #    print 0 

    def simSlopeInt(self, burnGens = 50, maxGens = 1000):
        if self.readQ:
            self.readQual(printq=False,trim=False)
        self.initPop()
   
        # burn in
        for i in range(0,burnGens):
            self.nextGen(printOut=False)

        # introduce disease
        num = 0
        while num < 100:
            irand = np.random.randint(5)
            jrand = np.random.randint(100)
            self.Fdgrid[irand,jrand] += 1
            num += 1 
        
        t_lag = 100 #np.random.randint(0,25) 

        #years = [t_lag,t_lag+2,t_lag+9,t_lag+10,t_lag+13]
        years = [t_lag+13]
        for i in range(0,t_lag+14):
            self.nextGen(printOut=False,opt=False)
 
            if i in years:
                abun_binned = np.zeros(12)
                prev_binned = np.zeros(12)
                for k in range(0,15):
                    for j in range(0,100):
                        abun = self.Fhgrid[k][j]+self.Fdgrid[k][j]+self.Vhgrid[k][j]+self.Vdgrid[k][j]
                        if np.isnan(abun) or abun == 0:
                            continue
                        prev = self.Fdgrid[k][j]
                        if abun > 11:
                            abun   = 11
                        abun_binned[abun] += self.Fdgrid[k][j]+self.Fhgrid[k][j]+self.Vdgrid[k][j]+self.Vhgrid[k][j]
                        prev_binned[abun] += self.Fdgrid[k][j]+self.Vdgrid[k][j]
                print (i+2005-t_lag, self.betaF, self.betaV, self.betaJ, t_lag,end=' ')
                for k in range(1,12):
                    if abun_binned[k] > 0:
                        prev_binned[k] /= (abun_binned[k]+0.)
                        print (prev_binned[k],end=' ')                               
                    else:
                        print ('NA', end=' ')
                print()

    def printSimData(self):

        print ("#", np.nansum(self.Fhgrid),  np.nansum(self.Fdgrid), np.nansum(self.Fdgrid)/(np.nansum(self.Fdgrid)+np.nansum(self.Fhgrid)))
        for i in range(0,15):
            for j in range(0,100):
                if math.isnan(self.qGridOrig[i,j]):
                    continue
                    #print i,j,'NA'
                else:
                    print (i,j,self.Fhgrid[i][j]+self.Fdgrid[i][j])
        return 
        
    def printQualData(self):

        print ("#", np.nansum(self.Fhgrid),  np.nansum(self.Fdgrid), np.nansum(self.Fdgrid)/(np.nansum(self.Fdgrid)+np.nansum(self.Fhgrid)))
        for i in range(0,15):
            for j in range(0,100):
                if math.isnan(self.qGridOrig[i,j]):
                    continue
                print (i,j,self.qGridOrig[i][j])
        return      
  
    def printRealData(self,year = False):
        self.readQual(printq=False,trim=False)
        self.readRealData()
        
        if not year:
            print ("#", np.nansum(self.Fh05),  np.nansum(self.Fd05), np.nansum(self.Fd05)/(np.nansum(self.Fd05)+np.nansum(self.Fh05)))
            for i in range(0,15):
                for j in range(0,100):
                    print (self.Fh05[i][j],end=' ')
                print()  
        
            print ("#", np.nansum(self.Fh07),  np.nansum(self.Fd07),np.nansum(self.Fd07)/(np.nansum(self.Fd07)+np.nansum(self.Fh07)))
            for i in range(0,15):
                for j in range(0,100):
                    print(self.Fh07[i][j],end=' ')
                print()  
        
            print("#", np.nansum(self.Fh14),  np.nansum(self.Fd14),np.nansum(self.Fd14)/(np.nansum(self.Fd14)+np.nansum(self.Fh14)))
            for i in range(0,15):
                for j in range(0,100):
                    print (self.Fh14[i][j],end=' ')
                print()  

            print ("#", np.nansum(self.Fh15),  np.nansum(self.Fd15),np.nansum(self.Fd15)/(np.nansum(self.Fd15)+np.nansum(self.Fh15)))
            for i in range(0,15):
                for j in range(0,100):
                    print (self.Fh15[i][j],end=' ')
                print()  
        
            print ("#", np.nansum(self.Fh18),  np.nansum(self.Fd18),np.nansum(self.Fd18)/(np.nansum(self.Fd18)+np.nansum(self.Fh18)))
            for i in range(0,15):
                for j in range(0,100):
                    print (self.Fh18[i][j],end=' ')
                print()  
            return        

            print("#", np.nansum(self.Fh19),  np.nansum(self.Fd19),np.nansum(self.Fd19)/(np.nansum(self.Fd19)+np.nansum(self.Fh19)))
            for i in range(0,15):
                for j in range(0,100):
                    print (self.Fh19[i][j],end=' ')
                print()
            return
        
        else:
            print ("#", np.nansum(self.Fh05),  np.nansum(self.Fd05),np.nansum(self.Fd05)/(np.nansum(self.Fd05)+np.nansum(self.Fh05)))
            for i in range(0,15):
                for j in range(0,100):
                    print (i,j,self.Fh05[i][j]+self.Fd05[i][j])
                print()

    def maxLReal(self):
    
        self.readQual(printq=False,trim=False)
        self.readRealData()

        self.Fhgrid0 = self.Fh14[:]
        self.Fdgrid0 = self.Fd14[:]

        self.Fhgrid = self.Fh15[:]
        self.Fdgrid = self.Fd15[:]

        #for i in range(0,15):
        #    for j in range(0,100):
        #        if self.Fhgrid[i][j] > 0 and self.Fhgrid0[i][j] == 0:
        #            print self.Fhgrid[i][j]-self.Fhgrid0[i][j], i, j
        #exit()

        #print np.nansum(self.F`hgrid),np.nansum(self.Fdgrid),np.nansum(self.Fhgrid0),np.nansum(self.Fhgrid0)

        def mLL(x):


            bF = x[0]
            bV = x[1]
            bJ = x[2]

            if bF < 0 or bV < 0 or bJ < 0:
                return 1000000

            self.compute_expectations_Real(1.-self.migSeeds,1.-self.migSpores,bF,bV,bJ)

            #for val in self.FhProb:
            #    for i in val:
            #        print i,
            #    print

            #print poisson.logpmf(self.Fhgrid,self.FhProb)

            Fh = np.nansum(poisson.logpmf(self.Fhgrid,self.FhProb))
            Fd = np.nansum(poisson.logpmf(self.Fdgrid,self.FdProb))
            
            #print np.nansum(self.FhProb), np.nansum(self.FdProb), np.nansum(self.Fhgrid), np.nansum(self.Fdgrid), Fh, Fd
            
            #exit()
            return -(Fh+Fd)

        b0 = [2,1,2]
        x0 = [np.random.random()*b0[0],np.random.random()*b0[1],np.random.random()*b0[2]]
        myOpt = minimize(mLL, x0, method='Nelder-Mead', tol=1e-3) #  bounds=((b0,1),(b1,1),(b2,0.5),(b3,0.1),(b4,0.5)))

        for thing in myOpt.x:
            print (thing,end=' ')
        print (mLL(myOpt.x))


