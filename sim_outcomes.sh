
for i in {0.05,0.1,0.2,0.5,1,2,5,10,20}; do python scripts/sim.py 1 $i > simsForPlot/sim.betaF.1.betaVJ.$i.txt; done
for i in {0.05,0.1,0.2,0.5,1,2,5,10,20}; do python scripts/sim.py $i 1 > simsForPlot/sim.betaF.$i.betaVJ.1.txt; done
for i in {0.05,0.1,0.2,0.5,1,2,5,10,20}; do python scripts/sim.py 0 $i > simsForPlot/sim.betaF.0.betaVJ.$i.txt; done
for i in {0.05,0.1,0.2,0.5,1,2,5,10,20}; do python scripts/sim.py $i 0 > simsForPlot/sim.betaF.$i.betaVJ.0.txt; done

python scripts/sim.py 1 10 > simsForPlot/sim.betaF.1.betaVJ.10.txt
python scripts/sim.py 1 20 > simsForPlot/sim.betaF.1.betaVJ.20.txt
python scripts/sim.py 1 5 > simsForPlot/sim.betaF.1.betaVJ.5.txt
python scripts/sim.py 10 0 > simsForPlot/sim.betaF.10.betaVJ.0.txt
python scripts/sim.py 10 1 > simsForPlot/sim.betaF.10.betaVJ.1.txt
python scripts/sim.py 20 0 > simsForPlot/sim.betaF.20.betaVJ.0.txt
python scripts/sim.py 5 1 > simsForPlot/sim.betaF.5.betaVJ.1.txt
python scripts/sim.py 20 1 > simsForPlot/sim.betaF.20.betaVJ.1.txt
