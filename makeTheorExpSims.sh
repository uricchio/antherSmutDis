for i in {1..100}; do python scripts/simLinTrans.py 1 1 > TheorExpSims/freq.1.dens.1.$i.txt; done
for i in {1..100}; do python scripts/simLinTrans.py 0 2.5 > TheorExpSims/freq.0.dens.1.$i.txt; done
for i in {1..100}; do python scripts/simLinTrans.py 6 0 > TheorExpSims/freq.1.dens.0.$i.txt; done

for i in {1..100}; do python scripts/simLinTransQual.py 1 1 > TheorExpSims/qual.freq.1.dens.1.$i.txt; done
for i in {1..100}; do python scripts/simLinTransQual.py 0 2.5 > TheorExpSims/qual.freq.0.dens.1.$i.txt; done
for i in {1..100}; do python scripts/simLinTransQual.py 6 0 > TheorExpSims/qual.freq.1.dens.0.$i.txt; done

for i in {1..100}; do cat TheorExpSims/freq.1.dens.1.$i.txt >> TheorExpSims/All.freq.1.dens.1.txt ; done
for i in {1..100}; do cat TheorExpSims/freq.0.dens.1.$i.txt >> TheorExpSims/All.freq.0.dens.1.txt ; done
for i in {1..100}; do cat TheorExpSims/freq.1.dens.0.$i.txt >> TheorExpSims/All.freq.1.dens.0.txt ; done

for i in {1..100}; do cat TheorExpSims/qual.freq.1.dens.1.$i.txt >> TheorExpSims/All.qual.freq.1.dens.1.txt ; done
for i in {1..100}; do cat TheorExpSims/qual.freq.0.dens.1.$i.txt >> TheorExpSims/All.qual.freq.0.dens.1.txt ; done
for i in {1..100}; do cat TheorExpSims/qual.freq.1.dens.0.$i.txt >> TheorExpSims/All.qual.freq.1.dens.0.txt ; done
