for i in {1..20000}; do python scripts/simRepGen.py >>  simsRepGenRandHigh/sim.$i.txt; done
for i in {1..5000}; do python scripts/simRepGen2.py >>  simsRepGenRandHigh2/sim.$i.txt; done
