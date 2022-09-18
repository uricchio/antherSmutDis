for i in {1..20000}; do python scripts/getAbun.py simsRepGenRandHigh/sim.$i.txt dens >> finalAbun/densLow.txt; done
for i in {1..20000}; do python scripts/getAbun.py simsRepGenRandHigh/sim.$i.txt freq >> finalAbun/freqLow.txt; done
