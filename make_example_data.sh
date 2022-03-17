python scripts/printQual.py 1 1 > qualExample.txt
python scripts/printSimExample.py 1 1 > simExample.txt
python scripts/printRealData.py > gridReal.txt

mv gridReal.txt exampleData
mv qualExample.txt exampleData
mv simExample.txt exampleData
