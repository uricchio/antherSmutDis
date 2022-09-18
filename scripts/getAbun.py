import sys

fh = open(sys.argv[1],'r')

depen = sys.argv[2]



keep = 0
for line in fh:
    data = line.strip().split()
    
    if depen == "freq" and float(data[1]) < 1.0 and float(data[2]) < 5.0:
        keep = 1
    elif depen == "dens" and float(data[2]) < 1.0 and float(data[1]) < 5.0:
        keep = 1
    break

if not keep:
    exit()

data = ''
for line in fh:
   data = line.strip()

data = data.split()

print(data[0],data[1],data[2],data[3],data[4])

