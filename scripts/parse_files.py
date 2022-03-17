import sys

fh = open(sys.argv[1], 'r')

for line in fh:
    data = line.strip().split()
    
    for thing in data[0:3]:
        print(thing,end=" ")
    print()
    for thing in data[3:]:
        print(thing,end=" ")
    print()
    break

for line in fh:
    print(line.strip())
