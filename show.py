#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

dir = "results/"
label = ["ions-density.txt","electrons-density.txt"]
num = len(label)

# Open density file
file = []*num
Content = []*num
CoList = []*num

# Read data from files
for n in range(num) :

    file[n] = open(dir+label[n],"rt")
    Content[n] = file[n].read()
    CoList[n] = Content[n].split("\n")


# Scan files
i = 0
j = 0

a = []
time = []

for line in CoList :
    # Get time
    if line[:4] == "t = " :
        time.insert(len(time),line[4:])
        a.insert(len(a),i)

    i+=1

size = a[1]-a[0]
s = len(CoList[a[0]+1:a[1]][0].split())

arry = np.zeros((size,s))
i=0

print(CoList[a[20]+1:a[21]])
for line in CoList[a[20]+1:a[21]] :
    aux = line.split()

    print(aux)

    for l in range(s) :
        arry[i][l] = float(aux[l])
    i+=1



plt.plot(arry[:-1,0],arry[:-1,1])

print(arry[:,0])
print(arry[:,1])
plt.show()
