import sys 

if len(sys.argv) != 3:
    print("Wrong arguments [Units] [Demand]")
    sys.exit()

import random
import numpy as np

units = int(sys.argv[1])
demand = int(sys.argv[2])

f = open("_data/input/txt/"+str(units)+"_"+str(demand)+".txt","w")
f2 = open("_data/input/dat/"+str(units)+"_"+str(demand)+".dat","w")

f.write("U " + str(units)+"\nD "+str(demand)) # Units and Demand

f2.write("param U := %d;\n" %units)
f2.write("param demand := %d;\n" %demand)

f.write("\nA") # A quadratic coefficient
f2.write("\nparam : A :=")
for i in range(0,units):
    num = random.randint(10,1000)/10

    f.write(" %.1f" % num)
    f2.write("\n\t%d %.1f" % (i+1, num))

f.write("\nB") # A quadratic coefficient
f2.write(";\n\nparam : B :=")
for i in range(0,units):
    num = random.randint(10,100)/100

    f.write(" %.2f" % num)
    f2.write("\n\t%d %.2f" % (i+1, num))

f.write("\nC") # A quadratic coefficient
f2.write(";\n\nparam : C :=")
for i in range(0,units):
    num = random.randint(10,100)/1000

    f.write(" %.3f" % num)
    f2.write("\n\t%d %.3f" % (i+1, num))

pmax = np.zeros(units)

f.write("\nPmin") # Pmin
f2.write(";\n\nparam : Pmin :=")
for i in range(0,units):
    pmax[i] = demand/units * (random.randint(100, demand/units*50)/100) 
    f.write(" %.2f" % (pmax[i]/(demand/units/2)))
    f2.write("\n\t%d %.2f" % (i+1, pmax[i]/(demand/units/2)))

f.write("\nPmax") # Pmin
f2.write(";\n\nparam : Pmax :=")
for i in range(0,units):
    f.write(" %.2f" % pmax[i])
    f2.write("\n\t%d %.2f" % (i+1, pmax[i]))

f2.write(";\n\nend;")

f.close()
f2.close()