# -*- coding: utf-8 -*-

import blueqat.wq as wq
import numpy as np
np.set_printoptions(edgeitems=70)
np.core.arrayprint._line_width = 500

import time
import sys
from os import path

a = wq.Opt()
# ------ for sa ------------ #
a.Ts  = 5     #default 5
a.R   = 0.95  #default 0.95
a.ite = 1000  #default 1000
# -------------------------- #

f = open('_data/input/txt/10_100.txt')

# -------------------- Parameters --------------------- #
Grids = 10    # number of grids/segments power gen will be divided into
U = int(f.readline().strip().split()[1]) # number of units
demand = int(f.readline().strip().split()[1]) # power demand

# Quadratic Coefficients
A = np.zeros(U)
B = np.zeros(U)
C = np.zeros(U)

# Minimum and Maximum production per unit
Pmin = np.zeros(U)
Pmax = np.zeros(U)
# ------------------------------------------------------ #

for line in f:
  fields = line.strip().split()
  if fields[0] == "A":
    for i in range(0,U):
      A[i] = float(fields[i+1])
  elif fields[0] == "B":
    for i in range(0,U):
      B[i] = float(fields[i+1])
  elif fields[0] == "C":
    for i in range(0,U):
      C[i] = float(fields[i+1])
  elif fields[0] == "Pmin":
    for i in range(0,U):
      Pmin[i] = float(fields[i+1])
  elif fields[0] == "Pmax":
    for i in range(0,U):
      Pmax[i] = float(fields[i+1])

f.close()

N = Grids + 2 # 1 extra grid to allow for max prod and another to allow for a shutdown
slack = 5
scope = U*(N) + 5
a.qubo = np.zeros((scope,scope), dtype=float) #creating QUBO matrix

def hi(i):
  return (Pmax[i]-Pmin[i])/Grids

def prod(i,k):
  return round((Pmin[i]+(k-1)*hi(i)))

def prod_slack(i):
  return -(2**i)

def cost(i,k):
  return A[i] + B[i]*prod(i,k) + C[i]*(prod(i,k)**2)

# --------------- objective function ---------------------
for i in range(0,U):
  for k in range(1,N):
    a.qubo[i*N+k][i*N+k] += cost(i,k)
    for m in range(k+1,N):          # extra func
      a.qubo[i*N+k][i*N+m] += 2*C[i]*prod(i,k)*prod(i,m)
# --------------------------------------------------------

delta_A = 30000
delta_B = 5

if len(sys.argv) == 3:
    delta_A = int(sys.argv[1])
    delta_B = int(sys.argv[2])

additive = 0

# restriction 1
for i in range(0,U):
  for k in range(0,N):
    a.qubo[i*N+k][i*N+k] += delta_A * 1  # 1**2 first factor +
    a.qubo[i*N+k][i*N+k] -= delta_A * 2 # -2*1 third factor -
    for kk in range(k+1,N):    
      a.qubo[i*N+k][i*N+kk] += delta_A * 2 # 2*1*1 second factor +
  additive += delta_A*(1) #1**2 constant added to additive

# restriction 2
for i in range(0,U):
  for k in range(1,N):
    a.qubo[i*N+k][i*N+k] += delta_B * (prod(i,k)**2) # first factor +
    a.qubo[i*N+k][i*N+k] += delta_B* -2 * demand * prod(i,k) # third factor -
    for ii in range(i,U):
      for kk in range(1,N):
        if ii == i and k >= kk: continue        
        a.qubo[i*N+k][ii*N+kk] += delta_B* 2 *prod(i,k) * prod(ii,kk) # second factor +
    for s in range(scope-slack, scope):
      a.qubo[i*N+k][s] += delta_B* 2 * demand * prod_slack(s - scope + slack) # slack variables | s - scope + slack -> range from 0 to 5 (slack variables)

for s in range(scope-slack, scope):
  a.qubo[s][s] += delta_B * (prod_slack(s - scope + slack)**2) # first factor +
  a.qubo[s][s] += delta_B* -2 * demand * prod_slack(s - scope + slack) # third factor -
  for s2 in range(s+1,scope):
    a.qubo[s][s2] += delta_B* 2 * prod_slack(s - scope + slack) * prod_slack(s2 - scope + slack)

additive += delta_B*(demand**2) # constant added to additive

def print_result(result, energy, time, mute):
  format_result = np.array_split(result[:scope-slack], U)
  produced = 0
  pcost = 0
  restrict1 = 0
  restrict2 = 0

  for i in range(0,U):
    temp_p = 0
    temp_c = 0
    count = 0
    for k in range(0,N):
      if format_result[i][k] == 1:
        count += 1
      if k == 0 :
        if mute != 1: print("[" + str(format_result[i][k]) + "|", end ="")
      else:
        if mute != 1: print(" " + str(format_result[i][k]) , end ="")
        if format_result[i][k] == 1:
          temp_p += prod(i,k)
          temp_c += cost(i,k)
    if count == 1:
      produced += temp_p
      pcost += temp_c
    elif count > 1:
      restrict1 += 1

    if produced != demand:
      restrict2 = produced - demand

    if mute != 1:
      if count > 1:
        print("]\t------Broke single unit restriction------")
      elif temp_p == 0 and temp_c == 0:
        print("]\t-----------------------------------------")
      else:
        print("]\tcost: {:10.3f}".format(temp_c) + "\tprod: {:10.3f}".format(temp_p))

  if mute != 1:    
    print("------- Slack -------\t\tTCost: {:9.3f}".format(pcost) + "\tTProd: {:9.3f}".format(produced) + "\tDemand: {:9.3f}".format(demand) + "\tRestrict1: {:d}".format(restrict1) + "\tRestrict2: {:.3f}".format(restrict2))
    slack_extra = 0
    for i in range(scope-slack,scope):
      print( "| {:d} ".format(result[i]), end ="")
      if result[i] == 1:
        slack_extra += 2**(i - scope + slack)
    print(f"|\tEnergy: {energy} \tSlackProd: {produced-slack_extra:.3f}\ntook {time:.3f}s\n")

  #if path.exists("_data/debug/"+str(U)+"_"+str(demand)+".txt") == False:
  #  f = open("_data/debug/"+str(U)+"_"+str(demand)+".txt","w")
  #  f.write("Solver Grids Delta_A Delta_B Pcost Time Restriction1 Restriction2\n")
  #f = open("_data/debug/"+str(U)+"_"+str(demand)+".txt","a")
  #f.write("%d %d %d %d %.2f %.3f %d %.3f\n" % (1, Grids, delta_A, delta_B, pcost, time, restrict1, restrict2))
 # f.close() # first term represents what type of program: 1 for fast sampler s_annealer, 2 for normal sampler s_annealer, 3 for s_quantum annealer, 4 for dwave

num_shots=100
min = -1
min_e = 0
for i in range(0,num_shots):
  timer = time.perf_counter()
  result = a.sa() #sampler="fast"

  if min == -1:
   min = result
   min_e = a.E

  if min_e[-1] > a.E[-1]:
    min = result
    min_e = a.E

  timer = time.perf_counter() - timer 
  print_result(result,a.E[-1],timer,0)

#result = a.sqa()
#for z in range(0,8):
  #print_result(result[z])

print("best result")

print_result(min,min_e[-1],0,0)

print("Additive: " + str(additive))
print("Energy w/ additive: " + str(min_e[-1] + additive))