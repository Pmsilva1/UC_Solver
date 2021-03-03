import numpy as np
import pandas as pd

#QUBO dict
from collections import defaultdict

# Dwave libs
from dwave.system import DWaveSampler, FixedEmbeddingComposite
import dwave.inspector
import dwave_networkx as dnx

U = 10 # number of units
Grids = 10 # number of grids
N = Grids + 1 

A = [5,3,5,10,20,35,10,50,100,25]
B = [0.1,0.2,0.4,0.2,0.3,0.2,0.3,0.1,0.2,0.3]
C = [0.01,0.02,0.3,0.2,0.01,0.02,0.03,0.05,0.09,0.05]

demand = 100
Pmin = [4,5,10,6,1,25,50,3,5,20]
Pmax = [25,20,30,20,10,50,100,30,45,60]

Qubo = defaultdict(float)

def hi(i):
  return (Pmax[i]-Pmin[i])/Grids

def prod(i,k):
  return Pmin[i]+(k)*hi(i)

def cost(i,k):
  return A[i] + B[i]*prod(i,k) + C[i]*(prod(i,k)**2)

# objective function
for i in range(0,U):
  for k in range(0,N):
    Qubo[str(i*N+k)][str(i*N+k)] = cost(i,k)

delta_A = 4500
delta_B = 5
additive = 0

# restriction 1
for i in range(0,U):
  for k in range(0,N):
    Qubo[str(i*N+k)][str(i*N+k)] += delta_A*(1)  # 1**2
    Qubo[str(i*N+k)][str(i*N+k)] += delta_A*(-2) # -2*1
    for kk in range(k+1,N):    
      Qubo[str(i*N+k)][str(i*N+kk)] += delta_A*(2) # 2*1*1
  additive += delta_A*(1) #1**2

# restriction 2
for i in range(0,U):
  for k in range(0,N):
    Qubo[str(i*N+k)][str(i*N+k)] += delta_B*(prod(i,k)**2)
    Qubo[str(i*N+k)][str(i*N+k)] += delta_B*(-2*demand*prod(i,k))
    for ii in range(i,U):
      for kk in range(0,N):
        if ii == i and k > kk: continue        
        Qubo[str(i*N+k)][str(ii*N+kk)] += delta_B*(2*prod(i,k)*prod(ii,kk))

additive += delta_B*(demand**2)

print(Qubo)