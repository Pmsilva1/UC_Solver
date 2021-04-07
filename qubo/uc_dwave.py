# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import json
from os import path

#QUBO dict
from collections import defaultdict

# Dwave libs
import dwave.system
from dwave.system import DWaveSampler, DWaveCliqueSampler, EmbeddingComposite, FixedEmbeddingComposite
from minorminer import find_embedding
from dwave.embedding.pegasus import find_clique_embedding
import dwave.inspector
import dwave_networkx as dnx
import networkx as nx


# --------------------- File Read --------------------- #
f = open('_data/input/txt/10_100.txt')
print("Reading File Parameters and Data")
# -------------------- Parameters --------------------- #
embedding_kind = 'regular'
#embedding_kind = 'clique'
Grids = 10    # number of grids/segments power gen will be divided into
Slack = 5
U = int(f.readline().strip().split()[1]) # number of units
Demand = int(f.readline().strip().split()[1]) # power demand

# Quadratic Coefficients
A = np.zeros(U)
B = np.zeros(U)
C = np.zeros(U)

# Minimum and Maximum production per unit
Pmin = np.zeros(U)
Pmax = np.zeros(U)
# ------------------------------------------------------ #
print(f"Problem space: {U} units w/ {Grids} grids and {Slack} slack variables (Demand: {Demand})")

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
# --------------------- File Read End --------------------- #

N = Grids + 2 # 1 extra grid to allow for max prod and another 1 to allow for a shutdown
scope = U*(N) + Slack

sampler = DWaveSampler()
#if embedding_kind == "regular":
#  sampler = DWaveSampler()
#elif embedding_kind == "clique":
 # sampler = DWaveCliqueSampler()  
# --------------------------- Find embedding or Open Embedding --------------------------- #
def get_embedding():
  embedding = []
  if path.exists(f'_data/embedding/{embedding_kind}_{U}_{Grids}_{Slack}.json') == False:
    print(f"Finding {embedding_kind} Embedding")
    if embedding_kind == "regular":
      embedding = find_embedding(nx.complete_graph(scope).edges(), sampler.edgelist)
    elif embedding_kind == "clique":
      embedding = find_clique_embedding(scope, target_graph=sampler.to_networkx_graph())
    #print(embedding)

    with open(f'_data/embedding/{embedding_kind}_{U}_{Grids}_{Slack}.json', "w") as ff:
            json.dump(embedding, ff)
  else:
    print(f"Importing saved {embedding_kind} Embedding")
  f = open(f'_data/embedding/{embedding_kind}_{U}_{Grids}_{Slack}.json')
  embedding = json.load(f)

    #dnx.draw_pegasus_embedding(sampler.to_networkx_graph(), embedding, unused_color=None)
    #plt.savefig(f'images/embedding_{embedding_type}N{N}q{q:.2f}B{B}P{P:.3f}.png')

  return embedding
# --------------------------- Find embedding or Open Embedding --------------------------- #

# ---------------------- /| Utility Functions |\ ---------------------- #
def hi(i):
  return (Pmax[i]-Pmin[i])/Grids

def prod(i,k):
  return (Pmin[i]+(k-1)*hi(i))

def prod_slack(i):
  return -(2**i)

def cost(i,k):
  return A[i] + B[i]*prod(i,k) + C[i]*prod(i,k)**2

def print_qubo(Qubo):
  for key, value in Qubo.items():
    print(key, ' : ', value)

def qubo(params):
  print(f"Creating Qubo dict w/ deltaP_A {params['deltaP_A']} & delta_B {params['delta_B']}")
  delta_B = params["delta_B"]
  deltaP_A = params["deltaP_A"]

  Qubo = defaultdict(float)
  additive = 0

  # objective function
  for i in range(0,U):
    for k in range(1,N):
      Qubo[str(i*N+k), str(i*N+k)] += float( cost(i,k) )
      for m in range(k+1,N):          # extra func
        Qubo[str(i*N+k), str(i*N+m)] += float( 2*C[i]*prod(i,k)*prod(i,m) )
   
  # restriction 2
  for i in range(0,U):
    for k in range(1,N):
      Qubo[str(i*N+k), str(i*N+k)] += float( delta_B*( prod(i,k)**2 ) )       # first term: sum xi **2
      Qubo[str(i*N+k), str(i*N+k)] += float( delta_B*( -2*Demand*prod(i,k) ) ) # third term: -2*C*sum xi
      for ii in range(i,U):
        for kk in range(1,N):
          if ii == i and k >= kk: continue        
          Qubo[str(i*N+k), str(ii*N+kk)] += float( delta_B*( 2*prod(i,k)*prod(ii,kk) ) ) # second term: 2*sum xi*xj, where j > i
      for s in range(scope-Slack, scope):
        Qubo[str(i*N+k), str(s)] += float( delta_B*( 2*prod(i,k)*prod_slack(s - scope + Slack) ) )# Slack variables | s - scope + Slack -> range from 0 to 5 (Slack variables)
  
  for s in range(scope-Slack, scope):
    Qubo[str(s), str(s)] += float( delta_B * (prod_slack(s - scope + Slack)**2) )# first factor +
    Qubo[str(s), str(s)] += float( delta_B* -2 * Demand * prod_slack(s - scope + Slack) )# third factor -
    for s2 in range(s+1,scope):
      Qubo[str(s), str(s2)] += float( delta_B* 2 * prod_slack(s - scope + Slack) * prod_slack(s2 - scope + Slack) )
  
  additive += delta_B*( Demand**2 )

  max_bias = 0
  for value in Qubo.values():
    if abs(value) >= max_bias:
      max_bias = value

  delta_A = max_bias * deltaP_A
  params["delta_A"] = delta_A
  max_bias_final = max_bias

  # restriction 1
  for i in range(0,U):
    for k in range(0,N):
      Qubo[str(i*N+k), str(i*N+k)] += float( delta_A*( 1 ) )  # 1**2
      Qubo[str(i*N+k), str(i*N+k)] += float( delta_A*( -2 ) ) # -2*1
      for kk in range(k+1,N):    
        Qubo[str(i*N+k), str(i*N+kk)] += float( delta_A*( 2 ) ) # 2*1*1
    additive += delta_A*(1) #1**2

  for value in Qubo.values():
    if abs(value) >= max_bias_final:
      max_bias_final = value

  #print_qubo(Qubo)

  params["additive"] = float( additive )
  params["max_bias_B"] = float( max_bias )
  params["max_bias_final"] = float( max_bias_final )

  print(f"Created Qubo dict w/ delta_A {params['delta_A']:.2f} & max_bias_B {params['max_bias_B']:.3f} & max_bias_final {params['max_bias_final']:.3f}")
  return Qubo

def print_result(result,params,time,mute,file):
  z = 0
  muted = mute
  for sample, energy, n_occurences, chain_break_freq in response.data():
    #print(sample)
    if muted != 1:
      if z == 10:
        muted = 1
      else:
        if z != 0: print("\n")
        z+=1

    produced = 0
    produced_broken = 0
    pcost = 0
    restrict1 = 0
    restrict2 = 0
    slack_on = 0

    for i in range(0,U):
      temp_p = 0
      temp_c = 0
      count = 0
      for k in range(0,N):
        if sample[str(i*N+k)] == 1:
          count += 1
        if k == 0 :
          if muted != 1: print(f"[{sample[str(i*N+k)]}|", end ="")
        else:
          if muted != 1: print(f" {sample[str(i*N+k)]}", end ="")
          if sample[str(i*N+k)] == 1:
            temp_p += prod(i,k)
            temp_c += cost(i,k)

      if count == 1:
        produced += temp_p
        pcost += temp_c
      else:
        restrict1 += 1
        produced_broken += temp_p
      if produced != Demand:
        restrict2 = produced - Demand
      if muted != 1:
        if count > 1:
          print("]\t------Broke single unit restriction------")
        elif temp_p == 0 and temp_c == 0:
          print("]\t-----------------------------------------")
        else:
          print(f"]\tcost: {temp_c:10.3f}\tprod: {temp_p:10.3f}")

    if muted != 1:    
      print(f"------- Slack -------\t\tTCost: {pcost:9.3f}\tTProd: {produced:9.3f}\tDemand: {Demand:9.3f}\tRestrict1: {restrict1:d}\tRestrict2: {restrict2:.3f}")
      slack_extra = 0
      for i in range(scope-Slack,scope):
        if muted != 1: print(f"| {sample[str(i)]} ", end ="")
        if sample[str(i)] == 1:
          slack_on += 1
          slack_extra += 2**(i - scope + Slack)
      print(f"|\t\t\t\tProd_check: {produced+produced_broken-slack_extra:.3f}")
      print(f"Energy: {energy}\tOcurrences: {n_occurences:4d} \tChain break freq: {chain_break_freq:.3f}")
      #print(f"took {time:.3f}s\n")

    if file == 1:
      if path.exists(f"_data/results/{U}_{Demand}.txt") == False:
        f = open(f"_data/results/{U}_{Demand}.txt","w")
        f.write("Solver Grids Slack DeltaP_A Delta_A Delta_B Anneal_Time Max_Bias_B Max_Bias Chain_Strength_P Chain_Strength Chain_Break Energy Pcost Time Restriction1 Restriction2 Restriction2_abs Restriction2T Slack_On\n")
      f = open(f"_data/results/{U}_{Demand}.txt","a")
      f.write(f"{4} {Grids} {Slack} {params['deltaP_A']} {params['delta_A']} {params['delta_B']} {params['anneal_time']} {params['max_bias_B']} {params['max_bias_final']} {params['chainstrength_P']} {params['chainstrength']} {chain_break_freq} {energy} {pcost} {time} {restrict1} {restrict2} {abs(restrict2)} {restrict2+produced_broken-slack_extra} {slack_on}\n")
      f.close() # first term represents what type of program: 1 for fast sampler s_annealer, 2 for normal sampler s_annealer, 3 for s_quantum annealer, 4 for dwave
# ---------------------- \| Utility Functions |/ ---------------------- #

# ---------------------- /| Main Program |\ ---------------------- #
params = defaultdict(float)
# Problem parameters
run_dwave = 1
params["deltaP_A"] = float( 1 )
params["delta_B"] = float( 1 )
# ------------------
embedding = get_embedding()
Qubo = qubo(params)


if run_dwave != 0:
  # Sampler/Annealing parameters
  params["chainstrength_P"] = float( 0.5 )
  params["chainstrength"] = float( params["chainstrength_P"] * params["max_bias_final"] )# Default 1 | Good rule of thumb is to have the same order of magnitude as Qubo.
  params["numruns"] = 10        # Default 1000
  params["anneal_time"] = 20      # Default 20
  #-----------------------------
  
  print(f"Running D-Wave w/ chainstrength {params['chainstrength_P']:.2f}|{params['chainstrength']:.2f} (%|final) & numruns {params['numruns']:d} & anneal_time {params['anneal_time']:d}")

  composite = FixedEmbeddingComposite(sampler, embedding=embedding)
  response = composite.sample_qubo(Qubo, annealing_time=params["anneal_time"], chain_strength=params["chainstrength"], num_reads=params["numruns"])

  print("Results:\n")
  print_result(response,params,time=0,mute=0,file=1)

  #dwave.inspector.show(response)

# ---------------------- \| Utility Functions |/ ---------------------- #

# minimum possible energy: -283987.43787980796 |\ 173.045734368

# maximum possible energy: -2466607.9822456012  |\  480.1178932

def extra():
  energy = 0
  costx = 0
  vars1 = [1, 6, 7, 8]
  vars2 = [2, 10, 7, 8]

  for i in range(len(vars1)):
    vars1[i] -= 1
    costx += cost(vars1[i],vars1[i])

  for k in range(len(vars2)):
    vars2[k] += 1

  for i in vars1:
      for k in vars2:
        energy += Qubo[str(i*N+k), str(i*N+k)]
        for m in vars2:
          if m < k: continue         
          energy += Qubo[str(i*N+k), str(i*N+m)]
  print(costx)
  print(energy)