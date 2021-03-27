import subprocess

for i in range(4000,10001,500):
    if i == 0: continue
    for k in range(0,21,1):
        if k == 0: continue
        subprocess.run(['python3', 'qubo/uc_qubo.py', str(i), str(k)])
        print("(%d,%d) Done" % (i,k))