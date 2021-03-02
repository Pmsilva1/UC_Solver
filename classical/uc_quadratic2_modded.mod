# Size of units and grids
param U;
param N;

# Params
param demand;		# Power demand
param Pmin{1..U};	# min prod of unit
param Pmax{1..U};	# max prod of unit

# Quadratic Coefficients
param A{1..U};
param B{1..U};
param C{1..U};
param prod{i in 1..U, k in 1..N+1} := (Pmin[i] + sum{j in 1..k}((Pmax[i] - Pmin[i]) / N ));

# Variables
var x{1..U,1..N+1} binary;
var wx{1..U,1..N+1,2..N+1} binary;

# Objective Function
minimize obj:
	sum{i in 1..U, k in 1..N+1} (x[i,k] * (A[i] + B[i]* prod[i,k] + C[i]* prod[i,k]**2)
	+ sum{m in k+1..N+1} 2 * C[i] * prod[i,k] * prod[i,m] * wx[i,k,m] );

# Restrictions
s.t. w1{i in 1..U, k in 1..N+1}:
	sum{m in k+1..N+1} wx[i,k,m] <= x[i,k];

s.t. w2{i in 1..U, m in 2..N+1}:
	sum{k in 1..N+1} wx[i,k,m] <= x[i,m];

s.t. w3{i in 1..U, k in 1..N+1, m in k+1..N+1}:
	wx[i,k,m] >= x[i,k] + x[i,m] - 1;	

s.t. cdemand:
	sum{i in 1..U, k in 1..N+1} prod[i,k]*x[i,k] = demand;

s.t. lock{i in 1..U}:
	sum{k in 1..N+1} x[i,k] <= 1;

printf: "-----------debug---------\n" > "result2_modded_debug.txt";
for{i in 1..U}{
	printf: "prod(%d):\t", i >> "result2_modded_debug.txt";
	printf{k in 1..N+1}: "%d\t",prod[i,k] >> "result2_modded_debug.txt";
	printf: "\tPmin: %d \tPmax: %d\n",Pmin[i], Pmax[i] >> "result2_modded_debug.txt";
}

#--------------------------------------------------------------------
solve;
#--------------------------------------------------------------------

# Prints data and solution

printf "-----------solution----------\n\n" > "result2_modded.txt";

printf "Units: %d   Grids: %d\n", U, N >> "result2_modded.txt";

printf: "total cost: %.2f\n", 
	sum{i in 1..U, k in 1..N+1} (x[i,k] * (A[i] + B[i]*prod[i,k] + C[i]*prod[i,k]**2) + sum{m in k+1..N+1} 2 * C[i] * prod[i,k] * prod[i,m] * wx[i,k,m] )
>> "result2_modded.txt";

printf: "produced: %.2f\ndemand: %.2f\n\n", sum{i in 1..U, k in 1..N+1} x[i,k]*prod[i,k], demand >> "result2_modded.txt";

for{i in 1..U}{
	printf{k in 1..N+1}: (if x[i,k] == 1 then "unit %d: On | %d | [%d]\n" else ""), i, prod[i,k], k >> "result2_modded.txt";
}
end;
