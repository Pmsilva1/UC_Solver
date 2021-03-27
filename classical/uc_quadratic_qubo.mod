# Size of units and grids
param U;
param N := 10;

# Params
param demand;		# Power demand
param Pmin{1..U};	# min prod of unit
param Pmax{1..U};	# max prod of unit

# Quadratic Coefficients
param A{1..U};
param B{1..U};
param C{1..U};
param prod{i in 1..U, k in 1..N+1} := (Pmin[i] + (k-1)*((Pmax[i] - Pmin[i]) / N ));

# Variables
var x{1..U, 1..N+1} binary;
var wx{1..U, 1..N+1, 2..N+1} binary;

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
	sum{i in 1..U, k in 1..N+1} prod[i,k]*x[i,k] >= demand;

s.t. binary_x{i in 1..U}:
	sum{k in 1..N+1} x[i,k] <= 1;

printf: "-----------debug---------\n" > "debug/quadqubo_debug.txt";
for{i in 1..U}{
	printf: "prod(%d):\t", i >> "debug/quadqubo_debug.txt";
	printf{k in 1..N+1}: "%d\t",prod[i,k] >> "debug/quadqubo_debug.txt";
	printf: "\tPmin: %d \tPmax: %d\n",Pmin[i], Pmax[i] >> "debug/quadqubo_debug.txt";
}

#--------------------------------------------------------------------
solve;
#--------------------------------------------------------------------

# Prints data and solution

printf "-----------solution----------\n\n" > "results/result_quadqubo.txt";

printf "Units: %d   Grids: %d\n", U, N >> "results/result_quadqubo.txt";

printf: "total cost: %.2f\n", 
	sum{i in 1..U, k in 1..N+1} (x[i,k] * (A[i] + B[i]*prod[i,k] + C[i]*prod[i,k]**2) + sum{m in k+1..N+1} 2 * C[i] * prod[i,k] * prod[i,m] * wx[i,k,m] )
>> "results/result_quadqubo.txt";

printf: "produced: %.2f\ndemand: %.2f\n\n", sum{i in 1..U, k in 1..N+1} x[i,k]*prod[i,k], demand >> "results/result_quadqubo.txt";

for{i in 1..U}{
	printf{k in 1..N+1}: (if x[i,k] == 1 then "unit %d: On | cost: %.2f | prod: %.2f [%d]\n" else ""), i, (A[i] + B[i]*prod[i,k] + C[i]*prod[i,k]**2) + sum{m in k+1..N+1} 2 * C[i] * prod[i,k] * prod[i,m] * wx[i,k,m],prod[i,k], k >> "results/result_quadqubo.txt";
}

printf{i in 1..U,k in 1..N+1, m in 2..N+1}: (if wx[i,k,m] != 0 then "wx %d%d%d is %d" else ""), i, k , m, wx[i,k,m] >> "results/result_quadqubo.txt";

end;
