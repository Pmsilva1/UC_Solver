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

# Objective Function
minimize obj:
	sum{i in 1..U, k in 1..N+1} x[i,k] * (A[i] + B[i]*prod[i,k] + C[i]*prod[i,k]**2);

# Restrictions
s.t. binary_x{i in 1..U}:
	sum{k in 1..N+1} x[i,k] <= 1;

s.t. cdemand:
	sum{i in 1..U, k in 1..N+1} prod[i,k]*x[i,k] = demand;

printf: "-----------debug---------\n" > "debug/lin_quadratic.txt";
for{i in 1..U}{
	printf: "prod(%d):\t", i >> "debug/lin_quadratic.txt";
	printf{k in 1..N+1}: "%.2f\t",prod[i,k] >> "debug/lin_quadratic.txt";
	printf: "\tPmin: %.2f \tPmax: %.2f\n",Pmin[i], Pmax[i] >> "debug/lin_quadratic.txt";
}

#--------------------------------------------------------------------
solve;
#--------------------------------------------------------------------

# Prints data and solution

printf "-----------solution----------\n\n" > "results/result_lin_quadratic.txt";

printf "Units: %d   Grids: %d\n", U, N >> "results/result_lin_quadratic.txt";

printf: "total cost: %.2f\n", sum{i in 1..U, k in 1..N+1} x[i,k] * (A[i] + B[i]*prod[i,k] + C[i]*prod[i,k]**2) >> "results/result_lin_quadratic.txt";
printf: "produced: %.2f\ndemand: %.2f\n\n", sum{i in 1..U, k in 1..N+1} x[i,k]*prod[i,k], demand >> "results/result_lin_quadratic.txt";

for{i in 1..U}{
	printf{k in 1..N+1}: (if x[i,k] == 1 then "unit %d: On | cost: %.2f | prod: %.2f [%d]\n" else ""), i, (A[i] + B[i]*prod[i,k] + C[i]*prod[i,k]**2), prod[i,k], k >> "results/result_lin_quadratic.txt";
}
end;
