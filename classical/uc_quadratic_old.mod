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

# Variables
var x{1..U} binary;
var y{1..U, 1..N+1} binary;
var prod{1..U};
var sq_prod{1..U} >= 0;

# Objective Function
minimize obj:
	sum{i in 1..U} ((1-x[i])*A[i] + B[i]*prod[i] + C[i]*sq_prod[i]);

# Restrictions
s.t. cprod{i in 1..U}:
	prod[i] = sum{k in 1..N+1} ((Pmin[i]+(k-1)*((Pmax[i] - Pmin[i]) / N))*y[i,k]);

s.t. sqc_{i in 1..U}:
	sq_prod[i] >= Pmin[i]*prod[i]*2 - Pmin[i]*Pmin[i]*(1-x[i]);
s.t. sqc__{i in 1..U}:
	sq_prod[i] >= Pmax[i]*prod[i]*2 - Pmax[i]*Pmax[i]*(1-x[i]);
s.t. sqc___{i in 1..U}:
	sq_prod[i] <= Pmax[i]*prod[i] + Pmin[i]*prod[i] - Pmax[i]*Pmin[i]*(1-x[i]);

s.t. vars{i in 1..U}:
	x[i] + sum{k in 1..N+1} y[i,k] = 1;

s.t. cdemand:
	sum{i in 1..U}prod[i] = demand;

#--------------------------------------------------------------------
solve;
#--------------------------------------------------------------------

# Prints data and solution

printf "-----------solution----------\n\n" > "results/result_quadratic.txt";

printf "Units: %d   Grids: %d\n", U, N >> "results/result_quadratic.txt";

printf: "total cost: %.2f\n", sum{i in 1..U} ((1-x[i])*A[i] + B[i]*prod[i] + C[i]*sq_prod[i]) >> "results/result_quadratic.txt";
printf: "produced: %.2f\ndemand: %.2f\n\n", sum{i in 1..U} prod[i], demand >> "results/result_quadratic.txt";

for{i in 1..U}{
	printf: (if x[i] == 0 then "unit %d: On | cost: %.2f [actual: %.2f] | prod: %.2f [%d]\n" else ""), i, ((1-x[i])*A[i] + B[i]*prod[i] + C[i]*sq_prod[i]), ((1-x[i])*A[i] + B[i]*prod[i] + C[i]*(prod[i]**2)),prod[i], sum{k in 1..N+1} k*y[i,k] >> "results/result_quadratic.txt";
}
end;
