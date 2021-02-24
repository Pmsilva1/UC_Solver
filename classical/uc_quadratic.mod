# Size of units and grids
param U;
param N;

# Params
param demand;		# Power demand
#param prod{1..U};	# Production of each unit
param Pmin{1..U};	# min prod of unit
param Pmax{1..U};	# max prod of unit

# Quadratic Coefficients
param A{1..U};
param B{1..U};
param C{1..U};

# Variables
var x{1..U} binary;
var y{1..U, 1..N+1} binary;
var cost{1..U};
var prod{1..U};
var sq_prod{1..U} >= 0;

# Objective Function
minimize obj:
	sum{i in 1..U} cost[i];

# Restrictions
s.t. cprod{i in 1..U}:
	prod[i] = sum{k in 1..N+1} ((Pmin[i]+(k-1)*((Pmax[i] - Pmin[i]) / N))*y[i,k]);

s.t. sqc_{i in 1..U}:
	sq_prod[i] >= Pmin[i]*prod[i]*2 - Pmin[i]*Pmin[i]*(1-x[i]);
s.t. sqc__{i in 1..U}:
	sq_prod[i] >= Pmax[i]*prod[i]*2 - Pmax[i]*Pmax[i]*(1-x[i]);
s.t. sqc___{i in 1..U}:
	sq_prod[i] <= Pmax[i]*prod[i] + Pmin[i]*prod[i] - Pmax[i]*Pmin[i]*(1-x[i]);

s.t. ccost{i in 1..U}:
	cost[i] = (1-x[i])*A[i] + B[i]*prod[i] + C[i]*sq_prod[i];

s.t. vars{i in 1..U}:
	x[i] + sum{k in 1..N+1} y[i,k] = 1;

s.t. cdemand:
	sum{i in 1..U}prod[i] = demand;

#--------------------------------------------------------------------
solve;
#--------------------------------------------------------------------

# Prints data and solution
printf "-----------data----------- \n" > "result.txt";

printf{i in 1..U}: "cost(%d): %d\tA:%.2f  B:%.2f  C:%.3f\n", i, cost[i], A[i], B[i], C[i] >> "result.txt"; 

printf "\n" >> "result.txt";

printf{i in 1..U}: "prod(%d): %g <= %d <= %g \n", i, Pmin[i], prod[i], Pmax[i] >> "result.txt";

printf "\n-----------solution----------\n\n" >> "result.txt";

printf "Units: %d   Grids: %d\n", U, N >> "result.txt";

printf: "total cost: %.2f\n", sum{i in 1..U} cost[i] >> "result.txt";
printf: "produced: %.2f\ndemand: %.2f\n\n", sum{i in 1..U} prod[i]*(1-x[i]), demand >> "result.txt";
printf{i in 1..U}: (if x[i] == 0 then "unit %d: On\n" else ""), i >> "result.txt";

end;
