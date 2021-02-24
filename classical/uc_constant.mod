# Size of units
param U;

# Params
param demand;		# Power demand
param prod{1..U};	# Production of each unit
param Pmin{1..U};	# min prod of unit
param Pmax{1..U};	# max prod of unit

# Quadratic Coefficients
param A{1..U};
param B{1..U};
param C{1..U};

# Variables
var x{1..U} binary;
var cost{1..U};

# Objective Function
minimize obj:
	sum{i in 1..U} cost[i];

# Restrictions
s.t. ccost{i in 1..U}:
	cost[i] = x[i]*(A[i] + B[i]*prod[i] + C[i]*prod[i]*prod[i]);

s.t. cdemand:
	sum{i in 1..U}(prod[i] * x[i]) = demand;

s.t. production1{i in 1..U}:
	Pmin[i]*x[i] <= prod[i]*x[i];

s.t. production2{i in 1..U}:
	prod[i]*x[i] <= Pmax[i]*x[i];

#--------------------------------------------------------------------
solve;
#--------------------------------------------------------------------

# Prints data and solution
printf "-----------data----------- \n" > "result.txt";

printf{i in 1..U}: "cost(%d): %d\tA:%.2f  B:%.2f  C:%.3f\n", i, cost[i], A[i], B[i], C[i] >> "result.txt"; 

printf "\n" >> "result.txt";

printf{i in 1..U}: "prod(%d): %g <= %d <= %g \n", i, Pmin[i], prod[i], Pmax[i] >> "result.txt";

printf "\n-----------solution----------\n\n" >> "result.txt";

printf "units: %d\n", U >> "result.txt";

printf: "total cost: %.2f\n", sum{i in 1..U} cost[i] >> "result.txt";
printf: "produced: %.2f\ndemand: %.2f\n\n", sum{i in 1..U} prod[i]*x[i], demand >> "result.txt";
printf{i in 1..U}: (if x[i] == 1 then "unit %d: On\n" else ""), i >> "result.txt";

end;
