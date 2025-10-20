pwd()

# for simulation
bigT = 100
bigN = 10
bigJ = 100
trim = 10

# for boundary of root finding
# a = 1.0+eps()
b = 2.0-eps()

# for root finding with initial guess 
a = 1.04+eps()

# target moment 
debt_to_liability = 0.122136
loan_to_asset = 0.713694
capital_to_deposit = 0.0913644
loan_rate = 1.0473  ## 2024년 기준

# initial guess for calibrated parameters
cM = 1.3e-5
cO = 0.2
cL = 0.1 # let operation cost for deposit being very small; 
# cL originally 0.5; with 0.5, optimization does not converge
ϵ = - 1.1 # (KEEP IT FIXED) loan demand elasticity, - 1.1 for April 
E = 200.0 # loan demand shifter
dBar = 1.0
σ = 0.9932 # (KEEP IT FIXED) 
