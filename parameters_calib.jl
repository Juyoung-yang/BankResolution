# for simulation
bigT = 100
bigN = 10
bigJ = 100
trim = 10

# for boundary of root finding
a = 1.0+eps()
b = 2.0-eps()

# target moment 
debt_to_liability = 0.122136
capital_to_deposit = 0.0913644
loan_to_asset = 0.713694

# initial guess for calibrated parameters
cM = 1.3e-5
cO = 0.2
cL = 0.5
ϵ = -0.5 # loan demand elasticity, 1.1 for April
E = 100.0 # loan demand shifter
dBar = 1.0
