# parameter values

qd = 0.9827
β = 0.9827
Rf = 0.0176
wr = 0.85
α = 0.07
ρ = 0.72
g = 0.17
ξ = 0.0214
cF = 0.0
dBar = 1.0
σ = 0.9932
τC = 0.23
z = 1.0
α1 = 0.3 # DELETE
α2 = 0.3 # DELETE
α3 = 0.4 # DELETE
δL = 10.0
δM = 60.0
δH = 100.0
cM = 1.3e-5
cO = 0.2
cL = 1.0
ϵ = -1.1 # loan demand elasticity 

H = [0.947368421052632 0.0526315789473684 0; 0.0235294117647059 0.929411764705882 0.0470588235294118; 0 0.0945945945945946 0.905405405405405]  # transition matrix: [iDelta, iDeltaPrime]
F = [0.9997 0.0002 0.0001; 0.005 0.99 0.005; 0 0.1 0.9] # transition matrix: [iLambda, iLambdaPrime]
M = 1.0

λL = -0.1428
λM = 0.2052
λH = 0.5

γ = 1.05
ϕ = 1.0


n_start = eps()
n_npts = 5 # 11
n_stop = 1000.0
l_start = 1.3e-5
l_npts = 5 # 11
l_stop = 100.0
s_start = eps()
s_npts = 5 # 11
s_stop = 100.0
b_start = eps()
b_npts = 5 # 11
b_stop = 100.0


# calibrated parameters: 8 parameters
# λH
# F[1,3]
# F[2,3]
# F[3,3]
# cM
# cO
# cL
# ϵ
