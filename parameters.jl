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
α1 = 0.3
α2 = 0.3
α3 = 0.4
δL = 10.0
δM = 60.0
δH = 100.0
cM = 1.3e-5
cO = 0.2
cL = 1.0

H = fill(1/3, 3, 3)  # transition matrix
Γ = fill(1/3, 3, 3) # transition matrix

λL = 0.0043
λM = 0.0226
λH = 0.5
γ = 1.05
ϕ = 1.0


n_start = eps()
n_npts = 11
n_stop = 100.0
l_start = 1.3e-5
l_npts = 11
l_stop = 100.0
s_start = eps()
s_npts = 11
s_stop = 100.0
b_start = eps()
b_npts = 11
b_stop = 100.0