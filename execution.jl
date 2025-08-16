pwd()
readdir()

include("C:\\Users\\master\\Desktop\\(2025) banking regulation\\BankResolution\\parameters.jl");
include("C:\\Users\\master\\Desktop\\(2025) banking regulation\\BankResolution\\main_v2.jl");

# 0. set the parameters
params = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,α1,α2,α3,δL,δM,δH,cM,cO,cL,ϵ,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)

# 1. solve the model under given params
@time sol = solve_model(params,true, 0.01, 1.9);
@time policy_ex = Get_PolicyFuncs(params, sol.eq.Iterobj_is, sol.Rl_star); # get policy functions from the solution

# 2. generate shock sequenc, simulate paths, and calculate moments
@time shocks = makeShockSequence(params, sol.eq.vFuncs, 100, 100, 10); # make shock sequence for simulation
paths = Initiate_Paths(shocks, params); # initiate paths for simulation

@time paths_ex = Simulate_paths(paths, policy_ex, sol.eq.vFuncs, params, sol.Rl_star, 10); # simulate paths and calculate moments
################################################################################################################################################

plot(paths_ex.failureSim[:,1,1]) # loan policy function for the first bank J at the first observation N
plot!(paths_ex.failureSim[:,1,3])
plot(paths_ex.lSim[:,1,1])
plot!(paths_ex.lSim[:,1,2])
paths_ex.failureSim[:,1,1]
plot(paths_ex.bSim[:,1,1])
paths_ex.lambdaSim[:,1,1]

sum(paths_ex.failureSim[:,10,3])
################################################################################################################################################
@time sol_false = solve_model(params,false, 0.01, 1.9);
@time policy_false = Get_PolicyFuncs(params, sol_false.eq.Iterobj_is, sol_false.Rl_star); # get policy functions from the solution
@time shocks_false = makeShockSequence(params, sol_false.eq.vFuncs, 100, 100, 10);
paths_false_exanti = Initiate_Paths(shocks, params); # initiate paths for simulation
@time paths_false = Simulate_paths(paths_false_exanti, policy_false, sol_false.eq.vFuncs, params, sol_false.Rl_star, 10); # simulate paths and calculate moments

mean(paths_false.lSim[:,1,1]) # mean failure rate for the first bank J at the first observation N
mean(paths_ex.lSim[:,1,1]) # mean failure rate for the first bank J at the first observation N

mean(paths_false.sSim[:,1,1])
mean(paths_ex.sSim[:,1,1])

mean(paths_false.bSim[:,1,1])
mean(paths_ex.bSim[:,1,1])

mean(paths_false.nSim[:,1,1])
mean(paths_ex.nSim[:,1,1])

mean(paths_false.assetPrimeSim[:,1,1])
mean(paths_ex.assetPrimeSim[:,1,1])

################################################################################################################################################
