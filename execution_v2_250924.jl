pwd()
readdir()

include("parameters.jl");
include("parameters_calib.jl");
include("main_v2_250922.jl");

params_cal = Initiate_Params_cal(bigT,bigN,bigJ,trim,a,b,debt_to_liability,capital_to_deposit,loan_to_asset);
target_moments = [params_cal.debt_to_liability, params_cal.loan_to_asset, params_cal.capital_to_deposit];   
target_moments  

using LaTeXStrings

Random.seed!(5877);

#=
ϵ = -0.1
E
cL = 0.7 # 1.0 originally
cM = 1.3e-5
cO = 0.2
=#

params = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
params_cal = Initiate_Params_cal(bigT,bigN,bigJ,trim,debt_to_liability,capital_to_deposit,loan_to_asset)

shocks = makeShockSequence(params, vFuncs, 100, 10, 100);
paths_false = Initiate_Paths(shocks, params);
paths_true = Initiate_Paths(shocks, params);

@time eq = VFI(params, 1.04, true, 1000, 1e-2)
@time eqq = stationary_distribution(params, 1.04, eq.vFuncs, eq.vFuncsNew, eq.Iterobj_is, 1000, 1e-4);
@time eqqq = solve_model_given_r2(1.04; Params = params, Regime = true)
@time policyyy = Get_PolicyFuncs(params, eq.Iterobj_is, 1.04); # get policy functions from the solution
@time sol_true = solve_model(params, true, 1.0+eps(), 2-eps());
@show policyyy.lPolicy[:,:,3]

eq.vFuncs.qBond
eq.vFuncs.qBond[:,1,3,2,:]

@show params.lGrid[10]+params.sGrid[10]
@show params.lGrid[19]+params.sGrid[1]



theme(:bright) # theme(:ggplot2)
plot(params.bGrid[7:15]./(params.lGrid[10]+params.sGrid[10]), eq.vFuncs.qBond[1,19,7:15,2,1], label ="low loan-to-safe asset", legend =:bottomleft, legendfontsize = 12, guidefontsize = 12, lw = 3, xlabel = "Bank debt/ total asset", ylabel = "Debt price")
plot!(params.bGrid[7:15]./(params.lGrid[15]+params.sGrid[5]), eq.vFuncs.qBond[19,1,7:15,2,1], label ="high loan-to-safe asset", lw = 3)
plot!(params.bGrid[7:15]./(params.lGrid[10]+params.sGrid[10]),eq.vFuncs.qBond[10,10,7:15,2,3], label ="high λ", lw = 3)
plot!(params.bGrid[7:15]./(params.lGrid[10]+params.sGrid[10]), eq_false.vFuncs.qBond[1,1,1,1,1] .*ones(size(eq.vFuncs.qBond[10,10,7:15,2,3])), label ="benchmark", ls =:dash, lw = 3)
savefig("plot_highres2.pdf")

plot(params.bGrid[7:15]./(params.lGrid[10]+params.sGrid[10]), eq.vFuncs.qBond[1,19,7:15,1,2], label ="low loan-to-safe asset", legend =:bottomleft, legendfontsize = 12, guidefontsize = 12, lw = 3, xlabel = "Bank debt/ total asset", ylabel = "Debt price")
plot!(params.bGrid[7:15]./(params.lGrid[10]+params.sGrid[10]), eq.vFuncs.qBond[1,19,7:15,2,2], label ="low loan-to-safe asset", legend =:bottomleft, legendfontsize = 12, guidefontsize = 12, lw = 3, xlabel = "Bank debt/ total asset", ylabel = "Debt price")

eq.vFuncsNew.X

plot(eq.vFuncs.qBond[10,10,9:13,2,1])
plot!(eq.vFuncs.qBond[10,10,9:13,2,2])
plot!(eq.vFuncs.qBond[10,10,9:13,2,3])
ones(size(eq.vFuncs.qBond[10,10,7:11,2,3]))

plot(eq.vFuncs.qBond[:,1,10,1,2])
plot!(eq.vFuncs.qBond[:,1,10,2,2])
plot!(eq.vFuncs.qBond[:,1,10,3,2])

#=
@time eq = VFI(params, 1.4, false, 1000, 1e-2);
@time eqq = stationary_distribution(params, 1.4, eq.vFuncs, eq.vFuncsNew, eq.Iterobj_is, 1000, 1e-2); # run the stationary distribution algorithm with the given parameters and regime
regime = false
@show solve_model_given_r_single(1+eps())
@show solve_model_given_r_single(2-eps())
=# 

@time sol_true = solve_model(params, true, 1.0+eps(), 2-eps());
@time policy_true = Get_PolicyFuncs(params, sol_true.eqq.Iterobj_is, sol_true.Rl_star); # get policy functions from the solution
@show sol_true.Rl_star
@show sol_true.eq.vFuncs.Γ
@show sum(policy_true.failure)
maximum(policy_true.lPolicy), minimum(policy_true.lPolicy)
aggre_loan_supply2(params, sol_true.eqq.vFuncs, sol_true.eqq.Iterobj_is)
aggre_loan_supply2(params, sol_true.eqq.vFuncs, sol_true.eqq.Iterobj_is)
policy_true.lPolicy[1,1,1], policy_true.lPolicy[1,2,1], policy_true.lPolicy[1,3,1]
mean(policy_true.lPolicy .== policy_false.lPolicy)
@show sol_true.eq.vFuncs.qBond[:,1,1,1,1], sol_true.eq.vFuncs.qBond[:,1,2,1,1], sol_true.eq.vFuncs.qBond[:,1,3,1,1]
plot(sol_true.eq.vFuncs.qBond[:,1,1,1,1])
plot!(sol_true.eq.vFuncs.qBond[:,1,2,1,1])
plot!(sol_true.eq.vFuncs.qBond[:,1,3,1,1])

include("main_v2_250922.jl");
@time eq111 = VFI(params, 1.04, true, 1000, 1e-2)
@time eq222 = VFI(params, 1.04, false, 1000, 1e-2)



endT = Int(params_cal.bigT - 3);
@time sol_false = solve_model(params, false, 1.0+eps(), 2-eps());
@time policy_false = Get_PolicyFuncs(params, sol_false.eqq.Iterobj_is, sol_false.Rl_star); # get policy functions from the solution
@time sim_false = simulate_and_moments(params,params_cal,sol_false.eq.vFuncs,policy_false,sol_false.Rl_star,false)
sol_false.eq.vFuncs.Γ
sol_false.Rl_star
policy_false.lPolicy
policy_false.failure
aggre_loan_supply2(params, sol_false.eqq.vFuncs, sol_false.eqq.Iterobj_is)
solve_model_given_r2(sol_false.Rl_star; Params = params, Regime = false)
sim_false.shocks.deltaIndSim
sim_false.paths.capitalToDeposit[:,:,1]
sim_false.paths.nPrimeSim
sim_false.moments
mean(sim_false.paths.debtToLiability[:, :, trim:endT])


@time simulation_false = Simulate_paths(paths_false, policy_false, sol_false.eqq.vFuncs, params, sol_false.Rl_star, 10)

aggre_loan_supply2(params, sol_false.eqq.vFuncs, sol_false.eqq.Iterobj_is)
aggre_loan_supply2(params, sol_false.eqq.vFuncs, sol_false.eqq.Iterobj_is)
policy_false.lPolicy[1,1,1], policy_false.lPolicy[1,2,1], policy_false.lPolicy[1,3,1]

@show sol_false.Rl_star
@show sol_false.eq.vFuncs.Γ
@show minimum(policy_false.nPrimePolicy)
@show sum(policy_false.failure)
@show policy_false.NAV
@show false_NoFailure = policy_false.nPrimePolicy .* (1 .- policy_false.failure)
@show false_Failure = policy_false.nPrimePolicy .* policy_false.failure
@show false_Failure[:, :, ]
aggre_loan_supply2(params, sol_false.eqq.vFuncs, sol_false.eqq.Iterobj_is)

sum(policy_false.failure[1,:,:,:]), sum(policy_false.failure[2,:,:,:]), sum(policy_false.failure[3,:,:,:])
sum(policy_false.failure[:,1,:,:]), sum(policy_false.failure[:,2,:,:]), sum(policy_false.failure[:,3,:,:])
sum(policy_false.failure[:,:,1,:]), sum(policy_false.failure[:,:,2,:]), sum(policy_false.failure[:,:,3,:]), sum(policy_false.failure[:,:,4,:]), sum(policy_false.failure[:,:,5,:])
sum(policy_false.failure[:,:,:,1]), sum(policy_false.failure[:,:,:,2]), sum(policy_false.failure[:,:,:,3])

mean(policy_false.lPolicy[1,:,:]), mean(policy_false.lPolicy[2,:,:]), mean(policy_false.lPolicy[3,:,:])
params.deltaGrid
params.nGrid
mean(policy_false.lPolicy[:,1,:]), mean(policy_false.lPolicy[:,2,:]), mean(policy_false.lPolicy[:,3,:])
mean(policy_false.lPolicy[:,:,1]), mean(policy_false.lPolicy[:,:,4]), mean(policy_false.lPolicy[:,:,10])
mean(policy_false.sPolicy[1,:,:]), mean(policy_false.sPolicy[2,:,:]), mean(policy_false.sPolicy[3,:,:])
mean(policy_false.sPolicy[:,1,:]), mean(policy_false.sPolicy[:,2,:]), mean(policy_false.sPolicy[:,3,:])
mean(policy_false.bPolicy[1,:,:]), mean(policy_false.bPolicy[2,:,:]), mean(policy_false.bPolicy[3,:,:])
mean(policy_false.bPolicy[:,1,:]), mean(policy_false.bPolicy[:,2,:]), mean(policy_false.bPolicy[:,3,:])

sol_false.eq.vFuncs.Γ[1,1,1], sol_false.eq.vFuncs.Γ[1,2,1], sol_false.eq.vFuncs.Γ[1,3,1]
sol_false.eq.vFuncs.Γ[2,1,1], sol_false.eq.vFuncs.Γ[2,2,1], sol_false.eq.vFuncs.Γ[2,3,1]
sol_false.eq.vFuncs.Γ[3,1,1], sol_false.eq.vFuncs.Γ[3,2,1], sol_false.eq.vFuncs.Γ[3,3,1]
sol_false.eq.vFuncs.Γ[1,1,2], sol_false.eq.vFuncs.Γ[1,2,2], sol_false.eq.vFuncs.Γ[1,3,2]
sol_false.eq.vFuncs.Γ[2,1,2], sol_false.eq.vFuncs.Γ[2,2,2], sol_false.eq.vFuncs.Γ[2,3,2]
sol_false.eq.vFuncs.Γ[3,1,2], sol_false.eq.vFuncs.Γ[3,2,2], sol_false.eq.vFuncs.Γ[3,3,2]

maximum(policy_false.lPolicy), minimum(policy_false.lPolicy)
maximum(policy_false.sPolicy), minimum(policy_false.sPolicy)
maximum(policy_false.bPolicy), minimum(policy_false.bPolicy)


g1 = sum(sol_false.eq.vFuncs.Γ[1,2,:]) 
g2 = sum(sol_false.eq.vFuncs.Γ[2,:,:])
g3 = sum(sol_false.eq.vFuncs.Γ[3,:,:])
sum(sol_false.eq.vFuncs.Γ[1,:,:]) * params.deltaGrid[1] + sum(sol_false.eq.vFuncs.Γ[2,:,:]) * params.deltaGrid[2] + sum(sol_false.eq.vFuncs.Γ[3,:,:]) * params.deltaGrid[3]
loanMarket = solve_model_given_r2(sol_false.Rl_star; Params = params, Regime = false)
sum(sol_false.eq.vFuncs.Γ[:,:,1]), sum(sol_false.eq.vFuncs.Γ[:,:,2])
sum(sol_false.eq.vFuncs.Γ[:,1,:]), sum(sol_false.eq.vFuncs.Γ[:,2,:]), sum(sol_false.eq.vFuncs.Γ[:,3,:])
aggre_loan_supply2(params, sol_false.eqq.vFuncs, sol_false.eqq.Iterobj_is)
##############################################################################


#=

solve_model_given_r_single = Rl -> solve_model_given_r(Rl; Params = params, Regime = regime)



@time A = solve_model_given_r2(1.0 + eps(); Params = params, Regime = false)
@time B = solve_model_given_r2(2.0 - eps(); Params = params, Regime = false)
@show params.ϵ, params.E

=# 

cL = 0.3;
params = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params.cL

@time A_lowCL = solve_simulate_and_moments(params, params_cal, false, 1.0 + eps(), 2.0 - eps());
@show A_lowCL.moments, A_lowCL.Rl_star

cL = 0.7;
params = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params.cL

@time A_highCL = solve_simulate_and_moments(params, params_cal, false, 1.0 + eps(), 2.0 - eps());
@show A_highCL.moments, A_highCL.Rl_star

cL = 0.5; cO = 0.01;
params = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params.cO, params.cL
@time A_lowCO = solve_simulate_and_moments(params, params_cal, false, 1.0 + eps(), 2.0 - eps());
@show A_lowCO.moments, A_lowCO.Rl_star

cL = 0.5; cO = 0.5;
params = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params.cO, params.cL
@time A_highCO = solve_simulate_and_moments(params, params_cal, false, 1.0 + eps(), 2.0 - eps());
@show A_highCO.moments, A_highCO.Rl_star

#=
@show B_lowEpsilon = solve_model_given_r2(2.0 - eps(); Params = params, Regime = false)

ϵ = -0.5; E = 150.0;
params = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,α1,α2,α3,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params.ϵ, params.E
@show A_highE = solve_model_given_r2(1.0 + eps(); Params = params, Regime = false)
@show B_highE = solve_model_given_r2(2.0 - eps(); Params = params, Regime = false)

ϵ = -0.5; E = 110.0;
params = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,α1,α2,α3,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params.ϵ, params.E
@show A_higherE = solve_model_given_r2(1.0 + eps(); Params = params, Regime = false)
@show B_higherE = solve_model_given_r2(2.0 - eps(); Params = params, Regime = false)

ϵ = -0.5; E = 120.0;
params = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,α1,α2,α3,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params.ϵ, params.E
@show A_higherE = solve_model_given_r2(1.0 + eps(); Params = params, Regime = false)
@show B_higherE = solve_model_given_r2(2.0 - eps(); Params = params, Regime = false)

ϵ = -0.5; E = 115.0;
params = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,α1,α2,α3,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params.ϵ, params.E
@show A_higherE = solve_model_given_r2(1.0 + eps(); Params = params, Regime = false)
@show B_higherE = solve_model_given_r2(2.0 - eps(); Params = params, Regime = false)
=#