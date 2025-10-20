pwd()
readdir()

include("parameters.jl");
include("parameters_calib.jl");
include("main_v2_noLoanConstraint.jl");
# include("main_v2_LoanConstraint.jl");



#####################################################################################################################
################# SET THE PARAMETERS ################################################################################ 
#####################################################################################################################


params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
# @show params.nGrid, params.lGrid, params.sGrid, params.bGrid

# b = 1.3+eps()
params_cal = Initiate_Params_cal(bigT,bigN,bigJ,trim,a,b,debt_to_liability,capital_to_deposit,loan_to_asset,loan_rate,cM,cO,cL,ϵ,E,dBar,σ);
@show params_cal.a, params_cal.b
@show params.τC
@show params_cal.ϵ, params_cal.σ

bigT = params_cal.bigT
bigN = params_cal.bigN
bigJ = params_cal.bigJ
trim = params_cal.trim

# shocks = makeShockSequence(params,vFuncs,bigT,bigN,bigJ);
# paths = Initiate_Paths(shocks, params);

Random.seed!(5877);

# initial parameter guess 
# @time Init = solve_simulate_and_moments(params,params_cal,false)


#####################################################################################################################
################# main code: calibration ############################################################################ 
#####################################################################################################################


using Sobol, Distributions, Distributed
@time calibration(params, params_cal, false)

#####################################################################################################################
#################  change parameters one by one #####################################################################

cM = 1.3e-5;
cO = 0.1;
cL = 0.1;
E = 170.0;
dBar = 1.0;

params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments



@time eq = solve_model(params_ex, false, params_cal.a, params_cal.b);

@show params_cal.a
@time eq_a = solve_model_given_r(params_cal.a; Params = params_ex, Regime = false);
@time eq_b = solve_model_given_r2(params_cal.b; Params = params_ex, Regime = false);
sum(eq_b.eq.vFuncs.Γ)


@time solve_simulate_and_moments(params_ex, params_cal, false)





#####################################################################################################################

 @show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
cM = cM * 1.2;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.ϵ, params_ex.σ
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false) 

@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
cM = 1.3e-5;
cL = 0.1;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.ϵ, params_ex.σ
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
@show ex.moments
# debt은 줄어들었으나 여전히 capital 이 높음; E 가 낮은게 loan rate을 낮추는데 큰 역할을 함 

@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
cM = 1.3e-4; # convex cost을 높이기 
cL = 0.1;
cO = 0.2;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.ϵ, params_ex.σ
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
@show ex.moments
@show ex.moments.other_moments
# debt이 줄어들고 capital 도 줄어드는 역할. cM을 계속해서 줄여보자

@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
cM = 1.0e-2; # 계속해서 convex cost을 높이기 
E = 170.0; # loan demand 도 줄여봄 
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.ϵ, params_ex.σ
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
@show ex.moments

@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
cM = 1.0e-2; # 계속해서 convex cost을 높이기 
# cM = 5.0e-2 까지 높이면 fail to converge in optimization
E = 190.0; # loan demand 도 줄여봄 
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.ϵ, params_ex.σ
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
@show ex.moments


# Oct 12 밤 새벽동안 돌리기 
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
cM = 1.0e-2; # 계속해서 convex cost을 높이기 
# cM = 5.0e-2 까지 높이면 fail to converge in optimization
E = 200.0; # loan demand 도 줄여봄 
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.ϵ, params_ex.σ
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
@show ex.moments

@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
cM = 2.0e-2; # 계속해서 convex cost을 높이기 
# cM = 5.0e-2 까지 높이면 fail to converge in optimization
E = 200.0; # loan demand 도 줄여봄 
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.ϵ, params_ex.σ
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
@show ex.moments

@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
cM = 3.0e-2; # 계속해서 convex cost을 높이기 
# cM = 5.0e-2 까지 높이면 fail to converge in optimization
E = 200.0; # loan demand 도 줄여봄 
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.ϵ, params_ex.σ
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
@show ex.moments

@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
cM = 1.0e-2; # 계속해서 convex cost을 높이기 
# cM = 5.0e-2 까지 높이면 fail to converge in optimization
E = 200.0; # loan demand 도 줄여봄 
cO = 0.3;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.ϵ, params_ex.σ
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
@show ex.moments

@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
cM = 1.3e-5;  
# cM = 5.0e-2 까지 높이면 fail to converge in optimization
E = 200.0; # loan demand 도 줄여봄 
cO = 10.0; # operation cost을 키워보기
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.ϵ, params_ex.σ
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
@show ex.moments


@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
cM = 1.3e-5;  
# cM = 5.0e-2 까지 높이면 fail to converge in optimization
E = 200.0; # loan demand 도 줄여봄 
cO = 100.0; # operation cost을 키워보기
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.ϵ, params_ex.σ
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
@show ex.moments


@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
cM = 1.3e-5;  
cL = 0.3;
# cM = 5.0e-2 까지 높이면 fail to converge in optimization
E = 200.0; # loan demand 도 줄여봄 
cO = 100.0; # operation cost을 키워보기
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.ϵ, params_ex.σ
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
@show ex.moments
# higher cL helps controling small capital(n prime), but still no bank failure happening 


@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
cM = 1.0e-2;  
cL = 0.3;
# cM = 5.0e-2 까지 높이면 fail to converge in optimization
E = 200.0; # loan demand 도 줄여봄 
cO = 100.0; # operation cost을 키워보기
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.ϵ, params_ex.σ
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
@show ex.moments
ex.eq


@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
cM = 1.0e-2;  
cL = 0.3;
# cM = 5.0e-2 까지 높이면 fail to converge in optimization
E = 200.0; # loan demand 도 줄여봄 
cO = 100.0; # operation cost을 키워보기
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.ϵ, params_ex.σ
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
@show ex.moments


# solve the model after adjusting psi function 
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
cM = 1.0e-2;  
cL = 0.3;
# cM = 5.0e-2 까지 높이면 fail to converge in optimization
E = 200.0; # loan demand 도 줄여봄 
cO = 100.0; # operation cost을 키워보기
dBar = 
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.ϵ, params_ex.σ
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
@show ex.moments
ex.simulation.paths.divSim

psi2(params::Params, d) = d >= 0 ? (d + params.dBar)^params.σ  - params.dBar^params.σ : 1 - exp(-d);

psi(params_ex, 100.0)
plot(-1000.0:0, x -> psi2(params, x))
ex.sol.eq.vFuncs.VF
#####################################################################################################################

@time eq = solve_model(params,false,params_cal.a,params_cal.b)
rl = eq.Rl_star
@time policy_woConst = Get_PolicyFuncs(params, eq.eqq.Iterobj_is, eq.eqq.vFuncs, eq.Rl_star)

policy_woConst.divPolicy # all dividences are negative
policy_woConst.taxPolicy # no tax because of too low Rl or too much loan supply (either because of low demand or negative d and extra loan supply)
policy_woConst.qBondPolicy
policy_woConst.failure
policy_woConst.lPolicy


E = 200.0
params_v1 = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_v1.E
@time vfi_v1 = solve_model(params_v1, false, params_cal.a, params_cal.b)
@show vfi_v1.Rl_star
@time policy_v1 = Get_PolicyFuncs(params_v1, vfi_v1.eqq.Iterobj_is, vfi_v1.eqq.vFuncs,vfi_v1.Rl_star)
policy_v1.divPolicy
policy_v1.nPrimePolicy
shocks_v1 = makeShockSequence(params_v1,vfi_v1.eq.vFuncs,bigT,bigN,bigJ);
paths_v1 = Initiate_Paths(shocks_v1, params);
Simulate_paths(paths_v1, policy_v1, vfi_v1.eq.vFuncs, params_v1, vfi_v1.Rl_star, trim, false)
mean(paths_v1.divSim)
vfi_v1.Rl_star
mean(paths_v1.failureSim) 
mean(paths_v1.nPrimeSim[paths_v1.failureSim .== 1])
mean(paths_v1.nPrimeSim[paths_v1.failureSim .== 0])
paths_v1.nPrimeSim[paths_v1.failureSim .== 1]
paths_v1.failureSim
paths_v1.capitalToDeposit
paths_v1.nPrimeSim
paths_v1.deltaSim
paths_v1.divSim



dBar = 2.0;
params_v2 = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_v2.dBar, params_v2.E
@time vfi_v2 = solve_model(params_v2, false, params_cal.a, params_cal.b)
@time policy_v2 = Get_PolicyFuncs(params_v2, vfi_v2.eqq.Iterobj_is, vfi_v2.eqq.vFuncs,vfi_v2.Rl_star)
policy_v2.divPolicy
shocks_v2 = makeShockSequence(params_v2, vfi_v2.eq.vFuncs, bigT, bigN, bigJ)
paths_v2 = Initiate_Paths(shocks_v2, params_v2)
Simulate_paths(paths_v2, policy_v2, vfi_v2.eq.vFuncs, params_v2, vfi_v2.Rl_star, trim, false)
mean(paths_v2.divSim)
mean(paths_v2.failureSim)
mean(paths_v2.nPrimeSim[paths_v2.failureSim .== 1])


dBar = 20.0;
params_v3 = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@time v3 = solve_simulate_and_moments(params_v3, params_cal, false)
mean(v3.simulation.paths.divSim)
mean(v3.simulation.paths.failureSim)
v3.Rl_star


dBar = 200.0; E = 150.0;
params_v4 = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δSSM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@time v4 = solve_simulate_and_moments(params_v4, params_cal, false);
mean(v4.simulation.paths.divSim)
mean(v4.simulation.paths.failureSim)
v4.Rl_star
v4.simulation.paths.divSim



 v= Initiate_vFunc(params)
shock = makeShockSequence(params,v,bigT, bigN, bigJ)
path = Initiate_Paths(shock, params)







ex = Initiate_PolicyFuncs(params)
ex.failure







policy_v2.divPolicy
policy_v2.lPolicy





policy_woConst.lPolicy 

policy_wConst.lPolicy

policy_woConst.failure
policy_wConst.failure
policy_woConst.nPrimePolicy
policy_wConst.nPrimePolicy

policy_woConst.NAV
policy_wConst.NAV

policy_woConst.sPolicy
policy_wConst.sPolicy

policy_woConst.bPolicy
policy_wConst.bPolicy

policy_woConst.lPolicy[3,:,:]
policy_woConst.sPolicy[3,:,:]
policy_woConst.bPolicy[:,:,:]
policy_woConst.failure[3,:,:, :]
policy_woConst.NAV[3,:,:,:]
policy_woConst.nPrimePolicy[3,:,:,:]

capital_const(l,s,b,delta) = (l + params.g*delta + s - delta - b)/(params.wr * l) 
params.α

capital_const.(policy_woConst.lPolicy[3,:,:],policy_woConst.sPolicy[3,:,:],policy_woConst.bPolicy[3,:,:],[params.deltaGrid[3]])
idel = 3; ilamb = 1; inInd =1

capital_const(policy_woConst.lPolicy[idel,ilamb,inInd],policy_woConst.sPolicy[idel,ilamb,inInd],policy_woConst.bPolicy[idel,ilamb,inInd],params.deltaGrid[idel])

(1-params.α*params.wr)*policy_woConst.lPolicy[idel,ilamb,inInd] 
policy_woConst.sPolicy[idel,ilamb,inInd] 
(1-params.g)*params.deltaGrid[idel] 
params.g
params.deltaGrid[idel] 

first(params.bGrid)
last(params.bGrid)

b_cap = min(last(params.bGrid),(1-params.α*params.wr)*policy_woConst.lPolicy[idel,ilamb,inInd] + policy_woConst.sPolicy[idel,ilamb,inInd]-(1-params.g)*params.deltaGrid[idel]) # b_start = (b_max+b_min)/2
l_cap = min(last(params.lGrid), (params.β - params.g) * params.deltaGrid[idel]) 
(1-params.α*params.wr)*policy_woConst.lPolicy[idel,ilamb,inInd] + policy_woConst.sPolicy[idel,ilamb,inInd]-(1-params.g)*params.deltaGrid[idel]

loanCon(l,delta) = (l + params.g*delta)/delta
loanCon.(policy_woConst.lPolicy[2,:,:],params.deltaGrid[2])

# after revising optimiztion code 
@time vfi_woConst_v2 = VFI(params, 1.1, false, 1000, 1e-3)
policy_wConst_v2 = Get_PolicyFuncs(params, vfi_woConst_v2.Iterobj_is, 1.1)

policy_wConst_v2.failure
policy_wConst_v2.nPrimePolicy
policy_wConst_v2.lPolicy[:,:,:]
policy_wConst_v2.lPolicy[2,:,:]
policy_wConst_v2.lPolicy[1,:,:]
policy_wConst_v2.sPolicy
policy_wConst_v2.bPolicy[:,:,:]
capital_const.(policy_wConst_v2.lPolicy[:,:,3],policy_wConst_v2.sPolicy[:,:,3],policy_wConst_v2.bPolicy[:,:,3],params.deltaGrid[:])
loanCon.(policy_wConst_v2.lPolicy[:,:,:],params.deltaGrid[:])
policy_wConst_v2.NAV


@time vfi_woConst_true = VFI(params, 1.1, true, 1000, 1e-3)
policy_woConst_true = Get_PolicyFuncs(params, vfi_woConst_true.Iterobj_is, 1.1)
policy_woConst_true.failure
policy_woConst_true.NAV
policy_woConst_true.lPolicy
policy_woConst_true.bPolicy
policy_woConst_true.sPolicy
vfi_woConst_true.vFuncs.qBond[:,:,:,2,3]


@time vfi_woConst_v3 = VFI(params, 1.1, false, 10, 1e-3)
vfi_woConst_v3.vFuncs.VF 
vfi_woConst_v3.vFuncsNew.VF
@time policy_woConst_v3 = Get_PolicyFuncs(params, vfi_woConst_v3.Iterobj_is, 1.1)
@time vfi_woConst_true_v3 = VFI(params, 1.1, true, 20, 1e-3)
@time policy_woConst_true_v3 = Get_PolicyFuncs(params, vfi_woConst_true_v3.Iterobj_is, 1.1)

vfi_woConst_true_v3.vFuncs.qBond[:,:,:,2,3]
policy_woConst_v3.lPolicy
policy_woConst_v3.sPolicy
policy_woConst_v3.bPolicy
vfi_woConst_v3.vFuncs.qBond
policy_woConst_true_v3.bPolicy
policy_woConst_v3.bPolicy
policy_woConst_true_v3.failure
policy_woConst_v3.failure

@time v_woConst_v4 = solve_model(params,false,0+eps(),1.1)
@time pol_woConst_v4 = Get_PolicyFuncs(params, v_woConst_v4.eq.Iterobj_is, v_woConst_v4.Rl_star)
shocks_v4 = makeShockSequence(params,v_woConst_v4.eq.vFuncs,bigT,bigN,bigJ);
paths_v4 = Initiate_Paths(shocks_v4, params);
Simulate_paths(paths_v4, pol_woConst_v4, v_woConst_v4.eq.vFuncs, params, v_woConst_v4.Rl_star, trim, false)


@time v_woConst_true_v4 = solve_model(params,true,0+eps(),1.1)
@time pol_woConst_true_v4 = Get_PolicyFuncs(params, v_woConst_true_v4.eq.Iterobj_is, v_woConst_true_v4.Rl_star)
shocks_true_v4 = makeShockSequence(params,v_woConst_true_v4.eq.vFuncs,bigT,bigN,bigJ);
paths_true_v4 = Initiate_Paths(shocks_true_v4, params);
Simulate_paths(paths_true_v4, pol_woConst_true_v4, v_woConst_true_v4.eq.vFuncs, params, v_woConst_true_v4.Rl_star, trim, true)

v_woConst_v4.eq.vFuncs.Γ
v_woConst_true_v4.eq.vFuncs.Γ

sum(paths_v4.failureSim)
sum(paths_true_v4.failureSim)

mean(paths_v4.qBondSim[trim+1:end, :, :])
mean(paths_true_v4.qBondSim[trim+1:end, :, :])

mean(paths_v4.bSim[trim+1:end, :, :])
mean(paths_true_v4.bSim[trim+1:end, :, :])

mean(paths_v4.lSim[trim+1:end, :, :])
mean(paths_true_v4.lSim[trim+1:end, :, :])

mean(paths_v4.leverageSim[trim+1:end, :, :])
mean(paths_true_v4.leverageSim[trim+1:end, :, :])

mean(paths_v4.assetPrimeSim[trim+1:end, :, :])
mean(paths_true_v4.assetPrimeSim[trim+1:end, :, :])

mean(paths_v4.assetSim[trim+1:end, :, :])
mean(paths_true_v4.assetSim[trim+1:end, :, :])

mean(paths_v4.divSim[paths_v4.failureSim.== 1][trim+1:end, :, :])
mean(paths_v4.divSim[paths_v4.failureSim.== 0][trim+1:end, :, :])

mean(paths_true_v4.divSim[paths_true_v4.failureSim .== 1][trim+1:end, :, :])
mean(paths_true_v4.divSim[paths_true_v4.failureSim .== 0][trim+1:end, :, :])

mean(paths_v4.nPrimeSim[paths_v4.failureSim.== 1][trim+1:end, :, :])
mean(paths_v4.nPrimeSim[paths_v4.failureSim .== 0][trim+1:end, :, :])
mean(paths_true_v4.nPrimeSim[paths_true_v4.failureSim .== 1][trim+1:end, :, :])
mean(paths_true_v4.nPrimeSim[paths_true_v4.failureSim .== 0][trim+1:end, :, :])


plot(paths_v4.nPrimeSim[trim+1:end-1, 5, 1])
plot!(paths_v4.failureSim[trim+1:end-1, 5, 1])
plot!(paths_v4.lSim[trim+1:end-1, 5, 1])
plot!(paths_v4.bSim[trim+1:end-1, 5, 1])
plot!(paths_v4.deltaSim[trim+1:end-1, 5, 1])
plot!(paths_v4.lambdaSim[trim+1:end-1, 5, 1])
plot!(paths_v4.qBondSim[trim+1:end-1, 5, 1])


plot(paths_true_v4.nPrimeSim[trim+1:end-1, 5, 1])
plot!(paths_true_v4.failureSim[trim+1:end-1, 5, 1])
plot!(paths_true_v4.lSim[trim+1:end-1, 5, 1])
plot(paths_true_v4.qBondSim[trim+1:end-1, 5, 1])

plot(paths_v4.deltaIndSim[trim+1:end-1, 5, 1])
plot(paths_v4.lambdaIndSim[trim+1:end-1, 5, 1])






vfi_diverge.vFuncs.VF - vfi_diverge.vFuncsNew.VF
# policy = Get_policyFuncs(params, eq.)

@time vfi = VFI(params, 0+eps(), false, 1000, 1e-3)
vfi.vFuncs.X
vfi.vFuncsNew.X
vfi.vFuncs.VF - vfi.vFuncsNew.VF
fieldnames(VFuncs)
@time sta = stationary_distribution(params,1.7,vfi.vFuncs,vfi.vFuncsNew,vfi.Iterobj_is,1000,1e-3)
sta.vFuncs.Γ
sta.vFuncsNew.Γ

@time sol = solve_model_given_r(1+eps(); Params = params, Regime = false)
solve_model_given_r_single = Rl -> solve_model_given_r(Rl; Params = params, Regime = false)
@time solve_model_given_r_single(params_cal.a)


sol.eq.vFuncs.VF - sol.eq.vFuncsNew.VF