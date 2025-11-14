pwd()
readdir()

include("parameters_v4_final.jl");
include("parameters_calib.jl");
include("main_v2_noLoanConstraint_capitalRequirement.jl");
include("main_v3_Reuse_dampning.jl");
Random.seed!(5877);

params_cal = Initiate_Params_cal(bigT,bigN,bigJ,trim,a,b,debt_to_liability,capital_to_deposit,loan_to_asset,loan_rate,cM,cO,cL,ϵ,E,dBar,σ);

lconstr = 2.4;
E = 220.0; # E = 185.0;
sigHat = 1.15;

params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@show params_ex.nGrid
@time eq = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
# sigHat = 1.2 일때 
@show eq.moments
sigHat = 1.15;
@show eq.moments
eq.vFuncs.VF

sigHat = 1.25;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq2 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq2.moments


# Nov 1st 

E = 200.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq2 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq2.moments


E = 200.0;
sigHat = 1.2;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq2 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq2.moments


E = 155.0;
sigHat = 1.15;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq3 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq3.moments

E = 190.0;
sigHat = 1.2;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq3 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq3.moments

E = 175.0;
sigHat = 1.2;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq3 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq3.moments

E = 185.0; # 더 높여야
sigHat = 1.2;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq3 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq3.moments

E = 189.0;
sigHat = 1.2;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq3 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq3.moments

E = 190.0;
sigHat = 1.2;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq3 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq3.moments

E = 189.5;
sigHat = 1.2;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq3 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq3.moments

E = 164.0;
sigHat = 1.15;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq3 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq3.moments

E = 165.0;
sigHat = 1.15;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq3 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq3.moments


E = 168.0;
sigHat = 1.15;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq3 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq3.moments

E = 174.0;
sigHat = 1.15;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq3 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq3.moments

E = 185.0;
sigHat = 1.15;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq3 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq3.moments


E = 179.0;
sigHat = 1.15;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq3 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq3.moments


E = 182.0;
sigHat = 1.15;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq3 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq3.moments


E = 184.0;
sigHat = 1.15;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq3 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq3.moments

lconstr = 2.4;
E = 183.0;
sigHat = 1.15;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq3 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq3.moments

E = 183.0;
sigHat = 1.15;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq3 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq3.moments


E = 182.0; ## Last try 
sigHat = 1.15;
l_stop = 100.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@show params_ex.lGrid
@time eq3 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq3.moments



E = 185.0;
sigHat = 1.15;
l_stop = 100.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@show params_ex.lGrid
@time eq3 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq3.moments


E = 210.0;
sigHat = 1.15;
l_stop = 100.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@show params_ex.lGrid
@show params_ex.cM, params_ex.cO, params_ex.cL
@time eq3 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq3.moments


###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################

vFuncs_false       = Initiate_vFunc(params_ex);
vFuncsNew_false    = Initiate_vFuncNew(params_ex);
Iterobj_is_false   = Initiate_MatrixIterObj_i(params_ex);

solve_model_given_r_single_false = Rl -> solve_model_given_r!(Rl;Params = params_ex,Regime = false,vFuncs = vFuncs_false,vFuncsNew = vFuncsNew_false,Iterobj_is = Iterobj_is_false)
Rl_safe = clamp(params_cal.a, 1.0 + eps(), 2.0 - eps())
@time Rl_star_false = find_zero(solve_model_given_r_single_false,Rl_safe,method = Roots.Secant();tol = 1e-2,maxevals = 100)

VFI!(params_ex, vFuncs_false, vFuncsNew_false, Iterobj_is_false, Rl_star_false, true, 1000, 1e-4)
stationary_distribution!(params_ex, Rl_star_false, vFuncs_false, vFuncsNew_false, Iterobj_is_false, 1000, 1e-4)

@time policy_false = Get_PolicyFuncs(params_ex, Iterobj_is_false, vFuncs_false, Rl_star_false, true);

@time simulation_false = simulate_and_moments(params_ex, params_cal, vFuncs_false, policy_false, Rl_star_false, true);

@show simulation_false.moments


### 최종 parameter: sigHat = 1.25; lconstr = 2.4; E = 250.0 with l_stop = 100.0
### 최종 parameter in Nov 1st: sigHat = 1.15; lconstr = 2.4; E = 200.0 with l_stop = 100.0

### Now go for counterfactual!
### solve line by line to debug
include("main_v3_Reuse_dampning.jl");

params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
params_cal = Initiate_Params_cal(bigT,bigN,bigJ,trim,a,b,debt_to_liability,capital_to_deposit,loan_to_asset,loan_rate,cM,cO,cL,ϵ,E,dBar,σ);


vFuncs_true       = Initiate_vFunc(params_ex);
vFuncsNew_true    = Initiate_vFuncNew(params_ex);
Iterobj_is_true   = Initiate_MatrixIterObj_i(params_ex);

solve_model_given_r_single_true = Rl -> solve_model_given_r!(Rl;Params = params_ex,Regime = true,vFuncs = vFuncs_true,vFuncsNew = vFuncsNew_true,Iterobj_is = Iterobj_is_true)

Rl_safe = clamp(params_cal.a, 1.0 + eps(), 2.0 - eps())
@time Rl_star = find_zero(solve_model_given_r_single_true,Rl_safe,method = Roots.Secant();tol = 1e-2,maxevals = 100)
@show Rl_star
# eq bond price : 1.039999

VFI!(params_ex, vFuncs_true, vFuncsNew_true, Iterobj_is_true, Rl_star, true, 1000, 1e-4)
stationary_distribution!(params_ex, Rl_star, vFuncs_true, vFuncsNew_true, Iterobj_is_true, 1000, 1e-4)

@time policy_true = Get_PolicyFuncs(params_ex, Iterobj_is_true, vFuncs_true, Rl_star, true);

@time simulation_true = simulate_and_moments(params_ex, params_cal, vFuncs_true, policy_true, Rl_star, true);

@show simulation_true.moments







policy_true.failure

vFuncs_true.qBond


vFuncs_true.qBond
vFuncsNew_true.qBond



######################################################################################################################################


















@time eq_counterfactual = solve_simulate_and_moments_reuse(params_ex, params_cal, true);

@show eq_counterfactual.moments

for il in eachindex(params_ex.lGrid)
    println(il)
end


params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq2_counterfactual = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq2_counterfactual.moments



@show eq.moments


eq.vFuncs.VF
@show eq.moments
@show eq.vFuncs.Γ
@show minimum(eq.paths.nSim[10:end-10,:,:])


lconstr = 2.4;
E = 220.0; # E = 185.0;
sigHat = 1.2;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@show params_ex.nGrid
@time eq2 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq2.moments





sum(eq.vFuncs.Γ[:,:,1])

eq.policy.failure
eq.policy.govSpend_bailout
eq.policy.govSpend_guarantee

params = params_ex




##################### counterfactual #######################

@time eq_true = solve_simulate_and_moments_reuse(params_ex, params_cal, true);











lconstr = 8.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat
@time eq = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
eq.moments

lconstr = 8.0;
E = 250.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
eq.moments

@show params_ex.E
lconstr = 2.0;
E = 160.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
eq.moments

@show params_ex.E
lconstr = 4.0;
E = 200.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
eq.moments

@show params_ex.E
lconstr = 3.0;
E = 180.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
eq.moments

@show params_ex.E
lconstr = 2.5;
E = 190.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
eq.moments

(params_ex.lconstr * 0.9801 - params_ex.g) * params_ex.deltaGrid[2]
eq.moments

@show params_ex.lconstr, params_ex.sigHat, params_ex.E
eq.moments

@show params_ex.lconstr, params_ex.sigHat, params_ex.E
eq.moments



########################################################################################################################################################



params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.sigHat

@time eq = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq.moments

@time eq = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show eq.moments
@show params_ex.sigHat

sigHat = 3.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.sigHat
@time eq = solve_simulate_and_moments_reuse(params_ex, params_cal, false);












@time eq_true = solve_simulate_and_moments_reuse(params_ex, params_cal, true);