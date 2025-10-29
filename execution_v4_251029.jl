pwd()
readdir()

include("parameters_v4_final.jl");
include("parameters_calib.jl");
include("main_v2_noLoanConstraint_capitalRequirement.jl");
include("main_v3_Reuse_dampning.jl");
Random.seed!(5877);

params_cal = Initiate_Params_cal(bigT,bigN,bigJ,trim,a,b,debt_to_liability,capital_to_deposit,loan_to_asset,loan_rate,cM,cO,cL,ϵ,E,dBar,σ);


lconstr = 8.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat
@time eq = solve_simulate_and_moments_reuse(params_ex, params_cal, false);





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