pwd()
readdir()

include("parameters_v4_final.jl");
include("parameters_calib.jl");
include("main_v2_noLoanConstraint_capitalRequirement.jl");
include("main_v3_Reuse_dampning.jl");
Random.seed!(5877);

params_cal = Initiate_Params_cal(bigT,bigN,bigJ,trim,a,b,debt_to_liability,capital_to_deposit,loan_to_asset,loan_rate,cM,cO,cL,ϵ,E,dBar,σ);

lconstr = 2.4;
E = 188.0; # E = 185.0;
sigHat = 2.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,ρ_bailout,g,ξ,cF,dBar,σ,τC,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.sigHat, params_ex.E
@time eq = solve_simulate_and_moments_reuse(params_ex, params_cal, false);

@show params_ex.lambdaGrid


params = params_ex














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