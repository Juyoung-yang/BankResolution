pwd()

include("parameters.jl");
include("parameters_calib.jl");
include("main_v2_noLoanConstraint.jl");

params = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
params_cal = Initiate_Params_cal(bigT,bigN,bigJ,trim,a,b,debt_to_liability,capital_to_deposit,loan_to_asset,cM,cO,cL,ϵ,E,dBar);

Random.seed!(5877);

# initial parameter guess 
@time Init = solve_simulate_and_moments(params,params_cal,false)

@show params.E
E = 150.0
params = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params.E
@time Init = solve_simulate_and_moments(params,params_cal,false)