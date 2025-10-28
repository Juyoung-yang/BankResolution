pwd()
readdir()

# include("parameters.jl");
include("parameters_v3.jl");
include("parameters_calib.jl");
# include("main_v2_noLoanConstraint.jl");
include("main_v2_noLoanConstraint_capitalRequirement.jl");
# include("main_v2_LoanConstraint.jl");
include("main_v3_Reuse_dampning.jl");


Random.seed!(5877);
#####################################################################################################################
################# SET THE PARAMETERS ################################################################################ 
#####################################################################################################################

include("main_v2_noLoanConstraint_capitalRequirement.jl");

E = 145.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
params_cal = Initiate_Params_cal(bigT,bigN,bigJ,trim,a,b,debt_to_liability,capital_to_deposit,loan_to_asset,loan_rate,cM,cO,cL,ϵ,E,dBar,σ);
@show params_ex.sigHat, params_ex.lconstr
@show params_ex.F[3,:]
@show params_ex.cL, params_ex.E
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@show params_ex.wr, params_ex.α

@time ex = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show ex.moments

E = 135.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@time ex = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show ex.moments








params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex1 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show ex1.moments

lconstr = 2.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.sigHat, params_ex.lconstr, params_ex.F[3,:]
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@show params_ex.λL, params_ex.λH
@time ex1 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show ex1.moments

sigHat = 1.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.sigHat, params_ex.lconstr, params_ex.F[3,:]
@time ex1 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show ex1.moments

lconstr = 1.2;
cL = -0.1;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr, params_ex.cL, params_ex.sigHat
@time ex1 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show ex1.moments

λL = 0.00275;
λM = 0.00530;
λH = 0.00745;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.λL, params_ex.λH
@time ex1 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show ex1.moments

@show params_ex.E
E = 140.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.E
@time ex1 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show ex1.moments

F[3,:] = [0.0 0.88 0.12]
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.F[3,:] 
@time ex1 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show ex1.moments

lconstr = 2.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr
@time ex1 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show ex1.moments

lconstr = 1.2;
F = [0.8065 0.1685 0.025; 0.1595 0.778 0.0625; 0.0 0.8812 0.1188]; # April
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.F, params_ex.lconstr, params_ex.λH
@time ex1 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show ex1.moments

@show params_ex.E
E = 150.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.E
@time ex1 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show ex1.moments
ex1.paths.failureSim[40:50, 1,1]
ex1.paths.loanToDeposit[30:40, 1,1]
ex1.paths.lSim[30:40, 1,1]
ex1.paths.deltaSim[50:60, 1,1]
ex1.paths.lambdaSim[30:40, 1,1]

# 은행 실패를 줄이려면 high deposit인 은행들이 과도하게 loan을 늘리지 말아야; loan supply cost를 늘려보자 
cM = 1.3e-3;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM
@time ex1 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show ex1.moments


cL = -0.003;
ρ = 0.1;
cM = 1.3e-4;
cO = 0.1;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cL, params_ex.ρ, params_ex.cM, params_ex.cO
@show params_ex.F, params_ex.sigHat, params_ex.lconstr
@time ex1 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show ex1.moments
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 3) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 1) .& (ex1.paths.failureSim[10:end-10,:,:] .== true))
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 3) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 1) .& (ex1.paths.failureSim[10:end-10,:,:] .== false))
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 2) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 1) .& (ex1.paths.failureSim[10:end-10,:,:] .== true))
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 2) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 3) .& (ex1.paths.failureSim[10:end-10,:,:] .== true))

# constraint의 alpha = 0.13, 실패율의 alpha = 0.8;
α = 0.08;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.α
@show params_ex.cL, params_ex.ρ, params_ex.cM, params_ex.cO
@show params_ex.F, params_ex.sigHat, params_ex.lconstr

@time ex1 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
ex1.moments

sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 3) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 1) .& (ex1.paths.failureSim[10:end-10,:,:] .== true)), sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 3) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 1) .& (ex1.paths.failureSim[10:end-10,:,:] .== false))
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 3) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 2) .& (ex1.paths.failureSim[10:end-10,:,:] .== true)), sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 3) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 2) .& (ex1.paths.failureSim[10:end-10,:,:] .== false))
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 3) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 3) .& (ex1.paths.failureSim[10:end-10,:,:] .== true)), sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 3) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 3) .& (ex1.paths.failureSim[10:end-10,:,:] .== false))

sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 2) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 1) .& (ex1.paths.failureSim[10:end-10,:,:] .== true)), sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 2) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 1) .& (ex1.paths.failureSim[10:end-10,:,:] .== false))
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 2) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 2) .& (ex1.paths.failureSim[10:end-10,:,:] .== true)), sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 2) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 2) .& (ex1.paths.failureSim[10:end-10,:,:] .== false))
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 2) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 3) .& (ex1.paths.failureSim[10:end-10,:,:] .== true)), sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 2) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 3) .& (ex1.paths.failureSim[10:end-10,:,:] .== false))

sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 1) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 1) .& (ex1.paths.failureSim[10:end-10,:,:] .== true)), sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 1) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 1) .& (ex1.paths.failureSim[10:end-10,:,:] .== false))
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 1) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 2) .& (ex1.paths.failureSim[10:end-10,:,:] .== true)), sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 1) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 2) .& (ex1.paths.failureSim[10:end-10,:,:] .== false))
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 1) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 3) .& (ex1.paths.failureSim[10:end-10,:,:] .== true)), sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 1) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 3) .& (ex1.paths.failureSim[10:end-10,:,:] .== false))

mean( ex1.paths.loanToDeposit[10:end-10, :, :][ex1.paths.deltaIndSim[10:end-10, :, :] .== 1 ])
mean( ex1.paths.loanToDeposit[10:end-10, :, :][ex1.paths.deltaIndSim[10:end-10, :, :] .== 2 ])
mean( ex1.paths.loanToDeposit[10:end-10, :, :][ex1.paths.deltaIndSim[10:end-10, :, :] .== 3 ])

mean( ex1.paths.loanToDeposit[10:end-10, :, :][ (ex1.paths.deltaIndSim[10:end-10, :, :] .== 3) .& (ex1.paths.lambdaIndSim[10:end-10, :, :] .== 3) ])
mean( ex1.paths.loanToDeposit[10:end-10, :, :][ (ex1.paths.deltaIndSim[10:end-10, :, :] .== 3) .& (ex1.paths.lambdaIndSim[10:end-10, :, :] .== 2) ])
mean( ex1.paths.loanToDeposit[10:end-10, :, :][ (ex1.paths.deltaIndSim[10:end-10, :, :] .== 3) .& (ex1.paths.lambdaIndSim[10:end-10, :, :] .== 1) ])

mean( ex1.paths.loanToDeposit[10:end-10, :, :][ (ex1.paths.deltaIndSim[10:end-10, :, :] .== 2) .& (ex1.paths.lambdaIndSim[10:end-10, :, :] .== 3) ])
mean( ex1.paths.loanToDeposit[10:end-10, :, :][ (ex1.paths.deltaIndSim[10:end-10, :, :] .== 2) .& (ex1.paths.lambdaIndSim[10:end-10, :, :] .== 2) ])
mean( ex1.paths.loanToDeposit[10:end-10, :, :][ (ex1.paths.deltaIndSim[10:end-10, :, :] .== 2) .& (ex1.paths.lambdaIndSim[10:end-10, :, :] .== 1) ])


cM = 1.3e-5;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM
@time ex11 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
ex11.moments
mean( ex11.paths.loanToDeposit[10:end-10, :, :][ex11.paths.deltaIndSim[10:end-10, :, :] .== 1 ])
mean( ex11.paths.loanToDeposit[10:end-10, :, :][ex11.paths.deltaIndSim[10:end-10, :, :] .== 2 ])
mean( ex11.paths.loanToDeposit[10:end-10, :, :][ex11.paths.deltaIndSim[10:end-10, :, :] .== 3 ])



@show params_ex.cO
cO = 0.2;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cO
@time ex12 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
ex12.moments


@show params_ex.cL
ρ = 0.01;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cL, params_ex.ρ
@time ex12 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
ex12.moments

ρ = 0.4;
wr = 0.428
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cL, params_ex.ρ
@show params_ex.wr
@time ex12 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show ex12.moments

mean( ex12.paths.loanToDeposit[10:end-10, :, :][ex12.paths.deltaIndSim[10:end-10, :, :] .== 1 ])
mean( ex12.paths.loanToDeposit[10:end-10, :, :][ex12.paths.deltaIndSim[10:end-10, :, :] .== 2 ])
mean( ex12.paths.loanToDeposit[10:end-10, :, :][ex12.paths.deltaIndSim[10:end-10, :, :] .== 3 ])


wr = 0.0
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.wr
@time ex2 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show ex2.moments

mean( ex2.paths.loanToDeposit[10:end-10, :, :][ex2.paths.deltaIndSim[10:end-10, :, :] .== 1 ])
mean( ex2.paths.loanToDeposit[10:end-10, :, :][ex2.paths.deltaIndSim[10:end-10, :, :] .== 2 ])
mean( ex2.paths.loanToDeposit[10:end-10, :, :][ex2.paths.deltaIndSim[10:end-10, :, :] .== 3 ])


lconstr = 2.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.lconstr
@time ex2 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show ex2.moments
# 0.4222072378290964, 0.6063174121654216, 0.23862158348918933, 0.32028409090909093, 1.077808976732123], other_moments = [-0.009073734145391654, 0.06409875552695055, -0.14830324642469547, 3.522686421651657, 19.869989351279408, 1.1830184258017915])

ex12.moments
















params_ex.cL
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 3) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 1) .& (ex1.paths.failureSim[10:end-10,:,:] .== true))
sum( (p_ex.deltaIndSim[10:end-10, :, :] .== 3) .& (p_ex.lambdaIndSim[11:end-9, :, :] .== 2) .& (p_ex.failureSim[10:end-10,:,:].== true))
sum( (p_ex.deltaIndSim[10:end-10, :, :] .== 3) .& (p_ex.lambdaIndSim[11:end-9, :, :] .== 3) .& (p_ex.failureSim[10:end-10,:,:].== true))

ρ = 0.1;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM




include("main_v2_noLoanConstraint_capitalRequirement.jl");

shocks_ex = makeShockSequence(params_ex,ex1.vFuncs,params_cal.bigT,params_cal.bigN,params_cal.bigJ);
p_ex = Initiate_Paths(shocks_ex, params_ex);
series_ex = Simulate_paths(p_ex, ex1.policy, ex1.vFuncs,params_ex,ex1.Rl_star, params_cal.trim, false)
p_ex.failureSim[61:71, 1,1]
p_ex.deltaSim[61:71, 1,1]
p_ex.lambdaSim[61:71, 1,1]
p_ex.nSim[45:55, 1,1]
p_ex.lSim[45:55, 1,1]
p_ex.sSim[45:55, 1,1]
p_ex.bSim[45:55, 1,1]
p_ex.deltaSim[45:55, 1,1]
p_ex.lambdaSim[45:55, 1,1]



sum( (p_ex.deltaIndSim[10:end-10, :, :] .== 1) .& (p_ex.lambdaIndSim[11:end-9, :, :] .== 1) .& (p_ex.failureSim[10:end-10,:,:] .== true))
sum( (p_ex.deltaIndSim[10:end-10, :, :] .== 1) .& (p_ex.lambdaIndSim[11:end-9, :, :] .== 2) .& (p_ex.failureSim[10:end-10,:,:].== true))
sum( (p_ex.deltaIndSim[10:end-10, :, :] .== 1) .& (p_ex.lambdaIndSim[11:end-9, :, :] .== 3) .& (p_ex.failureSim[10:end-10,:,:].== true))
sum( (p_ex.deltaIndSim[10:end-10, :, :] .== 2) .& (p_ex.lambdaIndSim[11:end-9, :, :] .== 1) .& (p_ex.failureSim[10:end-10,:,:] .== true))
sum( (p_ex.deltaIndSim[10:end-10, :, :] .== 2) .& (p_ex.lambdaIndSim[11:end-9, :, :] .== 2) .& (p_ex.failureSim[10:end-10,:,:].== true))
sum( (p_ex.deltaIndSim[10:end-10, :, :] .== 2) .& (p_ex.lambdaIndSim[11:end-9, :, :] .== 3) .& (p_ex.failureSim[10:end-10,:,:].== true))
sum( (p_ex.deltaIndSim[10:end-10, :, :] .== 3) .& (p_ex.lambdaIndSim[11:end-9, :, :] .== 1) .& (p_ex.failureSim[10:end-10,:,:] .== true))
sum( (p_ex.deltaIndSim[10:end-10, :, :] .== 3) .& (p_ex.lambdaIndSim[11:end-9, :, :] .== 2) .& (p_ex.failureSim[10:end-10,:,:].== true))
sum( (p_ex.deltaIndSim[10:end-10, :, :] .== 3) .& (p_ex.lambdaIndSim[11:end-9, :, :] .== 3) .& (p_ex.failureSim[10:end-10,:,:].== true))
sum( (p_ex.deltaIndSim[10:end-10, :, :] .== 3) .& (p_ex.lambdaIndSim[11:end-9, :, :] .== 1) .& (p_ex.failureSim[10:end-10,:,:] .== false))
sum( (p_ex.deltaIndSim[10:end-10, :, :] .== 3) .& (p_ex.lambdaIndSim[11:end-9, :, :] .== 2) .& (p_ex.failureSim[10:end-10,:,:].== false))
sum( (p_ex.deltaIndSim[10:end-10, :, :] .== 3) .& (p_ex.lambdaIndSim[11:end-9, :, :] .== 3) .& (p_ex.failureSim[10:end-10,:,:].== false))
sum( (p_ex.deltaIndSim[10:end-10, :, :] .== 2) .& (p_ex.lambdaIndSim[11:end-9, :, :] .== 3) .& (p_ex.failureSim[10:end-10,:,:].== false))


params_ex.nGrid[p_ex.nInitialIndSim]

ex1.paths.failureSim[20:30, 1, 1]
sum(ex1.paths.failureSim[10:end-10, :,:] .> -1)
sum(ex1.paths.failureSim[10:end-10, :,:] .== 1)
ex1.paths.failureSim[35:45, 1,1]
ex1.paths.lSim[25:35, 1,1]
ex1.paths.deltaSim[35:45, 1,1]
ex1.paths.lSim[35:45, 1,1]
ex1.paths.bSim[35:45, 1,1]
ex1.paths.sSim[35:45, 1,1]
ex1.paths.leverageSim[35:45, 1,1]
ex1.paths.nPrimeSim[45, 1,1]
n_success(ex1.paths.lSim[45,1,1], ex1.paths.sSim[45,1,1],ex1.paths.bSim[45,1,1],ex1.paths.lambdaSim[46,1,1],ex1.paths.deltaSim[45,1,1],params_ex,ex1.Rl_star)
ex1.paths.nSim[46,1,1]
params_ex.nGrid


NAV(l,s,b,lambda,delta,params,Rl) = Rl*( (1-lambda)*l + params.g*delta) + (1+Rl)*s-delta-b 
tax(l,s,b,lambda,delta,params,Rl) = params.τC*max(0.0, (Rl-1)*((1-lambda)*l +params.g*delta) + params.Rf*s - params.Rf *(delta + b)) 
n_success(l,s,b,lambda,delta,params,Rl) = NAV(l,s,b,lambda,delta,params,Rl) -tax(l,s,b,lambda,delta,params,Rl)







ex1.paths.deltaSim[25:35, 1,1] .* ex1.paths.failureSim[25:35, 1,1]
ex1.paths.lambdaSim[26:36, 1,1] .* ex1.paths.failureSim[25:35, 1,1]
ex1.paths.lambdaSim[25:35, 1,1]
count(==(0.1796), ex1.paths.lambdaSim[11:end-9,:,:] .* ex1.paths.failureSim[10:end-10,:,:])
count(==(0.0047), ex1.paths.lambdaSim[11:end-9,:,:] .* ex1.paths.failureSim[10:end-10,:,:])
count(==(0.0074), ex1.paths.lambdaSim[11:end-9,:,:] .* ex1.paths.failureSim[10:end-10,:,:])
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 1) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 1) .& (ex1.paths.failureSim[10:end-10,:,:] .== true))
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 1) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 2) .& (ex1.paths.failureSim[10:end-10,:,:].== true))
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 1) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 3) .& (ex1.paths.failureSim[10:end-10,:,:].== true))
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 2) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 1) .& (ex1.paths.failureSim[10:end-10,:,:] .== true))
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 2) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 2) .& (ex1.paths.failureSim[10:end-10,:,:].== true))
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 2) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 3) .& (ex1.paths.failureSim[10:end-10,:,:].== true))
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 3) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 1) .& (ex1.paths.failureSim[10:end-10,:,:] .== true))
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 3) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 2) .& (ex1.paths.failureSim[10:end-10,:,:].== true))
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 3) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 3) .& (ex1.paths.failureSim[10:end-10,:,:].== true))
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 3) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 1) .& (ex1.paths.failureSim[10:end-10,:,:] .== false))
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 3) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 2) .& (ex1.paths.failureSim[10:end-10,:,:].== false))
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 3) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 3) .& (ex1.paths.failureSim[10:end-10,:,:].== false))
sum( (ex1.paths.deltaIndSim[10:end-10, :, :] .== 2) .& (ex1.paths.lambdaIndSim[11:end-9, :, :] .== 3) .& (ex1.paths.failureSim[10:end-10,:,:].== false))


mask = ex1.paths.failureSim[10:end-10,:,:] .== 1
noFail = ex1.paths.failureSim[10:end-10,:,:] .== 0
mean( ex1.paths.loanToAsset[10:end-10, :, :][mask]), mean( ex1.paths.loanToAsset[10:end-10, :, :][noFail])
mean( ex1.paths.nPrimeSim[10:end-10, :, :][mask]), mean( ex1.paths.nPrimeSim[10:end-10, :, :][noFail])
mean(ex1.paths.capitalToDeposit[10:end-10, :, :][mask]), mean(ex1.paths.capitalToDeposit[10:end-10, :, :][noFail])


dist = kron(params_ex.H, params_ex.F)
dist ./= sum(dist)














params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
params_ex.sigHat
psi.([params_ex], -1.0:1.0:10.0)
sigHat = 2.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
params_ex.sigHat
psi.([params_ex], -1.0:1.0:10.0)


#####################################################################################################################
################# main code: calibration ############################################################################ 
#####################################################################################################################


using Sobol, Distributions, Distributed
@time calibration(params, params_cal, false)

#####################################################################################################################
#################  change parameters one by one #####################################################################




params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
# 예대율 있을 때
ex.moments
# 예대율 없을 때
ex.moments
# 예대율 규제 조정하고 신주발행 비용 증가: 대출은 증가하였지만 여전히 실패가 너무 많이 발생; 비용 증가는 배당금 지급 축소에 크게 영향을 미치지 않음
ex.moments
# 이자율을 맞추기 위해 대출수요를 줄이자
E = 150.0
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
# 여전히 대출 금리가 높아서 수요를 더 줄여야 
E = 130.0
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
ex.paths.nPrimeSim[params_cal.trim:20, 1, 1]
ex.paths.deltaSim[params_cal.trim:20, 1, 1]
ex.paths.lSim[params_cal.trim:20, 1, 1]
ex.paths.bSim[params_cal.trim:20, 1, 1]
ex.paths.divSim[params_cal.trim:20, 1, 1]
# 예대율 규제를 더욱 완화해보기 -> 은행 실패가 25%로 일어나고 (즉, 어느정도는 대출 규모가 있어야) 
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
E = 150.0
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
ex.paths.nPrimeSim[params_cal.trim:20, 1, 1] ./ ex.paths.deltaSim[params_cal.trim:20, 1, 1]
ex.paths.failureSim[params_cal.trim:20, 1, 1]
mean(ex.paths.nPrimeSim[params_cal.trim:80, :, :] ./ ex.paths.deltaSim[params_cal.trim:80, :, :])
# 아직 자기자본이 살짝 높음, 25% 은행 실패는 모델의 lambda shock 이 엄청 크기 때문 
# 증자와 자기자본을 좀 줄여보자, 증자 cost를 높여보자 
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
# 유동성 cost를 negative로 해서 debtToLiability를 낮춰보자 
cL = -0.2
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
# deposit으로 fund rasing 하는게 오히려 이득이니 dividend가 positive로 나오기 시작! Good
mean(ex.paths.nPrimeSim[params_cal.trim:80, :, :] ./ ex.paths.deltaSim[params_cal.trim:80, :, :])
# 자기자본도 꽤나 낮음 great!
# 이제 debt to liability을 낮추고 (혹은 debt을 낮추고), loan to asset을 높이고 (save asset 투자를 줄이고) 
cM = 1.3e-7
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
mean(ex.paths.sSim[params_cal.trim:80, :, :] ./ ex.paths.deltaSim[params_cal.trim:80, :, :])
mean(ex.paths.bSim[params_cal.trim:80, :, :] ./ ex.paths.deltaSim[params_cal.trim:80, :, :])
mean(ex.paths.lSim[params_cal.trim:80, :, :] ./ ex.paths.deltaSim[params_cal.trim:80, :, :])
# 예금이자율을 80% 수준으로 낮추고 (혹은 유동성 비용을 더 낮추고), 대출을 늘리기 위해 loan cost을 낮추자
##### 다시 배당금이 양수인 parameter를 찾자 
cL = -0.25;
cM = 1.3e-5;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
# 다시 배당금이 양수인 parameter를 찾자 
dBar = 1.0e-100;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
# dBar를 조정해도 그닥 효과는 없음..
# 배당금은 양수이고 크기고 작지만, loan이 작고 debt은 큼
dBar = 1.0;
cL = -0.25;
cO = 0.01;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
# what if very high operation cost?
cO = 10.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
cO = 100.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
# 오히려 operation cost를 낮춰보자 -> div이 양수가 될 것 
cO = 0.001;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
# 신주발행 곡률 높이기: 2.0 에서 2.5로
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
# 현실적인 이자율을 위해 demand 낮추기
E = 110.0
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments 
E = 145.0
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
# parameters: cL, cO, cM, E, div params (dBar, sigma, square)
# target moments: 4개 <- E, square + alpha 2 개
@time ex_true = solve_simulate_and_moments(params_ex, params_cal, true);

# 예대율 규제 parameter가 1.2 일때
params.β - params.g
cM = 1.3e-5;
cO = 0.1;
cL = 0.1;
E = 145.0;
dBar = 1.0e-100;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments

# 예대율 규제 parameter가 1.5일때
cM = 1.3e-5;
cO = 0.1;
cL = -0.25;
E = 145.0;
dBar = 1.0e-100;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments

# 예대율 규제가 1.3일떄
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments

# 예대율 규제가 1.2일떄
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
# 배당률은 괜춘스 -> 예대율 규제는 1.2로 고정 
# cM을 April 와 같이 해도 되는지; 배당에 대한 인센티브가 줄어들긴 함 
dBar = 1.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
# okay, keep cM, σ as April 
# maybe I can adjust coveture: 1 -> 1.5
include("main_v2_noLoanConstraint.jl");
dBar = 1.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
# adjust loan demand a little bit 
E = 140.0
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
# demand 가 좀 작다, 늘리자
E = 143.0
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
## great! finalize parameter to 
# cO = 0.1
# cL = -0.25
# E = 143.0
# Phi-power = 1.5

F[3,:] = [0.0 0.88 0.12]
cL = -0.25
cO = 0.1
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
params_ex.F
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
# failure rule: unsatiafication of capital requirement
ex.moments

# # failure rule: zero NAV
ex.moments
# even with markov process, no change in failure rate; rather lower loan and higher borrowing 

include("main_v3_Reuse_dampning.jl")
# 시간을 단축하기 위해 reuse + damping을 넣은 함수 별도제작 
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@show params_ex.F
@time ex = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show ex.moments

F = [0.625 0.2875 0.0875; 0.22449 0.500 0.27551; 0.075949 0.341772 0.582279] # from 고정이하여신비율 
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@show params_ex.F
@time ex = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show ex.moments

E = 150.0
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@show params_ex.F
@time ex = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show ex.moments

# 자기자본으로 바꿨을 때 균형 결과 
cL = -0.25;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@show params_ex.F
@time ex1 = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show ex1.moments
ex1.vFuncs.Γ

cL = -0.15
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments_reuse(params_ex, params_cal, false);
@show ex.moments

















A = kron(params_ex.F, params_ex.H)
using LinearAlgebra
ρ = maximum(abs.(eigvals(A)))
println("Spectral radius ρ(A) = ", ρ)

F = [0.625 0.2875 0.0875; 0.22449 0.500 0.27551; 0.075949 0.341772 0.582279] # from 고정이하여신비율 
AA = kron(F, params_ex.H)
ρρ = maximum(abs.(eigvals(AA)))
println("Spectral radius ρ(A) = ", ρρ)





###########################################################################################################################
###########################################################################################################################
###########################################################################################################################


# now turn to counterfactural
@time ex_counter = solve_simulate_and_moments(params_ex, params_cal, true);





















dBar = 100000000.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments



dBar = 1.0e-10;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments


dBar = 1.0;
E = 200.0;
cM = 1.0e-2;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
# cM 을 너무 높여버리니 capital 은 줄어들었는데 loan은 증가하고 bank failure를 안함 

dBar = 1.0;
E = 200.0;
σ = 0.5
cM = 1.0e-2;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar, params_ex.σ
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
ex.moments



σ = 0.01;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar, params_ex.σ
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
# 배당금 곡률을 바꾼다 한들 크게 moment가 변하지 않음

dBar = 500.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar, params_ex.σ
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments


cM = 1.0e-5;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar, params_ex.σ
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments

cM = 5.0e-4;
σ = 0.9932 # (KEEP IT FIXED) 
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar, params_ex.σ
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments

cM = 1.0e-5;
σ = 0.9932 # (KEEP IT FIXED) 
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar, params_ex.σ
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments





cL = 0.3;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar, params_ex.σ
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
# bamk fail 이 계속 안나오네 


cM = 1.3e-5;
cO = 0.1;
cL = 0.1;
E = 200.0;
dBar = 1000.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar, params_ex.σ
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments
# 바로  bank failure 가 일어남 

cM = 1.3e-5;
cO = 0.1;
cL = 0.1;
E = 200.0;
dBar = 10.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar, params_ex.σ
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments

cM = 1.3e-5;
cO = 0.1;
cL = 0.1;
E = 200.0;
dBar = 1.0;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar, params_ex.σ
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments

cM = 1.3e-5;
cO = 0.1;
cL = 0.1;
E = 200.0;
dBar = 1.0e-3;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar, params_ex.σ
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments


# liquidity cost 를 높여보자 
cM = 1.3e-5;
cO = 0.1;
cL = 0.4;
E = 200.0;
dBar = 1.0e-3;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar, params_ex.σ
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments


# operation cost를 높여보자 
cM = 1.3e-5;
cO = 10.0;
cL = 0.1;
E = 200.0;
dBar = 1.0e-3;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar, params_ex.σ
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments

# moderate operation cost 
cM = 1.3e-5;
cO = 1.0;
cL = 0.1;
E = 200.0;
dBar = 1.0e-3;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar, params_ex.σ
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments

endT =  Int(ex.paths.bigT - 3)
mean(ex.paths.sSim[trim:endT, :, :])

mean(ex.paths.lSim[trim:endT, :, :])

mean(ex.paths.bSim[trim:endT, :, :])

mean(ex.paths.deltaSim[trim:endT, :, :])

mean(ex.paths.divSim[trim:endT, :, :])
mean(ex.paths.capitalToDeposit[trim:endT, :, :])
mean(ex.paths.failureSim[trim:endT, :, :])

mean(ex.paths.sSim[trim:endT, :, :])/mean(ex.paths.assetSim[trim:endT, :, :])
mean(ex.paths.bSim[trim:endT, :, :])/mean(ex.paths.liabilitySim[trim:endT, :, :])
mean(ex.paths.deltaSim[trim:endT, :, :])/mean(ex.paths.lSim[trim:endT, :, :])
mean(ex.paths.lSim[trim:endT, :, :] ./ ex.paths.deltaSim[trim:endT, :, :])

ex.paths.lSim[trim:endT, 1, 1]
ex.paths.deltaSim[trim:endT, 1, 1]

params_ex.dBar
psi(params::Params, d) = d >= 0 ? (d + 0.001)^params.σ  - 0.001^params.σ : 1 - exp(-d);
psi2(params::Params, d) = d >= 0 ? (d + 0.0000000000001)^params.σ  - 0.0000000000001^params.σ : 1 - exp(-d)^5;
psi_ex(d) = d >= 0 ? (d+1.0)^0.9 - 1.0^0.9 : 1-exp(-d)^1.0
xx = 1.0e-100
yy = 0.9
zz = 1.1
psi_ex2(d) = d >= 0 ? (d+xx)^yy - xx^yy : 1-exp(-d)^zz

plot(-10:1:10, a -> psi_ex(a))
plot!(-10:1:10, a -> psi_ex(a))

@show psi_ex.(-10:1:10)'
@show psi_ex2.(-10:1:10)'
@show psi_ex2.(-10:1:10)'
@show psi_ex2.(-10:1:10)'




plot(-7:1:10, y -> psi(params_ex,y))

plot!(-7:1:10, y -> psi2(params_ex,y))

cM = 1.3e-5;
cO = 1.0;
cL = 0.1;
E = 200.0;
dBar = 1.0e-3;
params_ex = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
@show params_ex.cM, params_ex.cO, params_ex.cL, params_ex.E, params_ex.dBar, params_ex.σ
@time ex = solve_simulate_and_moments(params_ex, params_cal, false);
ex.moments


params_ex.β - params_ex.g











































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