pwd()
readdir()

include("C:\\Users\\master\\Desktop\\(2025) banking regulation\\BankResolution\\parameters.jl");
include("C:\\Users\\master\\Desktop\\(2025) banking regulation\\BankResolution\\main_v2.jl");

include("parameters.jl");
include("main_v2_.jl");
include( "main_v2_250922.jl")


# 0. set the parameters
Random.seed!(5877);
params = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,α1,α2,α3,δL,δM,δH,cM,cO,cL,ϵ,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)

# 1-1. solve the model under given params: counterfactual
@time sol = solve_model(params,true, 1.0+eps(), 2-eps()); # find Rl_star from 1.0 to 1.9 (since it is gross return)
@show sol.Rl_star
@show sol.eqq.vFuncsNew.Γ[:, :,1]
@time policy_ex = Get_PolicyFuncs(params, sol.eq.Iterobj_is, sol.Rl_star); # get policy functions from the solution

# 1-2. solve the model under given params: benchmark
@time sol_false = solve_model(params,false, 1.0+eps(), 2-eps());
@time policy_false = Get_PolicyFuncs(params, sol_false.eqq.Iterobj_is, sol_false.Rl_star); # get policy functions from the solution
@show sol_false.Rl_star

# 1-1. solve the model under given params: counterfactual
@show sol_false.eq.vFuncs.VF == sol_false.eqq.vFuncs.VF
@show sum(sol_false.eq.vFuncs.X)
@show typeof(sol_false.eq.Iterobj_is[1,1]), fieldnames(typeof(sol_false.eq.Iterobj_is[1,1]))
@show typeof(sol_false.eq.Iterobj_is[1,1]), size(sol_false.eq.Iterobj_is[1,1].solution)
@show sol_false.eq.Iterobj_is[1,1].solution[3] .> sol_false.eq.Iterobj_is[1,3].solution[3]
@show sol_false.eq.Iterobj_is[1,1].solution[2] .> sol_false.eq.Iterobj_is[3,1].solution[2]

@show fieldnames(typeof(policy_false)), size(policy_false.lPolicy)

mean(policy_false.lPolicy[:, 1,:]), mean(policy_false.lPolicy[:, 2, :]), mean(policy_false.lPolicy[:, 3,:])
mean(policy_false.failure)
mean(policy_false.sPolicy[:, 1,:]), mean(policy_false.sPolicy[:, 2, :]), mean(policy_false.sPolicy[:, 3,:])
mean(policy_false.bPolicy[:, 1,:]), mean(policy_false.bPolicy[:, 2, :]), mean(policy_false.bPolicy[:, 3,:])

policy_false.bPolicy[:, :,2]
NAV(l::T,s::T,b::T,lambda::T,δ::T) where {T<:Real,S<:Integer,F<:Bool} = 1.2*( (1-lambda)*l + params.g*δ) + (1+params.Rf)*s - δ - b; # 순자산, net asset value 
idel = 1; ilam = 3; inn = 1;
NAV.(policy_false.lPolicy[idel, ilam, inn], policy_false.sPolicy[idel, ilam, inn], policy_false.bPolicy[idel, ilam, inn], params.lambdaGrid, params.deltaGrid[idel])
1.2^params.ϵ, sol_false.Rl_star^params.ϵ, 1.0^params.ϵ

aggre_loan_supply2(params, sol_false.eq.vFuncs, sol_false.eq.Iterobj_is)

function plot_the_graph()
    # optimal l 
    A = [sol_false.eq.Iterobj_is[i,1].solution[2][1] for i in 1:3] # L
    B = [sol_false.eq.Iterobj_is[i,1].solution[5][1] for i in 1:3] # S

    # optimal b
    pl_lowDelta = plot(A, label = "sin(x)")
    plot!(B, label = "sin(x)", linestyle = :dash)

    pl_highDelta = plot()
    plot!() 


    plot(pl_lowDelta, pl_highDelta, pl_lowLambda, pl_highLambda, layout = (2,2), label = ["lowDelta" "highDelta" "lowLambda" "highLambda"])
end
plot_the_graph()


@show sol_false.eq.Iterobj_is[1,1].solution[2]
@show sol_false.eq.Iterobj_is[1,3].solution[2]
@show sol_false.eq.Iterobj_is[1].solution




@show sol_false.eq.Iterobj_is[2].solution
@show sol_false.eq.Iterobj_is[1].solution_index

@show sol_false.eqq.vFuncsNew.Γ[:, :,2]
@show sol_false.eqq.Iterobj_is[1]

@show policy_ex.lPolicy
@show policy_false.lPolicy

@show policy_ex.failure
@show policy_false.failure

@show sol.eq.vFuncsNew.qBond
@show sol_false.eq.vFuncs.qBond

function qBond_condiState2(params::Params{T,S}, Rl::T, il::S ,is::S, ib::S, iDelta::S) where {T<:Real,S<:Integer} 
        
        l, s, b, Delta = params.lGrid[il], params.sGrid[is], params.bGrid[ib], params.deltaGrid[iDelta]
        lambdaStar(l,s,b,Delta) = (Rl*(l+params.g*Delta)+(1+params.Rf)*s-Delta-b)/(Rl*l)
    
        function return_temp_underFailure(l,s,b,Delta,Lambda)::T
            b == 0 && return zero(T)
            term = Rl*((1-Lambda)*l+params.g*Delta)+(1+params.Rf)*s 
            num = max(min(b,max(zero(T),params.cF*term-Delta)),term-Delta-params.α*params.wr*Rl*(1-Lambda)*l)
            return num/b
        end

        @show lambdaStar_val = lambdaStar(l,s,b,Delta) # calculate lambdaStar for the given state (l,s,b,Delta)
        @show return_temp_underFailure.(l,s,b,Delta,params.lambdaGrid)
        return_temp_underFailure_val = ntuple(i -> return_temp_underFailure(l,s,b,Delta,params.lambdaGrid[i]), 3) # calculate return_temp_underFailure for the given state (l,s,b,Delta)
        @show return_temp_underFailure_val
        qBond_temp = zeros(T, 3)

        if lambdaStar_val < params.λL # Fail all the time  
            println("1")
            qBond_temp .= return_temp_underFailure_val
        elseif (params.λL <= lambdaStar_val) && (lambdaStar_val < params.λM)
            println("2")
            qBond_temp .= (one(T), return_temp_underFailure_val[2], return_temp_underFailure_val[3])
        elseif (params.λM <= lambdaStar_val) && (lambdaStar_val< params.λH)
            println("3")
            qBond_temp .= (one(T), one(T), return_temp_underFailure_val[3])
        else # No failure at all
            println("4")
            qBond_temp .= one(T)
        end

        x = params.β * (params.F * qBond_temp)
        @show x

        @show qBond_temp
        println("beta")
        @show params.β
        println("F")
        @show params.F

        
        return params.β * (params.F * qBond_temp)
end

x = qBond_condiState(params, 5.77882, 1, 5, 5, 3)

################################################################################################################################################



################################################################################################################################################

# 2. generate shock sequenc, simulate paths, and calculate moments
@time shocks = makeShockSequence(params, sol.eq.vFuncs, 100, 100, 10); # make shock sequence for simulation
paths = Initiate_Paths(shocks, params); # initiate paths for simulation

@time paths_ex = Simulate_paths(paths, policy_ex, sol.eq.vFuncs, params, sol.Rl_star, 10); # simulate paths and calculate moments
@time paths_false = Simulate_paths(paths, policy_false, sol_false.eq.vFuncs, params, sol_false.Rl_star, 10); # simulate paths and calculate moments
################################################################################################################################################

plot(paths_ex.failureSim[:,1,1]) # loan policy function for the first bank J at the first observation N
plot!(paths_ex.failureSim[:,1,3])
plot(paths_ex.lSim[:,2,1])
plot!(paths_false.lSim[:,2,1])
plot!(paths_ex.lSim[:,1,2])
plot(paths_ex.failureSim[:,1,1])
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
