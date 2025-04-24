pwd()



using Pkg
using LinearAlgebra
using Statistics
using QuantEcon
using JuMP
using Plots
using Ipopt
using Test
using Random, Statistics
using Interpolations
using Distributions


struct Params{T<:Real,S<:Integer}
    qd::T
    β::T
    Rf::T
    wr::T
    α::T
    ρ::T
    g::T
    ξ::T
    cF::T
    dBar::T
    σ::T
    τC::T
    z::T
    α1::T
    α2::T
    α3::T
    δL::T
    δM::T
    δH::T
    cM::T
    cO::T
    cL::T

    H  # transition matrix
    Γ  # transition matrix

    λL::T
    λM::T
    λH::T
    γ::T
    ϕ::T

    n_npts::S
    nGrid::Array{T,1}
    deltaGrid::Array{T,1}
    lambdaGrid::Array{T,1}
    l_npts::S
    s_npts::S
    b_npts::S
    lGrid::Array{T,1}
    sGrid::Array{T,1}   
    bGrid::Array{T,1}



    #---------------------------#
   # nstd::T
   # ygrid_npts::S
   # inflateEndPoints::Bool
   # bgrid_npts::S
   # bgrid_start::T
   # bgrid_stop::T
    
   # pay::T
   # ygrid_nodef::Array{T,1}
   # ygrid_def::Array{T,1}
   # bgrid::StepRangeLen{T,Base.TwicePrecision{T},Base.TwicePrecision{T}}
   # Pz::Array{T,2}
   # Ycons::Array{T,2}
   # DebtIncrease::Array{T,2}
end

struct vFuncs{T,S} where {T<:Real,S<:Integer}
    VF::Array{T,3} # value function, [Delta* Lambda* n]
    qBond::Array{T,5} # bank bond price schedule, [bPrime* lPrime* sPrime* Delta* Lambda]
    Rl::T # loan interest rate
    X::Array{T,4} # bank's failure decision 
end

struct vFuncsNew{T,S} where {T<:Real,S<:Integer} # storing updated outcome after VFI
    VF::Array{T,3} # value function, [Delta* Lambda* n]
    qBond::Array{T,5} # bank bond price schedule, [bPrime* lPrime* sPrime* Delta* Lambda]
    Rl::T # loan interest rate
    X::Array{T,4} # bank's failure decision 
    diffs::Array{T,3} # difference between VF VF_new
end

struct IterObj_i{T,S} where {T<:Real,S<:Integer} # objects given state (iDelta, iLambda), objects used in VFI_i
    EV::Array{T,3} # EV[l,s,b]
    G::Array{T,4} # G(l,s,b,n)
    solution::Array{T,1} # optimizer = [(l,s,b)]
    solution_index::Array{T,1} # optimizer = [(il,is,ib)]
    failure::Array{T,1} # failure decision = [(fail or not fail)]
end

v = zeros(3,3,n_npts);
pol_l = zeros(3,3,n_npts);
pol_s = zeros(3,3,n_npts);
pol_b = zeros(3,3,n_npts);

# flow utility of banks
psi(params::Params, d) = d >= 0 ? (d + params.dBar)^params.sigma  - params.dBar^params.sigma : 1 - e^(-d)

#=
function ModelParams_rating{T,S}(γc::T,σy::T,ρy::T,ξ::T,Rf::T,λ::T,z::T,β::T,d0::T,d1::T,σeps::T,ρeps::T,μeps::T,nstd::T,ygrid_npts::S,inflateEndPoints::Bool,bgrid_npts::S,bgrid_start::T,bgrid_stop::T,rJ::T,rF::T,x::T) where {T<:Real,S<:Integer}
        
    if bgrid_start!=0.0
        error("Grid for debt has to start at 0.0 !!")
    end

    # Define default costs
    ϕ(y) = max(0,d0*y+d1*y^2);

    pay = λ+(1.0-λ)*z;
    mc = tauchen(ygrid_npts, ρy, σy);
    Pz_temp = mc.p;
    ygrid_nodef_temp = mc.state_values;
    ygrid_nodef = exp.(ygrid_nodef_temp); # output in normal times
    Pz = Pz_temp';
 
    ygrid_def = ygrid_nodef-ϕ.(ygrid_nodef); #output in default
    bgrid = range(bgrid_start,stop=bgrid_stop,length=bgrid_npts);
    
    Ycons = ygrid_nodef'.-pay.*[bgrid;];
    DebtIncrease = (bgrid.-(1-λ).*bgrid');

    pam = ModelParams_rating{T,S}(γc,σy,ρy,ξ,Rf,λ,z,β,d0,d1,σeps,ρeps,μeps,nstd,ygrid_npts,inflateEndPoints,bgrid_npts,bgrid_start,bgrid_stop,pay,ygrid_nodef,ygrid_def,bgrid,Pz,Ycons,DebtIncrease,rJ,rF,x)
    return pam
    end
=#

function Initiate_vFunc(params::Params{T,S}) where Params{T<:Real,S<:Integer}
    VF = zeros(3, 3, params.n_npts); # expected value function
    qBond = zeros(3, 3, params.l_npts, params.s_npts, params.b_npts); # bank bond price schedule
    Rl = 0; # loan interest rate
    X = zeros(3, 3, params.n_npts, 3); # bank's failure decision 

    return vFuncs(VF, qBond, Rl, X)
end

###################################################

## function for iterating on the value function given Rl, w
function VFI(params::Params{T,S}, Rl::T, regime::F, maxiter::S, tol::T) where {T<:Real, S<:Integer, F<:Bool}

    function psi(d)
        return d >= 0 ? (d + params.dBar)^params.sigma  - params.dBar^params.sigma : 1 - e^(-d)
    end

    # 1. initiate value functions 
    vFuncs = Initiate_vFunc(params);
    vFuncsNew = Initiate_vFuncsNew(params); 
    Iterobj_is;

    Pz_trans 

    # 1-2. set the qBond schedule based on regime
    if regime == false # ordinary regime
        vFuncs.qBond[:, :, :, :, :] .= params.β; # set qBond to be 1 for all states 
    else # special regime
        qBond_specialRegime(params,vFuncs,Rl) # set qBond to be 1 for all states 
    end
    
    # 2. set numbers for iteration and set the iteration
    iter, maxdiffs, diffs = 1, one(T), zeros(S);

    # 3. start the iteration
    while iter <= maxiter && maxdiffs > tol

        # 3-1. set the value function for the next iteration

        # 3-2. The main iteration loop 
        Threads.@threads for iDelta in eachindex(params.deltaGrid) # 3: number of states for Delta
                 Threads.@threads for iLambda in eachindex(params.lambdaGrid)  # 3: number of states for Lambda
                    VFI_i(params, vFuncs, Rl, Iterobj_is[iDelta, iLmabda], iDelta, iLambda, regime); # VFI for a given exogenous state (iDelta, iLambda)
                end
        end

        Update_vFuncs_Diffs(vFuncs, vFuncsNew, params, Pz_trans, diffs); # after the optimization, calculate the differene 
        maxdiffs = maximum(diffs);

        if mod(iter, 200) == 0
            println("iter=", iter, ", maxdiffs=", maxdiffs); # report the iteration progress 
        end

        # 3-3. if the difference is not small enough, do the iteration again
        iter += 1;
    end

    # 4. finish the iteration and return the value function 
    return (vFuncs = vFuncs, Iterobj_is = Iterobj_is);  
end

## function for value function updating for a given exogenous state (iDelta, iLambda) and Rl
function VFI_i(params::Params{T,S}, vFuncs::vFuncs{T,S}, vFuncsNew, Rl::T, iterObj_i::IterObj_i{T,S}, iDelta::S, iLambda::S, regime::F) where {T<:Real, S<:Integer, F<:Bool}
    δ = params.deltaGrid[iDelta]; # get the state for Delta
    λ = params.lambdaGrid[iLambda]; # get the state for Lambda

    NAV(l::T,s::T,b::T,lambda::T)::T = Rl*( (1-lambda)*l + params.g*δ) + (1+params.Rf)*s - δ - b # net asset value of the bank, lambda here is lambda prime 
    tax(l::T,s::T,b::T,lambda::T)::T = params.τC * max(0, (Rl-1)*((1-lambda)*l +params.g*δ) + params.Rf*s - params.Rf *(δ + b)) # tax on the bank's asset value
    n_failure(l::T, lambda::T)::T = params.α * params.wr * Rl * (1-lambda)*l # next period asset conditional on bank failure 
    n_success(l::T,s::T,b::T,lambda::T)::T = NAV(l,s,b,lambda) - tax(l,s,b,lambda) # next period asset conditional on bank success 

    # 1) construct iterObj_i.EV[l,s,b] via gen_EV function 
    # 1-1) for a given (il, is, ib), construct V[delta prime, lambda prime] conditional on (delta prime, lambda prime)
    function inter_v_temp(params::Params{T,S}, vFuncs::vFuncs{T,S}, delta::T,lambda::T,n::T)::T
        nRange = params.nGrid # ex) 0:0.5:2
        lambdaRange = params.lambdaGrid
        deltaRange = params.deltaGrid
        vfVals = [vFuncs.VF[iDelta, iLambda, iN] for iN in eachindex(nRange), iLambda in eachindex(lambdaRange), iDelta in eachindex(deltaRange)] # [nDelta, nLambda], evaluated value functions at (delta prime, lambda prime) when choosing l,s,b

        itp = interpolate((nRange, lambdaRange, deltaRange), vfVals, BSpline(Cubic(Line(OnGrid())))) # Gridded(Linear())
        return itp(n,lambda,delta)
    end

    function gen_V_temp(l::T,s::T,b::T,params::Params{T,S},vFuncs::vFuncs{T,S},regime::F)::Array{T,2} # generate interpolated value of V (evaluated value functions) at (delta prime, lambda prime) when choosing l,s,b
        V_temp = Array{T}(undef, length(params.deltaGrid), length(params.lambdaGrid))  
        
        @inbounds for (iλ, λprime) in pairs(params.lambdaGrid)
            for (iδ, δprime) in pairs(params.deltaGrid)
                nav = NAV(l,s,b,λprime)
                failure = nav <= zero(T)

                if failure 
                        n_temp = n_failure(l,λprime)
                        V_temp[iδ, iλ] = regime ? zero(T) : (1 - params.ρ) * inter_v_temp(params, vFuncs, δprime, λprime, n_temp)
                else
                        n_temp = n_success(l,s,b,λprime)
                        V_temp[iδ, iλ] = inter_v_temp(params,vFuncs,δprime,λprime,n_temp) # interpolating v at (delta, lambda, n(iLambda))
                end
            end
            end

        return V_temp # [nDelta, nLambda] object 
    end

    # 1-2) multiply by transition matrix to get EV \in R: gen_EV_temp function 
    function gen_EV_temp(l::T,s::T,b::T,params::Params{T,S},vFuncs::vFuncs{T,S},regime::F)::T # incorporate failure array and conditional asset array to get ex-post asset array 
    
        V_temp = gen_V_temp(l,s,b,params,vFuncs, regime) # [nDelta, nLambda], evaluated value functions at (delta prime, lambda prime) when choosing l,s,b
        EV_conditional = dot(params.Markov_Delta[iDelta, :], V_temp * params.Markov_Lambda[:, iLambda]) # the return is a scalar 
        return EV_conditional 
    end

    # 1-3) iterate over all (il, is, ib) to get EV[L,S,B]: gen_EV function 
    function gen_EV!(params::Params{T,S},vFuncs::vFuncs{T,S},iterObj_i::IterObj_i{T,S},regime::F)

        EV = iterObj_i.EV 
        @inbounds for (il, l) in pairs(params.lGrid)
                    for (is, s) in pairs(params.sGrid)
                         for (ib, b) in pairs(params.bGrid)
                              EV[il,is,ib] = gen_EV_temp(l,s,b,params,vFuncs,regime) 
                        end
                     end
                end
        return iterObj_i
    end

    # 2) Construct iterObj_i.G[l,s,b,n] = flow utility[l,s,b,n] + β * EV[l,s,b]
    # 3) set the optimization model given state (iDelta, iLmambda, iN)
    # 4) store solution in iterObj_i.solution, iterObj_i.failure
    # 5) update vFuncs.VF
    function solve_bank_problem(vFuncs, iDelta::S, iLambda::S, iN::S) where {T<:Real, S<:Integer}
        model = Model(Ipopt.Optimizer)
        set_optimizer_attribute(model, "tol", 1e-8)
        l_min, l_max = first(params.lGrid), last(params.lGrid)
        s_min, s_max = first(params.sGrid), last(params.sGrid)
        b_min, b_max = first(params.bGrid), last(params.bGrid)
        @variable(model, l_min <= l <= l_max)
        @variable(model, s_min <= s <= s_max)
        @variable(model, b_min <= b <= b_max)

  #       div(l,s,b) = params.nGrid[iN] + (1-params.cL)*params.β*params.deltaGrid[iDelta] + vFuncs.qBond(l,s,b,iDelta,iLambda)*b - l - params.g * params.deltaGrid[iDelta] -s-params.cM*l^2 -params.cO 

        @NLobjective(model, Max, itp_G) # ** interpolated G as a function of (l,s,b), ordering of arguments?? 

        @NLconstraint(model, (l + params.g*params.deltaGrid[iDelta] +s - params.deltaGrid[iDelta]-b)/(params.wr*l) >= params.α) # constratint 1
        @NLconstraint(model, (l + params.g*params.deltaGrid[iDelta]) <= params.β*params.deltaGrid[iDelta]) # constratint 2
        optimize!(model)

        if is_solved_and_feasible(model) == false
            println("model is infeasible")
            break
        else
            return value(l), value(s), value(b) # returning the optimal values of l, s, b
        end
    end

    function update_VF_with_solution(sol,params::Params{T,S},vFuncs::vFuncs{T,S},regime)::T
        l,s,b = sol[1], sol[2], sol[3]
        ev = gen_EV_temp(l,s,b,params,vFuncs,regime) 
        div_val = div(l,s,b)
        return div_val + params.β* ev
    end

    # for each n in nGrid, find the optimal (l, s, b), store iterobj_i.solution, iterobj_i.failure, and update vFuncs.VF
    Threads.@threads for (iN, n) in pairs(params.nGrid) # for each n 
        
        div(l,s,b) = params.nGrid[iN] + (1-params.cL)*params.β*params.deltaGrid[iDelta] + vFuncs.qBond(l,s,b,iDelta,iLambda)*b - l - params.g * params.deltaGrid[iDelta] -s-params.cM*l^2 -params.cO 
        div = n .+ (1-params.cL)*params.β*δ .+ vFuncs.qBond(:, :, :, iDelta, iLambda)*B - L .- params.g*δ .- S - params.cM*L.^2 .- params.cO # (nL, nS, nB) matrix
        gen_EV!(params,vFuncs,iterObj_i,regime)
        G = div .+ params.β* iterObj_i.EV # G(l,s,b,n) = flow utility[l,s,b,n] + β * EV[l,s,b], (nL, nS, nB) matrix

        # G above will be feeded into solve_bank_problem for interpolating v at (l,s,b) and then solve for (l,s,b)
        sol = solve_bank_problem(vFuncs, iDelta, iLambda, iN); # find the solution, a vector with 3 elements (l,s,b)
        iterObj_i.solution[iN] .= sol; # store the solution in iterObj_i.solution
        iterObj_i.failure[iN] .= NAV(sol[1], sol[2], sol[3], params.lambdaGrid) .<= 0 ? 1 : 0; # store the failure decision in iterObj_i.failure
        vFuncsNew.VF[iDelta, iLambda, iN] = update_VF_with_solution(sol,params,vFuncs,regime)
    end
end    


# after VFI, get policy function for (l, s, b) and bank failure
function GetPolicyFuncs(vFuncs, params::Params{T,S}) where {T<:Real, S<:Integer}
    v
end

####### DONE FOR CONSTRUCTION ########  

function qBond_condiState(params::Params{T,S}, Rl::T, il::S ,is::S, ib::S, iDelta::S) where {T<:Real,S<:Integer} 
        
    l, s, b, Delta = params.lGrid[il], params.sGrid[is], params.bGrid[ib], params.deltaGrid[iDelta]

    lambdaStar(l,s,b,Delta) = (Rl*(l+params.g*Delta)+(1+params.Rf)*s-Delta-b)/(Rl*l)
    
    function return_temp_underFailure(l,s,b,Delta,Lambda)::T
        b == 0 && return zero(T)
        term = Rl*((1-Lambda)*l+params.g*Delta)+(1+params.Rf)*s 
        num = max(min(b,max(zero(T),params.cL*term-Delta)),term-Delta-params.α*params.wr*Rl*(1-Lambda)*l)
        return num/b
    end

    lambdaStar_val = lambdaStar(l,s,b,Delta) # calculate lambdaStar for the given state (l,s,b,Delta)
    return_temp_underFailure_val = ntuple(i -> return_temp_underFailure(l,s,b,Delta,params.lambdaGrid[i]), 3) # calculate return_temp_underFailure for the given state (l,s,b,Delta)
    qBond_temp = zeros(T, 3)

    if lambdaStar_val < params.λL # Fail all the time  
        return qBond_temp .= return_temp_underFailure_val
    elseif (params.λL <= lambdaStar_val) && (lambdaStar_val < params.λM)
        return qBond_temp .= (one(T), return_temp_underFailure_val[2], return_temp_underFailure_val[3])
    elseif (params.λM <= lambdaStar_val) && (lambdaStar_val< params.λH)
        return qBond_temp .= (one(T), one(T), return_temp_underFailure_val[3])
    else # No failure at all
        return qBond_temp .= one(T)
    end

    return params.β * (params.Markov_Lambda * qBond_temp)
end


function qBond_specialRegime(params::Params{T,S}, vFuncs::vFuncs{T,S}, Rl::T) where {T<:Real,S<:Integer}
    for il in eachindex(params.lGrid)
        for is in eachindex(params.sGrid)
            for ib in eachindex(params.bGrid)
                for iDelta in eachindex(params.deltaGrid)
                    vFuncs.qBond[il,is,ib,iDelta,:] .= qBond_condiState(params,Rl,il,is,ib,iDelta)
                end
            end
        end
    end
    return vFuncs
end


function Update_vFuncs_Diffs(vFuncs::vFuncs{T,S}, vFuncsNew::vFuncsNew{T,S}, params::Params{T,S}, diffs::Array{T,1}) where {T<:Real, S<:Integer}; # after the optimization, calculate the differene

     # 1. calculate the difference of VF_0 and VF_1 
     vFuncsNew.diffs = vFuncsNew.VF - vFuncs.VF; # difference between VF and VF_new
     diffs = norm(vFuncsNew.diffs, Inf); # infinity norm of the difference of VF

     # 2. update vFuncsNew objects to vFuncs
     vFuncs.VF .= vFuncsNew.VF; # update EVF
end


###################################### OLD CODE ################

 # function gen_n_temp(l::T,s::T,b::T,params::Params{T,S})::Vector{T} # incorporate failure array and conditional asset array to get ex-post asset array 
 #   n_temp = Vector{T}(undef, length(params.lambdaGrid)) # temporary variable for the next period asset value
 #
 #   @inbounds for (iλ, λprime) in pairs(params.lambdaGrid)
 #       nav = NAV(l,s,b,λprime)
 #       n_temp[iλ] = nav < zero(T) ? n_failure(l, λprime) : n_success(l,s,b, λprime) 
 #   end
 #   return n_temp

