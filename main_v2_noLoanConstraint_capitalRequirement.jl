
pwd()


# include("C:\\Users\\master\\Desktop\\(2025) banking regulation\\BankResolution\\main.jl")
# include("C:\\Users\\master\\Desktop\\(2025) banking regulation\\BankResolution\\parameters.jl")
# include("/Users/juyoungyang_kdi/BankResolution/parameters.jl")
# include("/Users/juyoungyang_kdi/BankResolution/main.jl")

using JuMP, Ipopt, Interpolations, MathOptInterface, LinearAlgebra, ForwardDiff, Plots, Ipopt, Roots, Statistics
using Profile, ProfileView, QuantEcon, StatsBase, Distributions, Random
using Sobol, Distributions, Distributed
const MOI = MathOptInterface

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
        δL::T
        δM::T
        δH::T
        cM::T
        cO::T
        cL::T
        ϵ::T
        E::T

        H  # transition matrix for delta: [delta, deltaPrime] so that sum(H[1,:]) = 1 
        F  # transition matrix for lambda: [lambdaPrime, lambda] so that sum(G[:,1]) = 1
        M::T # bank mass, policy tool 

        λL::T
        λM::T
        λH::T
        γ::T
        ϕ::T

        sigHat::T ## curvature for div, especially for stock issuance 
        lconstr::T ## deposit to loan constraint

        deltaGrid::Array{T,1}
        lambdaGrid::Array{T,1}
        nGrid::Array{T,1}
        lGrid::Array{T,1}
        sGrid::Array{T,1}   
        bGrid::Array{T,1}
     #---------------------------#
     # v = zeros(3,3,n_npts);
     # pol_l = zeros(3,3,n_npts);
     # pol_s = zeros(3,3,n_npts);
     # pol_b = zeros(3,3,n_npts);


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

struct VFuncs{T<:Real,S<:Integer}
        VF::Array{T,3} # value function, [Delta* Lambda* n]
        qBond::Array{T,5} # bank bond price schedule, [bPrime* lPrime* sPrime* Delta* Lambda]
        X::Array{T,4} # bank's failure decision 
        Γ::Array{T,3} # bank mass for each state, [Delta* Lambda* n]
end

struct VFuncsNew{T<:Real,S<:Integer} # storing updated outcome after VFI
        VF::Array{T,3} # value function, [Delta* Lambda* n]
        qBond::Array{T,5} # bank bond price schedule, [bPrime* lPrime* sPrime* Delta* Lambda]
        X::Array{T,4} # bank's failure decision 
        diffs::Array{T,3} # difference between VF VF_new
        Γ::Array{T,3} # bank mass for each state, [Delta* Lambda* n]
end

struct IterObj_i{T<:Real,S<:Integer} # objects given state (iDelta, iLambda), objects used in VFI_i
        EV::Array{T,3} # EV[l,s,b]
        G::Array{T,4} # G(l,s,b,n)
        solution::Array{Array{T,1},1} # optimizer = [(l,s,b)]
        solution_index::Array{Array{T,1},1} # optimizer = [(il,is,ib) for each n] 
        failure::Array{Array{T,1},1} # failure decision = [(fail or not fail vector for each lambda prime) for each n]
end

function Initiate_Params(qd::T,β::T,Rf::T,wr::T,α::T,ρ::T,g::T,ξ::T,cF::T,dBar::T,σ::T,τC::T,z::T,δL::T,δM::T,δH::T,cM::T,cO::T,cL::T,ϵ::T,E::T,H::Array{T,2},F::Array{T,2},M::T,λL::T,λM::T,λH::T,γ::T,ϕ::T,sigHat::T,lconstr::T,n_start::T,n_npts::S,n_stop::T,l_start::T,l_npts::S,l_stop::T,s_start::T,s_npts::S,s_stop::T,b_start::T,b_npts::S,b_stop::T) where {T<:Real,S<:Integer}
        deltaGrid = [δL,δM,δH] # Define a Tuple, immutable 
        lambdaGrid = [λL,λM,λH] # regular array, mutable
        nGrid = range(n_start,stop=n_stop,length=n_npts)
    
        lGrid = range(l_start,stop=l_stop,length=l_npts)
        sGrid = range(s_start,stop=s_stop,length=s_npts)
        bGrid = range(b_start,stop=b_stop,length=b_npts)
    
        pam = Params{T,S}(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,δL,δM,δH,cM,cO,cL,ϵ,E,H,F,M,λL,λM,λH,γ,ϕ,sigHat,lconstr,deltaGrid,lambdaGrid,nGrid,lGrid,sGrid,bGrid)
        return pam
end

function Initiate_vFunc(params::Params{T,S}) where {T<:Real,S<:Integer}
        VF = ones(3, 3, length(params.nGrid)); # expected value function
        qBond = zeros(length(params.lGrid), length(params.sGrid), length(params.bGrid), 3, 3); # bank bond price schedule
        X = zeros(3, 3, length(params.nGrid), 3); # bank's failure decision 
        Γ = ones(3, 3, length(params.nGrid))/sum(ones(3, 3, length(params.nGrid))); # bank mass for each state, [Delta* Lambda* n]
 
        return VFuncs{T,S}(VF, qBond, X, Γ)
end

function Initiate_vFuncNew(params::Params{T,S}) where {T<:Real,S<:Integer}
        VF = zeros(3, 3, length(params.nGrid)); # expected value function
        qBond = zeros(length(params.lGrid), length(params.sGrid), length(params.bGrid), 3, 3); # bank bond price schedule
        X = zeros(3, 3, length(params.nGrid), 3); # bank's failure decision 
        diffs = zeros(3, 3, length(params.nGrid)); # difference between VF and VF_new
        Γ = ones(3, 3, length(params.nGrid))/sum(ones(3, 3, length(params.nGrid))); # bank mass for each state, [Delta* Lambda* n]

        return VFuncsNew{T,S}(VF, qBond, X, diffs, Γ)
end

function Initiate_IterObj_i(params::Params{T,S}) where {T<:Real,S<:Integer}
        EV = zeros(length(params.lGrid), length(params.sGrid), length(params.bGrid)); # expected value function conditional on a choice of (l,s,b)
        G = zeros(length(params.lGrid), length(params.sGrid), length(params.bGrid), length(params.nGrid)); # bank bond price schedule

        nN = length(params.nGrid)
        solution        = Vector{Vector{T}}(undef, nN)
        solution_index  = Vector{Vector{T}}(undef, nN)
        failure = Vector{Vector{T}}(undef, nN)
        for i in 1:nN
            solution[i] = zeros(T, 3)
            solution_index[i] = zeros(T, 3)
            failure[i] = zeros(T, 3)
        end
        return IterObj_i{T,S}(EV, G, solution, solution_index, failure)
end

function Initiate_MatrixIterObj_i(params::Params{T,S}) where {T<:Real,S<:Integer}
        IterObj_is = Matrix{IterObj_i{T,S}}(undef, length(params.deltaGrid), length(params.lambdaGrid));

        for iDelta in 1:length(params.deltaGrid)
            for iLambda in 1:length(params.lambdaGrid)
                IterObj_is[iDelta, iLambda] = Initiate_IterObj_i(params);
            end
        end

        return IterObj_is; # return a vector {IterObj_iy[iy]} s.t. a vector of (delta, lambda)-contingent struct IterObj_iy
end


#######################################################################################
psi(params::Params, d) = d >= 0 ? (d + params.dBar)^params.σ  - params.dBar^params.σ : 1 - exp(-d)^params.sigHat; # 1.2
# psi(params::Params, d) = d >= 0 ? (d + params.dBar)^params.σ  - params.dBar^params.σ : 0.0;
#=
function inter_v_temp(params::Params{T,S}, vFuncs::VFuncs{T,S}, delta::T,lambda::T,n::T)::T  where {T<:Real, S<:Integer}
        nRange = range(params.nGrid[1], stop = params.nGrid[end], length = length(params.nGrid)) # ex) 0:0.5:2
        lambdaRange = params.lambdaGrid
        deltaRange = params.deltaGrid
        vfVals = [vFuncs.VF[iDelta, iLambda, iN] for iN in eachindex(nRange), iLambda in eachindex(lambdaRange), iDelta in eachindex(deltaRange)] # [nDelta, nLambda], evaluated value functions at (delta prime, lambda prime) when choosing l,s,b
        itp = interpolate((nRange, lambdaRange, deltaRange), vfVals, Gridded(Linear())) # , BSpline(Cubic(Line(OnGrid())))
        itp_ext = extrapolate(itp, Line())

        return itp_ext(n, lambda, delta)
end

function gen_V_temp(l::T,s::T,b::T,params::Params{T,S},vFuncs::VFuncs{T,S},regime::F)::Array{T,2}  where {T<:Real, S<:Integer, F<:Bool} # generate interpolated value of V (evaluated value functions) at (delta prime, lambda prime) when choosing l,s,b
        V_temp = Array{T}(undef, length(params.deltaGrid), length(params.lambdaGrid))  
        
        @inbounds for (iλ, λprime) in pairs(params.lambdaGrid)
            for (iδ, δprime) in pairs(params.deltaGrid)
                nav = NAV(l,s,b,λprime)
                failure = nav <= zero(T)
                # @show iδ, iλ

                if failure 
                   #  println("failure at (iδ, iλ): ", iδ, iλ)
                        n_temp = n_failure(l,λprime)
                        V_temp[iδ, iλ] = regime ? zero(T) : params.ρ * inter_v_temp(params, vFuncs, δprime, λprime, n_temp)
                else
                   #  println("success at (iδ, iλ): ", iδ, iλ)
                        n_temp = n_success(l,s,b,λprime)
                        V_temp[iδ, iλ] = inter_v_temp(params,vFuncs,δprime,λprime,n_temp) # interpolating v at (delta, lambda, n(iLambda))
                end
            end
            end

        return V_temp # [nDelta, nLambda] object 
end

    # 1-2) multiply by transition matrix to get EV \in R: gen_EV_temp function 
function gen_EV_temp(l::T,s::T,b::T,params::Params{T,S},vFuncs::VFuncs{T,S},regime::F)::T  where {T<:Real, S<:Integer, F<:Bool} # incorporate failure array and conditional asset array to get ex-post asset array 
    
        V_temp = gen_V_temp(l,s,b,params,vFuncs, regime) # [nDelta, nLambda], evaluated value functions at (delta prime, lambda prime) when choosing l,s,b
        @show size(V_temp), size(params.H[iDelta, :]), size( params.F[iLambda,:]), iDelta, iLambda
        EV_conditional = dot(params.H[iDelta, :], V_temp * params.F[iLambda,:]) # the return is a scalar 
        return EV_conditional 
end

    # 1-3) iterate over all (il, is, ib) to get EV[L,S,B]: gen_EV function 
function gen_EV!(params::Params{T,S},vFuncs::VFuncs{T,S},iterObj_i::IterObj_i{T,S},regime::F) where {T<:Real, S<:Integer, F<:Bool}

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

function NAV(l::T,s::T,b::T,lambda::T)::T where {T<:Real} 
        return Rl*( (1-lambda)*l + params.g*δ) + (1+params.Rf)*s - δ - b # net asset value of the bank, lambda here is lambda prime 
end

function tax(l::T,s::T,b::T,lambda::T)::T where {T<:Real} 
        return params.τC * max(0, (Rl-1)*((1-lambda)*l +params.g*δ) + params.Rf*s - params.Rf *(δ + b)) # tax on the bank's asset value
end

function n_failure(l::T, lambda::T)::T where {T<:Real} 
        return params.α * params.wr * Rl * (1-lambda)*l # next period asset conditional on bank failure 
end

function n_success(l::T,s::T,b::T,lambda::T)::T where {T<:Real} 
        return NAV(l,s,b,lambda) - tax(l,s,b,lambda) # next period asset conditional on bank success 
end


function gen_V_temp(l::T,s::T,b::T,params::Params{T,S},vFuncs::VFuncs{T,S},regime::F)::Array{T,2}  where {T<:Real, S<:Integer, F<:Bool} # generate interpolated value of V (evaluated value functions) at (delta prime, lambda prime) when choosing l,s,b
        V_temp = Array{T}(undef, length(params.deltaGrid), length(params.lambdaGrid))  
        @show l,s,b
        @inbounds for (iλ, λprime) in pairs(params.lambdaGrid)
            for (iδ, δprime) in pairs(params.deltaGrid)
                @show nav = NAV(l,s,b,λprime)
                failure = nav <= zero(T)

                if failure 
                    println("failure at (iδ, iλ): ", iδ, iλ)
                        @show n_temp = n_failure(l,λprime)
                        @show V_temp[iδ, iλ] = regime ? zero(T) : (1 - params.ρ) * inter_v_temp(params, vFuncs, δprime, λprime, n_temp)
                else
                    println("success at (iδ, iλ): ", iδ, iλ)
                        @show n_temp = n_success(l,s,b,λprime)
                        @show V_temp[iδ, iλ] = inter_v_temp(params,vFuncs,δprime,λprime,n_temp) # interpolating v at (delta, lambda, n(iLambda))
                end
            end
            end

        return V_temp # [nDelta, nLambda] object 
end
=#

function VFI_i(params::Params{T,S}, vFuncs::VFuncs{T,S}, vFuncsNew::VFuncsNew{T,S}, Rl::T, iterObj_i::IterObj_i{T,S}, iDelta::S, iLambda::S, regime::F) where {T<:Real, S<:Integer, F<:Bool}
        δ = params.deltaGrid[iDelta]; # get the state for Delta
        λ = params.lambdaGrid[iLambda]; # get the state for Lambda

        NAV(l::T,s::T,b::T,lambda::T)::T = Rl*( (1-lambda)*l + params.g*δ) + (1+params.Rf)*s - δ - b # net asset value of the bank, lambda here is lambda prime 
        tax(l::T,s::T,b::T,lambda::T)::T = params.τC * max(0, (Rl-1)*((1-lambda)*l +params.g*δ) + params.Rf*s - params.Rf *b - params.Rf*0.8*δ) # tax on the bank's asset value
        n_failure(l::T, lambda::T)::T = params.α * params.wr * Rl * (1-lambda)*l # next period asset conditional on bank failure 
        n_success(l::T,s::T,b::T,lambda::T)::T = NAV(l,s,b,lambda) - tax(l,s,b,lambda) # next period asset conditional on bank success 
        # div_func(l::T,s::T,b::T,n::T)::T = n + (1-params.cL)*params.β*δ + vFuncs.qBond[l,s,b,iDelta,iLambda]*b - l - params.g * δ -s-params.cM*l^2 -params.cO 
        # 자기자본규제 미충족시 은행이 실패한다고 가정할 때
        capitalRequire(l::T,s::T,b::T,lambda::T)::T = (1-params.α*params.wr)*(1-lambda)*l + s - b - (1-params.g)*params.deltaGrid[iDelta] # bank fails if negative

        function div_func_with_solution(params::Params{T,S}, vFuncs::VFuncs{T,S},iDelta::S,iLambda::S,l::T,s::T,b::T,n::T)::T
            # need to interpolate qBond at the optimum level (l,s,b) given state (iDelta, iLambda) with n 

            lRange = range(first(params.lGrid), last(params.lGrid), length(params.lGrid))
            sRange = range(first(params.sGrid), last(params.sGrid), length(params.sGrid))
            bRange = range(first(params.bGrid), last(params.bGrid), length(params.bGrid))
            qBond_Matrix = @view vFuncs.qBond[:,:,: , iDelta, iLambda]

            itp = interpolate(qBond_Matrix, BSpline(Quadratic(Line(OnGrid()))))
            itp = Interpolations.scale(itp, lRange, sRange, bRange) # scale the interpolation to the grid ranges
            itp_ext = extrapolate(itp, Flat())
            qBond_interpolated = itp_ext(l, s, b) # qBond[l,s,b,iDelta,iLambda] interpolated at (l,s,b)

            div_interpolated = n + (1-params.cL)*params.β*δ + qBond_interpolated*b - l - params.g * δ -s-params.cM*l^2 -params.cO 
            return div_interpolated 
        end 

        # 1) construct iterObj_i.EV[l,s,b] via gen_EV function 
        # 1-1) for a given (il, is, ib), construct V[delta prime, lambda prime] conditional on (delta prime, lambda prime)
        function inter_v_temp(params::Params{T,S}, vFuncs::VFuncs{T,S}, delta::T,lambda::T,n::T)::T
            nRange = range(params.nGrid[1], stop = params.nGrid[end], length = length(params.nGrid)) # ex) 0:0.5:2
            lambdaRange = params.lambdaGrid
            deltaRange = params.deltaGrid
            vfVals = [vFuncs.VF[iDelta, iLambda, iN] for iN in eachindex(nRange), iLambda in eachindex(lambdaRange), iDelta in eachindex(deltaRange)] # [nDelta, nLambda], evaluated value functions at (delta prime, lambda prime) when choosing l,s,b

            itp = interpolate((nRange, lambdaRange, deltaRange), vfVals, Gridded(Linear())) # , BSpline(Cubic(Line(OnGrid())))
            itp_ext = extrapolate(itp, Line())
            return itp_ext(n, lambda, delta)
        end

        function gen_V_temp(l::T,s::T,b::T,params::Params{T,S},vFuncs::VFuncs{T,S},regime::F)::Array{T,2} # generate interpolated value of V (evaluated value functions) at (delta prime, lambda prime) when choosing l,s,b
            V_temp = Array{T}(undef, length(params.deltaGrid), length(params.lambdaGrid))  
        
            @inbounds for (iλ, λprime) in pairs(params.lambdaGrid)
                for (iδ, δprime) in pairs(params.deltaGrid)
                    # nav = NAV(l,s,b,λprime)
                    # failure = nav <= zero(T)
                    # 자기자본규제 충족여부가 은행의 실패여부를 결정할 때 
                    capitalRe = capitalRequire(l, s, b, λprime);
                    failure = capitalRe <= zero(T)

                    if failure 
                        n_temp = n_failure(l,λprime)
                        V_temp[iδ, iλ] = regime ? zero(T) : params.ρ * inter_v_temp(params, vFuncs, δprime, λprime, n_temp) # ρ being failure panelty
                    else
                        n_temp = n_success(l,s,b,λprime)
                        V_temp[iδ, iλ] = inter_v_temp(params,vFuncs,δprime,λprime,n_temp) # interpolating v at (delta, lambda, n(iLambda))
                    end
                end
            end

            return V_temp # [nDelta, nLambda] object 
        end

        # 1-2) multiply by transition matrix to get EV \in R: gen_EV_temp function 
        function gen_EV_temp(l::T,s::T,b::T,params::Params{T,S},vFuncs::VFuncs{T,S},regime::F)::T # incorporate failure array and conditional asset array to get ex-post asset array 
            V_temp = gen_V_temp(l,s,b,params,vFuncs, regime) # [nDelta, nLambda], evaluated value functions at (delta prime, lambda prime) when choosing l,s,b
            EV_conditional = dot(params.H[iDelta, :], V_temp * params.F[iLambda, :]) # the return is a scalar 
            return EV_conditional # expected value 
        end

        # 1-3) iterate over all (il, is, ib) to get EV[L,S,B]: gen_EV function 
        function gen_EV!(params::Params{T,S},vFuncs::VFuncs{T,S},iterObj_i::IterObj_i{T,S},regime::F)
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
        function solve_bank_problem(params::Params{T,S},vFuncs::VFuncs{T,S},iDelta::S,iLambda::S,iN::S,G::Array{T,3}; warm_start::Union{Nothing,NTuple{3,Float64}}=nothing) where {T<:Real, S<:Integer}

            lRange = range(first(params.lGrid), last(params.lGrid), length(params.lGrid))
            sRange = range(first(params.sGrid), last(params.sGrid), length(params.sGrid))
            bRange = range(first(params.bGrid), last(params.bGrid), length(params.bGrid))
    
            Gmin = minimum(G)
            Gmax = maximum(G)
            Gnorm = (G .- Gmin) ./ max(Gmax - Gmin, eps(T)) # scaling 
            G_itp_rev = interpolate(Gnorm, BSpline(Linear())) # , BSpline(Cubic(Line(OnGrid())))
            G_itp_rev = Interpolations.scale(G_itp_rev, lRange, sRange, bRange)
            G_itp_ext_rev = extrapolate(G_itp_rev, Flat()) # , Line()
            G_interp_rev(l, s, b) = G_itp_ext_rev(l, s, b)
    
            model = Model(Ipopt.Optimizer)
            set_optimizer_attribute(model, "tol", 1e-4)
            set_optimizer_attribute(model, "acceptable_tol", 5e-3)
            set_optimizer_attribute(model, "acceptable_constr_viol_tol", 1e-4)
            set_optimizer_attribute(model, "acceptable_dual_inf_tol", 1e-4)
            set_optimizer_attribute(model, "acceptable_compl_inf_tol", 1e-4)
            set_optimizer_attribute(model, "acceptable_obj_change_tol", 1e-6)
            set_optimizer_attribute(model, "max_iter", 5000)
            set_optimizer_attribute(model, "acceptable_iter", 5)
            set_optimizer_attribute(model, "print_level", 0)
            set_optimizer_attribute(model, "hessian_approximation", "limited-memory")
            set_optimizer_attribute(model, "mu_strategy", "adaptive")
            set_optimizer_attribute(model, "nlp_scaling_method", "gradient-based")
            set_optimizer_attribute(model, "watchdog_shortened_iter_trigger", 8)  # 기본 10
            set_optimizer_attribute(model, "watchdog_trial_iter_max", 3) 
           # set_optimizer_attribute(model, "max_wall_time", 150.0) # 시간 방어막: iteration_limit 대신 acceptable로 떨어지게 유도
            set_optimizer_attribute(model, "print_options_documentation", "yes")

            l_min, l_max = first(params.lGrid), last(params.deltaGrid) # (last(params.lGrid)+first(params.lGrid))/2, last(params.lGrid)
            s_min, s_max = first(params.sGrid), last(params.sGrid)
            b_min, b_max = first(params.bGrid), last(params.bGrid)
            ## alternative constraints: 예대율 규제
            # l_cap = min(l_max, (params.β - params.g) * params.deltaGrid[iDelta]) # l_start = (l_max+l_min)/2
            l_cap = min(l_max, params.lconstr * params.deltaGrid[iDelta]) # l_start = (l_max+l_min)/2

            # b_cap = min(b_max, (1 - params.α * params.wr) * l_max + s_max - (1 - params.g) * params.deltaGrid[iDelta]) # b_start = (b_max+b_min)/2

            # start: warm-start 있으면 사용, 없으면 중앙값
           #  l0 = warm_start === nothing ? (l_min + l_max)/2 : warm_start[1] # 예대율 규제가 없을 때
             l0 = warm_start === nothing ? (l_min + l_cap)/2 : warm_start[1] # 예대율 규제가 들어갈 때
            s0 = warm_start === nothing ? (s_min + s_max)/2 : warm_start[2]
           #  b0 = warm_start === nothing ? (b_min + b_cap)/2 : warm_start[3]
            b0 = warm_start === nothing ? (b_min + b_max)/2 : warm_start[3]

            
            function G_interp_rev_safe(l, s, b)
               l_safe = clamp(l, l_min, l_cap) # 예대율 규제가 들어갈 때
             #   l_safe = clamp(l, l_min, l_max) # 예대율 규제가 없을 때
                s_safe = clamp(s, s_min, s_max)
                b_safe = clamp(b, b_min, b_max)
                return G_itp_ext_rev(l_safe, s_safe, b_safe)
            end

            JuMP.register(model, :G_interp_rev, 3, G_interp_rev_safe; autodiff = true)

            # println("c1 residual = ", (1-params.α*params.wr)*l_start + s_start - b_start - (1-params.g)*params.deltaGrid[iDelta])
            # println("c2 residual = ", l_start - (params.β - params.g)*params.deltaGrid[iDelta])
          #  @variable(model, l_min <= l <= l_max, start = clamp(l0, l_min, l_max)) # 예대율 규제가 없을 때
            @variable(model, l_min <= l <= l_cap, start = clamp(l0, l_min, l_cap)) # 예대율 규제가 들어갈 때
            @variable(model, s_min <= s <= s_max, start = clamp(s0, s_min, s_max))
            @variable(model, b_min <= b <= b_max, start = clamp(b0, b_min, b_max))
            # println("G_interp at start = ", G_interp_rev(l_start, s_start, b_start))
            # grad = ForwardDiff.gradient(u -> G_interp_rev(u[1], u[2], u[3]), [l_start, s_start, b_start])
            # println("gradient at start = ", grad)
    
            if any(isnan, G) || any(isinf, G)
                @warn "NaN or Inf detected in G at iDelta=$iDelta, iLambda=$iLambda, iN=$iN"
            end

            if !all(isfinite, [params.α, params.β, Rl])
                @warn "Non-finite parameter detected: α=$(params.α), β=$(params.β), Rl=$Rl"
            end

            @NLobjective(model, Max, G_interp_rev(l, s, b))
            @constraint(model, (1-0.13*params.wr)*l + s - b >= (1-params.g)*params.deltaGrid[iDelta]) # 자기자본규제

            # @constraint(model, (1-params.α*params.wr)*l + s - b >= (1-params.g)*params.deltaGrid[iDelta]) # 자기자본규제
           # @constraint(model, l <= (params.β - params.g)*params.deltaGrid[iDelta]) # constratint 2
    
            try
                optimize!(model)
                catch e
                @warn "Solver threw an error: $e. 그래도 상태를 확인합니다."
            end
    
            # @show is_solved_and_feasible(model)
            stat = termination_status(model)
            # println("termination_status = ", stat)
    
            if stat == MathOptInterface.OPTIMAL || stat == MathOptInterface.LOCALLY_SOLVED || stat == MathOptInterface.ALMOST_LOCALLY_SOLVED
                return [value(l), value(s), value(b)] # optimizer
            elseif stat == MathOptInterface.ITERATION_LIMIT && primal_status(model) == MathOptInterface.FEASIBLE_POINT
            # 안전장치: 제약 위반이 작으면 '사실상 해'로 채택 (임계치는 필요시 조정)
            # 간단 체크 예시 (정밀 필요시 제약 잔차 직접 계산)
                return (value(l), value(s), value(b))
            elseif termination_status(model) == MathOptInterface.NUMERICAL_ERROR
                @warn "Numerical error detected."
                @info "Values at failure:" l=value(l) s=value(s) b=value(b)
                @show l_min, l_max, s_min, s_max, b_min, b_max
                @info "Objective value:" obj=objective_value(model)
            else
                error("Optimization did not converge: status=$stat, primal=$(primal_status(model))")
            end
        end

        function update_VF_with_solution(sol,params::Params{T,S},vFuncs::VFuncs{T,S},regime,n)::T
            l,s,b = sol[1], sol[2], sol[3]
            ev = gen_EV_temp(l,s,b,params,vFuncs,regime) 
            div_val = div_func_with_solution(params, vFuncs, iDelta, iLambda, l,s,b,n)
            return div_val + params.β* ev
        end
    
        # for each n in nGrid, find the optimal (l, s, b), store iterobj_i.solution, iterobj_i.failure, and update vFuncs.VF
        for (iN, n) in pairs(params.nGrid) # for each n / Threads.@threads 
    
            L = reshape(params.lGrid, :, 1, 1)    # nL × 1 × 1       
            SS = reshape(params.sGrid, 1, :, 1)   # 1  × nS × 1
            B = reshape(params.bGrid, 1, 1, :)    # 1  × 1  × nB
            QB = @view vFuncs.qBond[:, :, :, iDelta, iLambda]  # nL × nS × nB
            divv = n .+ (1 - params.cL) * params.β * δ .+ QB .* B .- L .- params.g * δ .- SS .- params.cM .* L.^2 .- params.cO
            
            # div(l,s,b) = n + (1-params.cL)*params.β*δ + vFuncs.qBond[l,s,b,iDelta,iLambda]*b - l - params.g * δ -s-params.cM*l^2 -params.cO 
            # div = n .+ (1-params.cL)*params.β*δ .+ vFuncs.qBond[:, :, :, iDelta, iLambda]*b - L .- params.g*δ .- S - params.cM*L.^2 .- params.cO # (nL, nS, nB) matrix
            gen_EV!(params,vFuncs,iterObj_i,regime)
            G = psi.([params], divv) .+ params.β * iterObj_i.EV # G(l,s,b,n) = flow utility[l,s,b,n] + β * EV[l,s,b], (nL, nS, nB) matrix with fixed n 
    
            #### psi ####
    
            # G above will be feeded into solve_bank_problem for interpolating v at (l,s,b) and then solve for (l,s,b)
            sol = solve_bank_problem(params, vFuncs, iDelta, iLambda, iN, G); # find the solution, a vector with 3 elements (l,s,b)
            iterObj_i.solution[iN] .= sol; # store the solution in iterObj_i.solution

            # iterObj_i.failure[iN] = [NAV(sol[1], sol[2], sol[3], lamPrime) < 0 for lamPrime in params.lambdaGrid] # if true, bank fails 
            # 자기자본규제가 실패여부를 결정지을 때 
            iterObj_i.failure[iN] = [capitalRequire(sol[1], sol[2], sol[3], lamPrime) < 0 for lamPrime in params.lambdaGrid] # if true, bank fails 

            vFuncsNew.VF[iDelta, iLambda, iN] = update_VF_with_solution(sol,params,vFuncs,regime,n)
        end
end

function VFI(params::Params{T,S}, Rl::T, regime::F, maxiter::S, tol::T) where {T<:Real, S<:Integer, F<:Bool}

        # 1. initiate value functions 
        vFuncs = Initiate_vFunc(params);
        vFuncsNew = Initiate_vFuncNew(params); 
        Iterobj_is = Initiate_MatrixIterObj_i(params); 

        # 1-2. set the qBond schedule based on regime
        if regime == false # ordinary regime
            vFuncs.qBond .= params.β; # bank bond is risk-free 
        end
    
        ## 2. set numbers for iteration and set the iteration  
        iter, maxdiffs, diffs = 1, one(T), zeros(T);

        if iter == 1
            println("VFI 시작 시 VF 평균값: ", mean(vFuncs.VF))
        end
        # 3. start the iteration
        while iter <= maxiter && maxdiffs > tol

            # 3-1. set the value function for the next iteration

            # 3-2. The main iteration loop 
            Threads.@threads for iDelta in 1:length(params.deltaGrid) # 3: number of states for Delta
                 Threads.@threads for iLambda in 1:length(params.lambdaGrid)  # 3: number of states for Lambda
                    VFI_i(params, vFuncs, vFuncsNew, Rl, Iterobj_is[iDelta, iLambda], iDelta, iLambda, regime); # VFI for a given exogenous state (iDelta, iLambda)
                end
            end

            if regime == false # ordinary regime
                vFuncsNew.qBond .= params.β; # bank bond is risk-free 
            else # special regime
                qBond_specialRegime(params,vFuncsNew,Rl) # set qBond to be 1 for all states 
            end

            diffs = Update_vFuncs_Diffs(vFuncs, vFuncsNew, params); # after the optimization, calculate the differene 
            maxdiffs = maximum(diffs);

            if mod(iter, 200) == 0
            println("VFI: iter=", iter, ", maxdiffs=", maxdiffs); # report the iteration progress 
            end

            # 3-3. if the difference is not small enough, do the iteration again
            iter += 1;
        end

        # 4. finish the iteration and return the value function 
        println("VFI 종료 시 VF 평균값: ", mean(vFuncsNew.VF))
        return (vFuncs = vFuncs, vFuncsNew = vFuncsNew, Iterobj_is = Iterobj_is);  
end

function VFI!(params::Params{T,S}, Rl::T, regime::F, maxiter::S, tol::T) where {T<:Real, S<:Integer, F<:Bool}

        iter, maxdiffs, diffs = 1, one(T), zeros(T);

        if iter == 1
            println("VFI 시작 시 VF 평균값: ", mean(vFuncs.VF))
        end

        # start the iteration
        while iter <= maxiter && maxdiffs > tol

            Threads.@threads for iDelta in 1:length(params.deltaGrid) # 3: number of states for Delta
                for iLambda in 1:length(params.lambdaGrid)
                    VFI_i(params, vFuncs, vFuncsNew, Rl, Iterobj_is[iDelta, iLambda], iDelta, iLambda, regime); # VFI for a given exogenous state (iDelta, iLambda)
                end
            end

            if regime == false # ordinary regime
                vFuncsNew.qBond .= params.β; # bank bond is risk-free 
            else # special regime
                qBond_specialRegime(params,vFuncsNew,Rl) # set qBond to be 1 for all states 
            end

            # 감쇠(damping) 업데이트로 진동 억제
            @. vFuncs.VF = (1 - damping) * vFuncs.VF + damping * vFuncsNew.VF
            @. vFuncs.qBond = (1 - damping) * vFuncs.qBond + damping * vFuncsNew.qBond
            @. vFuncs.X = (1 - damping) * vFuncs.X + damping * vFuncsNew.X

            # 수렴도 계산
            diffs = abs.(vFuncs.VF .- vFuncsNew.VF)
            maxdiffs = maximum(diffs)

            if mod(iter, 200) == 0
                println("VFI!: iter=$iter, maxdiffs=$maxdiffs")
            end

            iter += 1
        end

    println("VFI! 종료 시 VF 평균값: ", mean(vFuncs.VF))
    return nothing
end

function Update_vFuncs_Diffs(vFuncs::VFuncs{T,S}, vFuncsNew::VFuncsNew{T,S}, params::Params{T,S}) where {T<:Real, S<:Integer}; # after the optimization, calculate the differene

        # 1. calculate the difference of VF_0 and VF_1 
        vFuncsNew.diffs .= vFuncsNew.VF - vFuncs.VF; # difference between VF and VF_new
        # @show findmax(vFuncsNew.diffs)[1], findmax(vFuncsNew.diffs)[2]
        row, col, hgt = findmax(vFuncsNew.diffs)[2][1], findmax(vFuncsNew.diffs)[2][2], findmax(vFuncsNew.diffs)[2][3]
        # @show vFuncsNew.VF[row, col, hgt]
        # @show vFuncs.VF[row, col, hgt]
        diffs = norm(vFuncsNew.diffs, Inf); # infinity norm of the difference of VF, already a scalar

        # 2. update vFuncsNew objects to vFuncs
      # copyto!(vFuncs.VF, vFuncsNew.VF); # update EVF
       #  copyto!(vFuncs.qBond, vFuncsNew.qBond) # update qBond
       #  copyto!(vFuncs.X, vFuncsNew.X) # update qBond

        # 2-2. update in a weighted fashion to prevent ossilation
        @. vFuncs.VF = (1 - 0.3) * vFuncs.VF + 0.3 * vFuncsNew.VF
        @. vFuncs.qBond = (1 - 0.3) * vFuncs.qBond + 0.3 * vFuncsNew.qBond
        @. vFuncs.X = (1 - 0.3) * vFuncs.X + 0.3 * vFuncsNew.X
        return diffs
end

# debt pricing in counterfactural regime
function qBond_condiState(params::Params{T,S}, Rl::T, il::S ,is::S, ib::S, iDelta::S) where {T<:Real,S<:Integer} 
        
        l, s, b, Delta = params.lGrid[il], params.sGrid[is], params.bGrid[ib], params.deltaGrid[iDelta]
        lambdaStar(l,s,b,Delta) = (Rl*(l+params.g*Delta)+(1+params.Rf)*s-Delta-b)/(Rl*l)
        # capitalRequire(l::T,s::T,b::T,lambda::T)::T = (1-params.α*params.wr)*(1-lambda)*l + s - b - (1-params.g)*params.deltaGrid[iDelta] # bank fails if negative

        function return_temp_underFailure(l,s,b,Delta,Lambda)::T
            b == 0 && return zero(T)
            term = Rl*((1-Lambda)*l+params.g*Delta)+(1+params.Rf)*s 
            num = max(min(b,max(zero(T),params.cF*term-Delta)),term-Delta-params.α*params.wr*Rl*(1-Lambda)*l)
            return num/b
        end

       lambdaStar_val = lambdaStar(l,s,b,Delta) # calculate lambdaStar for the given state (l,s,b,Delta)
       # @show return_temp_underFailure.(l,s,b,Delta,params.lambdaGrid)
       return_temp_underFailure_val = ntuple(i -> return_temp_underFailure(l,s,b,Delta,params.lambdaGrid[i]), 3) # calculate return_temp_underFailure for the given state (l,s,b,Delta)
        qBond_temp = zeros(T, 3)

        if lambdaStar_val < params.λL # Fail all the time  
          #  println("1")
            qBond_temp .= return_temp_underFailure_val
        elseif (params.λL <= lambdaStar_val) && (lambdaStar_val < params.λM)
           # println("2")
            qBond_temp .= (one(T), return_temp_underFailure_val[2], return_temp_underFailure_val[3])
        elseif (params.λM <= lambdaStar_val) && (lambdaStar_val< params.λH)
           # println("3")
            qBond_temp .= (one(T), one(T), return_temp_underFailure_val[3])
        else # No failure at all
          #  println("4")
            qBond_temp .= one(T)
        end

       # @show qBond_temp
       # println("beta")
       # @show params.β
        #println("F")
        #@show params.F


        return params.β * (params.F * qBond_temp)
end

function qBond_specialRegime(params::Params{T,S}, vFuncsNew::VFuncsNew{T,S}, Rl::T) where {T<:Real,S<:Integer}
        for il in eachindex(params.lGrid)
            for is in eachindex(params.sGrid)
                for ib in eachindex(params.bGrid)
                    for iDelta in eachindex(params.deltaGrid)
                        vFuncsNew.qBond[il,is,ib,iDelta,:] .= qBond_condiState(params,Rl,il,is,ib,iDelta)
                    end
                end
            end
        end
end

function Update_vFuncs_Diffs_stationary(vFuncs::VFuncs{T,S}, vFuncsNew::VFuncsNew{T,S}, params::Params{T,S}) where {T<:Real, S<:Integer}

    # 1. calculate the difference of Γ_0 and Γ_1
    A = vFuncsNew.Γ - vFuncs.Γ; # difference between Γ and Γ_new
    diffs = norm(A, Inf); # infinity norm of the difference of Γ

     # 2. update vFuncsNew objects to vFuncs
    # vFuncs.Γ .= vFuncsNew.Γ; # update Γ
    copyto!(vFuncs.Γ, vFuncsNew.Γ) # copy Γ_new to Γ

    return diffs
end

function Update_stationary_dist(params::Params{T,S}, vFuncs::VFuncs{T,S}, vFuncsNew::VFuncsNew{T,S}, Iterobj_is::Matrix{IterObj_i{T,S}}, iDelta::S, iLambda::S, iN::S, iLambdaPrime::S, Rl::T) where {T<:Real,S<:Integer}
    lambda, δ = params.lambdaGrid[iLambdaPrime], params.deltaGrid[iDelta]; # get the state for Lambda Prime
    l, s, b = params.lGrid[iLambdaPrime], params.sGrid[iLambdaPrime], params.bGrid[iLambdaPrime]

    NAV(l::T,s::T,b::T,lambda::T)::T = Rl*( (1-lambda)*l + params.g*δ) + (1+params.Rf)*s - δ - b # net asset value of the bank, lambda here is lambda prime 
    tax(l::T,s::T,b::T,lambda::T)::T = params.τC * max(0, (Rl-1)*((1-lambda)*l +params.g*δ) + params.Rf*s - params.Rf *b - params.Rf * 0.8 * δ) # tax on the bank's asset value
    n_failure(l::T, lambda::T)::T = params.α * params.wr * Rl * (1-lambda)*l # next period asset conditional on bank failure 
    n_success(l::T,s::T,b::T,lambda::T)::T = NAV(l,s,b,lambda) - tax(l,s,b,lambda) # next period asset conditional on bank success 

    # 1. determining failure or success and consequential nPrime
    if Iterobj_is[iDelta, iLambda].failure[iN][iLambda] == 1 # if fail in (iDelta, iLambda, in, iLambda)
        nPrime = n_failure(l, lambda)
    elseif Iterobj_is[iDelta, iLambda].failure[iN][iLambda] == 0 # if success in (iDelta, iLambda, in, iLambda)
        nPrime = n_success(l, s, b, lambda)
    else
        error("unhandled failure code: ")
    end

    # 2. assigning mass weight 
    if searchsortedlast(params.nGrid, nPrime) == 0 # if n is smaller than the minimum nGrid
        vFuncsNew.Γ[:, iLambdaPrime, 1] .+= params.H[iDelta, :] .* vFuncs.Γ[iDelta, iLambda, iN] .* params.F[iLambda, iLambdaPrime]
    elseif searchsortedlast(params.nGrid, nPrime) == length(params.nGrid) # if n is larger than the maximum nGrid
        vFuncsNew.Γ[:, iLambda, length(params.nGrid)] .+= params.H[iDelta, :] .* vFuncs.Γ[iDelta, iLambda, iN] .* params.F[iLambda, iLambdaPrime]
    else # set in between
        iNn, iNnPrime = searchsortedlast(params.nGrid, nPrime), searchsortedfirst(params.nGrid, nPrime)
        vFuncsNew.Γ[:,iLambdaPrime,iNnPrime] .+= params.H[iDelta, :] .* ((nPrime - params.nGrid[iNn]) / (params.nGrid[iNnPrime] - params.nGrid[iNn])) * vFuncs.Γ[iDelta, iLambda, iN] .* params.F[iLambda, iLambdaPrime]
        vFuncsNew.Γ[:,iLambdaPrime,iNn] .+= params.H[iDelta, :] .* ((params.nGrid[iNnPrime] - nPrime) / (params.nGrid[iNnPrime] - params.nGrid[iNn])) * vFuncs.Γ[iDelta, iLambda, iN] .* params.F[iLambda, iLambdaPrime]
    end

end

function stationary_distribution(params::Params{T,S}, Rl::T, vFuncs::VFuncs{T,S},vFuncsNew::VFuncsNew{T,S},Iterobj_is::Matrix{IterObj_i{T,S}}, maxiter::S, tol::T) where {T<:Real,S<:Integer}

    iter, maxdiffs, diffs = 1, one(T), zeros(T);

    while iter <= maxiter && maxdiffs > tol

       # println("stationary_distribution: iter = ", iter, ", maxdiffs = ", maxdiffs); # report the iteration progress
        Threads.@threads for iDelta in 1:length(params.deltaGrid) # 3: number of states for Delta
            for iLambda in 1:length(params.lambdaGrid)  # 3: number of states for Lambda
                for iN in 1:length(params.nGrid)
                    for iLambdaPrime in 1:length(params.lambdaGrid)
                        Update_stationary_dist(params, vFuncs, vFuncsNew, Iterobj_is, iDelta, iLambda, iN, iLambdaPrime, Rl)
                    end
                end
            end
        end
        
        vFuncsNew.Γ ./= sum(vFuncsNew.Γ) # normalize the distribution
        diffs = Update_vFuncs_Diffs_stationary(vFuncs, vFuncsNew, params); # after updating the mass, calculate the difference
        maxdiffs = maximum(diffs);

        if mod(iter, 200) == 0
            println("stationary_distribution: iter=", iter, ", maxdiffs=", maxdiffs); # report the iteration progress 
        end

        # 3-3. if the difference is not small enough, do the iteration again
        iter += 1;
    end

    return (vFuncs = vFuncs, vFuncsNew = vFuncsNew, Iterobj_is = Iterobj_is);  
end

# aggregate loan supply, including government guaranteed loan 
function aggre_loan_supply(params::Params{T,S}, vFuncs::VFuncs{T,S}, Iterobj_is::Matrix{IterObj_i{T,S}}) where {T<:Real,S<:Integer}
        l_given_state = [Iterobj_is[a,b].solution[c][1] for a in 1:3, b in 1:3, c in 1:length(params.nGrid)] # 3x3xn_npts matrix 
        public_given_state = [params.g * params.deltaGrid[a] for a in 1:3] # 3x1 matrix, government guaranteed loan

        loan_supply = (l_given_state .+ public_given_state) .* vFuncs.Γ; # 3x3xn_npts matrix, total loan supply
        loan_supply = sum(loan_supply) * params.M # aggregate loan supply, 1xn_npts matrix
    return loan_supply
end

function aggre_loan_supply2(params::Params{T,S}, vFuncs::VFuncs{T,S}, Iterobj_is::Matrix{IterObj_i{T,S}}) where {T<:Real,S<:Integer}
        l_given_state = [Iterobj_is[a,b].solution[c][1] for a in 1:3, b in 1:3, c in 1:length(params.nGrid)] # 3x3xn_npts matrix 
        public_given_state = [params.g * params.deltaGrid[a] for a in 1:3] # 3x1 matrix, government guaranteed loan

        loan_supply = (l_given_state .+ public_given_state) .* vFuncs.Γ; # 3x3xn_npts matrix, total loan supply
        loan_supply = sum(loan_supply) * params.M # aggregate loan supply, 1xn_npts matrix

        loan_supply_private = sum(l_given_state.* vFuncs.Γ) * params.M
        loan_supply_public = sum(public_given_state.* vFuncs.Γ) * params.M

    return (loan_supply =  loan_supply, loan_supply_private = loan_supply_private, loan_supply_public = loan_supply_public)
end

################################################################################################################################################
###### solving for equilibrium given parameters and regime ########

function solve_model_given_r(Rl::T; Params::Params{T,S}, Regime::F) where {T<:Real,S<:Integer,F<:Bool} 

    println("SUB BLOCK BEGAINS -- VFI calculation begins:   Regime is ", Regime) 
    eq = VFI(Params, Rl, Regime, 1000, 1e-2); # run the VFI algorithm with the given parameters and regime

    println("SUB BLOCK BEGAINS -- staionary distribution calculation begins") 
    eqq = stationary_distribution(Params, Rl, eq.vFuncs, eq.vFuncsNew, eq.Iterobj_is, 1000, 1e-4); # run the stationary distribution algorithm with the given parameters and regime
    excess_loan_supply = aggre_loan_supply(Params, eqq.vFuncs, eqq.Iterobj_is) - Rl^(Params.ϵ)*Params.E
    # excess_loan_supply = aggre_loan_supply(Params, eqq.vFuncs, eqq.Iterobj_is) - Rl^(Params.ϵ)
    return excess_loan_supply

    println("SUB BLOCK ENDED") 
end

function solve_model_given_r2(Rl::T; Params::Params{T,S}, Regime::F) where {T<:Real,S<:Integer,F<:Bool} 

    println("SUB BLOCK BEGAINS -- VFI calculation begins:   Regime is ", Regime) 
    @time eq = VFI(Params, Rl, Regime, 1000, 1e-3); # run the VFI algorithm with the given parameters and regime

    println("SUB BLOCK BEGAINS -- staionary distribution calculation begins") 
    @time eqq = stationary_distribution(Params, Rl, eq.vFuncs, eq.vFuncsNew, eq.Iterobj_is, 1000, 1e-4); # run the stationary distribution algorithm with the given parameters and regime
    loan_supply = aggre_loan_supply2(Params, eqq.vFuncs, eqq.Iterobj_is) 
    loan_demand = Rl^(Params.ϵ)*Params.E
    excess_loan_supply = loan_supply.loan_supply - loan_demand
    # excess_loan_supply = aggre_loan_supply(Params, eqq.vFuncs, eqq.Iterobj_is) - Rl^(Params.ϵ)
    return (eq = eq, eqq = eqq, loan_supply = loan_supply, loan_demand = loan_demand, excess_loan_supply = excess_loan_supply)

    println("SUB BLOCK ENDED") 
end

# (FINAL) output: equilibrium loan rate
function solve_model(params::Params{T,S}, regime::F, a::T, b::T) where {T<:Real,S<:Integer,F<:Bool}
    
    solve_model_given_r_single = Rl -> solve_model_given_r(Rl; Params = params, Regime = regime)
    println("MAIN BLOCK BEGAINS:   Regime is ", regime) 

    # @show solve_model_given_r_single(a)
    # @show solve_model_given_r_single(b)
    # @assert solve_model_given_r_single(a) * solve_model_given_r_single(b) < 0 "No root in interval"
    
    # println("end points condition satisfied, proceeding to find root") 
    # @show A = solve_model_given_r_single(a)
    # @show B = solve_model_given_r_single(b)
    # return (A = A, B = B)
    # Rl_star = find_zero(solve_model_given_r_single, (a,b), Bisection(); tol=1e-2, maxevals=1000);
    Rl_safe = clamp(a, 1.0+eps(), 2.0-eps())
    Rl_star = find_zero(solve_model_given_r_single, Rl_safe, method=Roots.Secant(); tol = 1e-2, maxevals = 100);
    # println("Rl_star", Rl_star) 

   # fa = solve_model_given_r_single(a)
   # fb = solve_model_given_r_single(b)
   # @assert fa * fb < 0 "The function must have opposite signs at the endpoints a and b."

    # obtain the equilibrium under equilibrium loan rate Rl_star
    eq =  VFI(params, Rl_star, regime, 1000, 1e-4);
    eqq = stationary_distribution(params, Rl_star, eq.vFuncs, eq.vFuncsNew, eq.Iterobj_is, 1000, 1e-4);

    # update-in-place to reuse the interim solution
    # VFI!(params, Rl_star, regime, 1000, 1e-4);
    return (Rl_star = Rl_star, eq = eq, eqq = eqq)
    
end

################################################################################################################################################
## code for policy function and the associated variables under the optimal choice (for example, bond price, tax, dividend, etc)
struct PolicyFuncs{T<:Real,S<:Integer}
    lPolicy::Array{T,3} # loan policy function (delta, lambda, n)
    sPolicy::Array{T,3} # savings policy function (delta, lambda, n)
    bPolicy::Array{T,3} # bond policy function (delta, lambda, n)
    failure::Array{Bool,4} # failure or not (delta, lambda, n, lambdaPrime)
    nPrimePolicy::Array{T,4} # next period net asset value (delta, lambda, n, lambdaPrime)
    NAV::Array{T,4} # net asset value (delta, lambda, n, lambdaPrime)
    divPolicy::Array{T,3} # dividend policy function (delta, lambda, n) under the optimal choice of (l,s,b)
    qBondPolicy::Array{T,3} # bond price policy function (delta, lambda, n)
    taxPolicy::Array{T,3} # tax policy function (delta, lambda, n)
end

function Initiate_PolicyFuncs(params::Params{T,S}) where {T<:Real,S<:Integer}

    lPolicy = zeros(T, length(params.deltaGrid), length(params.lambdaGrid), length(params.nGrid)); # loan policy function (delta, lambda, n)
    sPolicy = zeros(T, length(params.deltaGrid), length(params.lambdaGrid), length(params.nGrid)); # savings policy function
    bPolicy = zeros(T, length(params.deltaGrid), length(params.lambdaGrid), length(params.nGrid)); # bond policy function
    failure = falses(length(params.deltaGrid), length(params.lambdaGrid), length(params.nGrid), length(params.lambdaGrid)); # failure or not (delta, lambda, n)
    nPrimePolicy = zeros(T, length(params.deltaGrid), length(params.lambdaGrid), length(params.nGrid), length(params.lambdaGrid)); # next period net asset value (delta, lambda, n, lambdaPrime)
    NAV = zeros(T, length(params.deltaGrid), length(params.lambdaGrid), length(params.nGrid), length(params.lambdaGrid)); # net asset value (delta, lambda, n, lambdaPrime)
    divPolicy = zeros(T, length(params.deltaGrid), length(params.lambdaGrid), length(params.nGrid)); # dividend policy function (delta, lambda, n) under the optimal choice of (l,s,b)
    qBondPolicy = zeros(T, length(params.deltaGrid), length(params.lambdaGrid), length(params.nGrid));
    taxPolicy = zeros(T, length(params.deltaGrid), length(params.lambdaGrid), length(params.nGrid));

    policyy = PolicyFuncs{T,S}(lPolicy, sPolicy, bPolicy, failure, nPrimePolicy, NAV, divPolicy, qBondPolicy, taxPolicy);
    return policyy;
end

function Get_PolicyFuncs(params::Params{T,S},Iterobj_is::Matrix{IterObj_i{T,S}},vFuncs::VFuncs{T,S},Rl::T) where {T<:Real,S<:Integer}

    NAV(l::T,s::T,b::T,lambda::T,Rl::T,δ::T)::T = Rl*( (1-lambda)*l + params.g*δ) + (1+params.Rf)*s - δ - b # 이익잉여금, lambda here is lambda prime 
    tax(l::T,s::T,b::T,lambda::T,Rl::T,δ::T)::T = params.τC * max(0, (Rl-1)*((1-lambda)*l +params.g*δ) + params.Rf*s - params.Rf *b - params.Rf* 0.8 *δ  ) # tax on the bank's asset value
    n_failure(l::T, lambda::T,Rl::T)::T = params.α * params.wr * Rl * (1-lambda)*l # next period asset conditional on bank failure 
    n_success(l::T,s::T,b::T,lambda::T,Rl::T,δ::T)::T = NAV(l,s,b,lambda,Rl,δ) - tax(l,s,b,lambda,Rl,δ) # next period asset conditional on bank success 
    # 자기자본규제 미충족시 은행이 실패한다고 가정할 때
    capitalRequire(l::T,s::T,b::T,lambda::T,δ::T)::T = (1-params.α*params.wr)*(1-lambda)*l + s - b - (1-params.g)*δ # bank fails if negative

    # Interpolate qBond for a given value of choice (l,s,b) and state (delta, lambda)
    function qBond_interpolated(l::T, s::T, b::T, deltaInd::S, lambdaInd::S, params::Params{T,S}, vFuncs::VFuncs{T,S})::T

        lRange = range(first(params.lGrid), last(params.lGrid), length(params.lGrid));
        sRange = range(first(params.sGrid), last(params.sGrid), length(params.sGrid));
        bRange = range(first(params.bGrid), last(params.bGrid), length(params.bGrid));
        

        q = vFuncs.qBond[:, :, :, deltaInd, lambdaInd];
        qmin = minimum(q);
        qmax = maximum(q);

        if qmin == qmax # if qBond is constant, return the constant value
            qnorm = qmin .* ones(size(q));
        else
            qnorm = (q .- qmin) ./ (qmax - qmin);
        end
        
        q_itp_rev = interpolate(qnorm, BSpline(Quadratic(Line(OnGrid()))));
        q_itp_rev = Interpolations.scale(q_itp_rev, lRange, sRange, bRange);
        q_itp_ext_rev = extrapolate(q_itp_rev, Line());
        q_interp_rev(ll, ss, bb) = q_itp_ext_rev(ll, ss, bb);

        q_value = q_interp_rev(l, s, b) * (qmax - qmin) + qmin; # interpolate qBond for a given choice (l,s,b)
        return q_value
    end

    # 1. initiate policy functions
    policy = Initiate_PolicyFuncs(params);

    # 3. calculate the policy functions (l,s,b) and failure or not
    for iDelta in 1:length(params.deltaGrid)
        for iLambda in 1:length(params.lambdaGrid)
            iterobj = Iterobj_is[iDelta, iLambda]; 

            for iN in 1:length(params.nGrid)

                policy.lPolicy[iDelta, iLambda, iN] = iterobj.solution[iN][1]; # loan policy function
                policy.sPolicy[iDelta, iLambda, iN] = iterobj.solution[iN][2]; # savings policy function
                policy.bPolicy[iDelta, iLambda, iN] = iterobj.solution[iN][3]; # bond policy function

                # policy.failure[iDelta, iLambda, iN, :] .= [NAV(iterobj.solution[iN][1], iterobj.solution[iN][2], iterobj.solution[iN][3], params.lambdaGrid[iLambdaPrime], Rl, params.deltaGrid[iDelta]) < 0 for iLambdaPrime in 1:length(params.lambdaGrid)]; # failure or not, denoted as 1 or 0
                policy.failure[iDelta, iLambda, iN, :] .= [capitalRequire(iterobj.solution[iN][1], iterobj.solution[iN][2], iterobj.solution[iN][3], params.lambdaGrid[iLambdaPrime], params.deltaGrid[iDelta]) < 0 for iLambdaPrime in 1:length(params.lambdaGrid)]; # failure or not, denoted as 1 or 0
                policy.nPrimePolicy[iDelta, iLambda, iN, :] .= [policy.failure[iDelta, iLambda, iN, iLambdaPrime] == 1 ? n_failure(iterobj.solution[iN][1], params.lambdaGrid[iLambdaPrime], Rl) : n_success(iterobj.solution[iN][1], iterobj.solution[iN][2], iterobj.solution[iN][3], params.lambdaGrid[iLambdaPrime], Rl,params.deltaGrid[iDelta]) for iLambdaPrime in 1:length(params.lambdaGrid)]
                policy.NAV[iDelta, iLambda, iN, :] .= [NAV(iterobj.solution[iN][1], iterobj.solution[iN][2], iterobj.solution[iN][3], params.lambdaGrid[iLambdaPrime], Rl, params.deltaGrid[iDelta]) for iLambdaPrime in 1:length(params.lambdaGrid)]
                q = qBond_interpolated(iterobj.solution[iN][1], iterobj.solution[iN][2], iterobj.solution[iN][3], iDelta, iLambda, params, vFuncs); # bond price under the optimal choice of (l,s,b) at state (delta, lambda)
                policy.qBondPolicy[iDelta, iLambda, iN] = q; # bond price policy function
                policy.divPolicy[iDelta, iLambda, iN] = params.nGrid[iN] + (1-params.cL)*params.β*params.deltaGrid[iDelta] + q*iterobj.solution[iN][3] - iterobj.solution[iN][1] - params.g*params.deltaGrid[iDelta] - iterobj.solution[iN][2] - params.cM*iterobj.solution[iN][1]^2 - params.cO; # dividend policy function under the optimal choice of (l,s,b)
                policy.taxPolicy[iDelta, iLambda, iN] = tax(iterobj.solution[iN][1], iterobj.solution[iN][2], iterobj.solution[iN][3], params.lambdaGrid[iLambda], Rl, params.deltaGrid[iDelta]); # tax policy function under the optimal choice of (l,s,b)
            end
        end
    end

    return policy;
end

################################################################################################################################################
## code for simulation and calculating moments
struct SimShocks{S<:Integer}
    deltaIndSim::Array{S,3} # N * J * T
    lambdaIndSim::Array{S,3} # N * J * T
    nInitialIndSim::Array{S,2} # initial net asset value (index) for each bank J -- N * J
end

struct SimPaths{T<:Real,S<:Integer}
    deltaIndSim::Array{S,3}
    lambdaIndSim::Array{S,3}
    nInitialIndSim::Array{S,2} # initial net asset value (index) for each bank J
    bigN::S # number of observations (N)
    bigJ::S # number of banks (J)
    bigT::S # number of periods (T, year)
    deltaSim::Array{T,3} # delta states for each bank J
    lambdaSim::Array{T,3} # lambda states for each bank J
    lSim::Array{T,3} # loan 
    sSim::Array{T,3} # safe-asset saving
    bSim::Array{T,3} # bond issuance, bank borrowing 
    failureSim::Array{S,3} # failure or not, denoted as 1 or 0
    nSim::Array{T,3} # retailed earning, 이익잉여금, current period 
    divSim::Array{T,3} # dividend(배당금) OR equity(자본); current period
    nIndSim::Array{S,3} # net asset value (index), current period 
    nPrimeSim::Array{T,3} # net asset value, next period with lambdaPrime
    nPrimeIndSim::Array{S,3} # net asset value (index), next period with lambdaPrime
    assetSim::Array{T,3} # asset (loan before lambda + savings), current period
    assetPrimeSim::Array{T,3} # asset, next period with lambdaPrime
    liabilitySim::Array{T,3} # liability (deposit + bank debt), current period
    Γ::Array{T,3} # distribution of states, current period
    ΓPrime::Array{T,3} # distribution of states, next period with lambdaPrime
    qBondSim::Array{T,3} # bond price under the choice of (l,s,b) at state (delta, lambda), current period 
    gSim::Array{T,3} # government guaranteed loan, current period
    govSpendingSim::Array{T,3} # government spending for bank bailout, next period with lambdaPrime
    aggGovSpendingSim::Array{T,2} # aggregate government spending for bank bailout, next period with lambdaPrime
    depositLoanRatioSim::Array{T,3} # deposit to loan ratio, current period
    leverageSim::Array{T,3} # leverage ratio, current period
    debtToLiability::Array{T,3} # for moment1
    capitalToDeposit::Array{T,3} # for moment2
    loanToAsset::Array{T,3} # for moment3
    loanToDeposit::Array{T,3}
end

# storing data moments
struct Params_cal{T<:Real,S<:Integer}
    bigT::S
    bigN::S
    bigJ::S
    trim::S
    a::T
    b::T
    debt_to_liability::T
    capital_to_deposit::T
    loan_to_asset::T
    loan_rate::T
    cM::T
    cO::T
    cL::T
    ϵ::T
    E::T
    dBar::T
    σ::T
end

function Initiate_Params_cal(bigT::S,bigN::S,bigJ::S,trim::S,a::T,b::T,debt_to_liability::T,capital_to_deposit::T,loan_to_asset::T,loan_rate::T, cM::T,cO::T,cL::T,ϵ::T,E::T,dBar::T,σ::T) where {T<:Real,S<:Integer}

    pam = Params_cal{T,S}(bigT,bigN,bigJ,trim,a,b,debt_to_liability,capital_to_deposit,loan_to_asset,loan_rate,cM,cO,cL,ϵ,E,dBar,σ);
    return pam
end

# draw bank types from Γ and simulate (delta, lambda) from Markov Chain 
function makeShockSequence(params::Params{T,S}, vFuncs::VFuncs, bigT::S, bigN::S, bigJ::S) where {T<:Real,S<:Integer}

    deltaIndSim = zeros(Int, bigT, bigJ, bigN); # delta states for each bank J
    lambdaIndSim = zeros(Int, bigT, bigJ, bigN); # lambda states for each bank J
    nInitialIndSim = zeros(Int, bigJ, bigN); # initial net asset value (index) for each bank J

    # 1. set up drawing indices from the stationary distribution Γ
    Γ = vFuncs.Γ; # stationary distribution of states, Γ[iDelta, iLambda, iN] = P(delta = iDelta, lambda = iLambda, n = iN)
    Γ ./= sum(Γ); 
    items = collect(CartesianIndices(Γ))
    w = Weights(vec(Γ))    
    idx = sample(items, w, bigN * bigJ) # draw indices from the stationary distribution Γ, size = (bigN * bigJ, 3)

    # 2. draw initial state for each bank J for each observation N
    for smallN = 1:bigN
        for smallJ = 1:bigJ # for bank J
            number = (smallN - 1) * bigJ + smallJ; # index for the bank J in the observation N
            deltaIndSim[1, smallJ, smallN] = idx[number][1]; 
            lambdaIndSim[1, smallJ, smallN] = idx[number][2]; 
            nInitialIndSim[smallJ, smallN] = idx[number][3]; 
        end
    end

    # 3. for each period T, draw (delta, lambda) from the markov process 
    H_Markov = [Categorical(params.H[i, :]) for i in 1:size(H,1)]
    F_Markov = [Categorical(params.F[i, :]) for i in 1:size(F,1)]

    # @show bigT
    # @show bigN
    for smallN = 1:bigN
        for smallJ = 1:bigJ # for bank J
            for smallT = 2:bigT # for each period T
                deltaIndSim[smallT, smallJ, smallN] = rand(H_Markov[Int(deltaIndSim[smallT-1, smallJ, smallN])]) # draw delta from the Markov chain
                lambdaIndSim[smallT, smallJ, smallN] = rand(F_Markov[Int(lambdaIndSim[smallT-1, smallJ, smallN])]) # draw lambda from the Markov chain                
            end
        end
    end

   return SimShocks{S}(deltaIndSim, lambdaIndSim, nInitialIndSim); # store the simulated shocks in SimShocks struct
end

function Initiate_Paths(shocks::SimShocks{S}, params::Params{T,S}) where {T<:Real,S<:Integer}

    deltaIndSim = shocks.deltaIndSim;
    lambdaIndSim = shocks.lambdaIndSim;
    nInitialIndSim = shocks.nInitialIndSim;

    bigT = size(deltaIndSim, 1); # number of periods (T, year)
    bigJ = size(deltaIndSim, 2); # number of observations (N)
    bigN = size(deltaIndSim, 3); # number of banks

    # below are variables needed to be filled
    deltaSim = zeros(T, bigT, bigJ, bigN);
    lambdaSim = zeros(T, bigT, bigJ, bigN);
    lSim = zeros(T, bigT, bigJ, bigN); # loan policy function
    sSim = zeros(T, bigT, bigJ, bigN); # savings policy function
    bSim = zeros(T, bigT, bigJ, bigN); # bond policy function
    failureSim = zeros(S, bigT, bigJ, bigN); # failure or not, denoted as 1 or 0
    nSim = zeros(T, bigT, bigJ, bigN); # retained earining (잉여금; 기업의 이익 중에서 배당금으로 배분하지 않고 재투자를 위해 남겨둔 금액); current period
    divSim = zeros(T, bigT, bigJ, bigN); # dividend (배당금); current period
    nIndSim = zeros(S, bigT, bigJ, bigN); # net asset value (index); current period
    nPrimeSim = zeros(T, bigT, bigJ, bigN); # retained earining (잉여금; 기업의 이익 중에서 배당금으로 배분하지 않고 재투자를 위해 남겨둔 금액); next period with lambdaPrime
    nPrimeIndSim = zeros(S, bigT, bigJ, bigN); # net asset value (index); next period with lambdaPrime
    assetSim = zeros(T, bigT, bigJ, bigN); # total asset; current period
    assetPrimeSim = zeros(T, bigT, bigJ, bigN); # total asset; next period with lambdaPrime
    liabilitySim = zeros(T, bigT, bigJ, bigN); # total liability; current period
    Γ = zeros(3, 3, length(params.nGrid)); # stationary distribution
    ΓPrime = zeros(3, 3, length(params.nGrid)); # next period stationary distribution
    qBondSim = zeros(T, bigT, bigJ, bigN); # bond price under the choice of (l,s,b) at state (delta, lambda), current period
    gSim = zeros(T, bigT, bigJ, bigN); # government guaranteed loan, current period
    govSpendingSim = zeros(T, bigT, bigJ, bigN); # government spending for bank bailout, next period with lambdaPrime
    aggGovSpendingSim = zeros(T, bigT, bigN); # aggregate government spending for bank bailout, next period with lambdaPrime
    depositLoanRatioSim = zeros(T, bigT, bigJ, bigN); # deposit to loan ratio, current period
    leverageSim = zeros(T, bigT, bigJ, bigN); # leverage ratio, current period
    debtToLiability = zeros(T, bigT, bigJ, bigN);
    capitalToDeposit = zeros(T, bigT, bigJ, bigN);
    loanToAsset = zeros(T, bigT, bigJ, bigN);
    loanToDeposit = zeros(T, bigT, bigJ, bigN);

    return SimPaths(deltaIndSim, lambdaIndSim, nInitialIndSim, bigN, bigJ, bigT, deltaSim, lambdaSim, lSim, sSim, bSim, failureSim, nSim, divSim, nIndSim, nPrimeSim, nPrimeIndSim, assetSim, assetPrimeSim, liabilitySim, Γ, ΓPrime, qBondSim, gSim, govSpendingSim, aggGovSpendingSim, depositLoanRatioSim, leverageSim, debtToLiability, capitalToDeposit, loanToAsset, loanToDeposit);
end

# after initiating paths, simulate the eq and calculate moments of interest
function Simulate_paths(paths::SimPaths{T,S}, policy::PolicyFuncs{T,S}, vFuncs::VFuncs{T,S}, params::Params{T,S}, Rl::T, trim::S,regime::F) where {T<:Real,S<:Integer,F<:Bool}

    bigT = size(paths.deltaIndSim, 1); # number of observations (N)
    bigJ = size(paths.deltaIndSim, 2); # number of banks
    bigN = size(paths.deltaIndSim, 3); # number of periods (T, year)

    # 1. for bank smallJ calculate the series of interest at time (smallT, smallN) or given state at (smallT, smallN)
    for smallN = 1:bigN
        for smallJ = 1:bigJ
            for smallT = 1:bigT-1

                calculate_series(paths, policy, vFuncs, params, Rl, smallT, smallN, smallJ,regime)
            end
        end
    end
    #return paths

    # 2. after obtaining simulation paths, trim some first parts and calculate moments
    moments = calculate_moments(paths, policy, vFuncs, params, Rl, trim);

    # 3. return key moments of interest
    return moments 
end

# for bank smallJ, calculate the variable of interest at time (smallT, smallN); given state (delta, lambda, n), calculate the series of interest
## for calibration, don't need to calculate Gamma 
function calculate_series(paths::SimPaths{T,S}, policy::PolicyFuncs{T,S}, vFuncs::VFuncs{T,S}, params::Params{T,S}, Rl::T, smallT::S, smallN::S, smallJ::S,regime::F) where {T<:Real,S<:Integer,F<:Bool}

    NAV(l::T,s::T,b::T,lambda::T,δ::T)::T = Rl*( (1-lambda)*l + params.g*δ) + (1+params.Rf)*s - δ - b; # 순자산, net asset value 
    totalAsset(l::T,s::T,lambda::T,δ::T)::T = Rl * ( (1-lambda)*l + params.g*δ) + (1+params.Rf)*s; # total asset of the bank, 총자산
    tax(l::T,s::T,b::T,lambda::T,δ::T)::T = params.τC * max(0, (Rl-1)*((1-lambda)*l +params.g*δ) + params.Rf*s - params.Rf * b - params.Rf * 0.8 * δ ); # tax on the bank's asset value
    n_success(l::T,s::T,b::T,lambda::T,δ::T)::T = NAV(l,s,b,lambda,δ) - tax(l,s,b,lambda,δ); # next period asset conditional on bank success 
    n_failure(l::T, lambda::T)::T = params.α * params.wr * Rl * (1-lambda)*l; # next period asset conditional on bank failure 
    govGuarantedSpending(lambda::T, δ::T)::T = Rl * (lambda + params.ξ) * params.g * δ; 
    capitalRequire(l::T,s::T,b::T,lambda::T,δ::T)::T = (1-params.α*params.wr)*(1-lambda)*l + s - b - (1-params.g)*δ; # bank fails if negative

    ### government spending for regime == true 
    # bHat_true(delta::T,lambda::T,l::T,s::T,b::T)::T =  # Rl*[(1-lambda)*l + params.g*delta] + (1+params.Rf)*s - delta; # when the bank survives, the amount of liability that needs to be bailed out by the government
    # vLenderB_true(delta::T,lambda::T,l::T,s::T,b::T)::T = 
    # theta_true(delta::T,lambda::T,l::T,s::T,b::T)::T = maximum(0, vLenderB_true(delta::T,lambda::T,l::T,s::T,b::T) - bHat_true(delta::T,lambda::T,l::T,s::T,b::T))
    # bailoutCost_true(delta::T,lambda::T,l::T,s::T,b::T)::T = theta_true(delta::T,lambda::T,l::T,s::T,b::T) -
    # government spending for regime == false
    # theta_false(delta::T,lambda::T,l::T,s::T,b::T)::T = delta + b - (1-params.α*params.wr)*Rl*[(1-lambda)*l + params.g*delta] - (1+params.Rf)*s; # when the bank fails, the amount of liability that needs to be bailed out by the government
    # bailoutCost_false(delta::T,lambda::T,l::T,s::T,b::T)::T = theta_false(delta::T,lambda::T,l::T,s::T,b::T) - 

    # Interpolate qBond for a given value of choice (l,s,b) and state (delta, lambda)
    function qBond_interpolated(l::T, s::T, b::T, deltaInd::S, lambdaInd::S, params::Params{T,S}, vFuncs::VFuncs{T,S})::T

        lRange = range(first(params.lGrid), last(params.lGrid), length(params.lGrid));
        sRange = range(first(params.sGrid), last(params.sGrid), length(params.sGrid));
        bRange = range(first(params.bGrid), last(params.bGrid), length(params.bGrid));
        

        q = vFuncs.qBond[:, :, :, deltaInd, lambdaInd];
        qmin = minimum(q);
        qmax = maximum(q);

        if qmin == qmax # if qBond is constant, return the constant value
            qnorm = qmin .* ones(size(q));
        else
            qnorm = (q .- qmin) ./ (qmax - qmin);
        end
        
        q_itp_rev = interpolate(qnorm, BSpline(Quadratic(Line(OnGrid()))));
        q_itp_rev = Interpolations.scale(q_itp_rev, lRange, sRange, bRange);
        q_itp_ext_rev = extrapolate(q_itp_rev, Line());
        q_interp_rev(ll, ss, bb) = q_itp_ext_rev(ll, ss, bb);

        q_value = q_interp_rev(l, s, b) * (qmax - qmin) + qmin; # interpolate qBond for a given choice (l,s,b)
        return q_value
    end

    # Interpolate (l,s,b) for a given n value and a given state (delta, lambda)
    function lsb_interpolated(params::Params{T,S}, policy::PolicyFuncs{T,S}, deltaInd::S, lambdaInd::S, n::T)::NTuple{3,T}
        # Interpolate l, s, b for a given n value and a given state (delta, lambda)

        nRange = range(first(params.nGrid), last(params.nGrid), length(params.nGrid));

        l_array = policy.lPolicy[deltaInd, lambdaInd, :]; # loan policy function for a given state (delta, lambda)
        # @show l_array # check the loan policy function for a given state (delta, lambda)
        l_array_min, l_array_max = minimum(l_array), maximum(l_array);
        if l_array_min == l_array_max # if lPolicy is constant, return the constant value
           #  println("lPolicy is constant, l_array_min = ", l_array_min)
            l_array_norm = l_array_min .* ones(length(l_array)); # normalize the loan policy function
        else
           #  println("lPolicy is not constant, l_array_min = ", l_array_min)
            l_array_norm = (l_array .- l_array_min) ./ (l_array_max - l_array_min);
        end
       #  @show l_array_min, l_array_max, l_array_norm # check the minimum, maximum, and normalized loan policy function
        l_itp = interpolate(l_array_norm, BSpline(Quadratic(Line(OnGrid()))));
        l_itp = Interpolations.scale(l_itp, nRange);    
        l_itp_ext = extrapolate(l_itp, Line());
        l_interp(n_val) = l_itp_ext(n_val);
        l_value = l_interp(n) * (l_array_max - l_array_min) + l_array_min; # interpolate loan policy function for a given n value
        


        s_array = policy.sPolicy[deltaInd, lambdaInd, :]; # savings policy function for a given state (delta, lambda)
        s_array_min, s_array_max = minimum(s_array), maximum(s_array);
        if s_array_min == s_array_max # if sPolicy is constant, return the constant value
            s_array_norm = s_array_min .* ones(length(s_array))
        else
            s_array_norm = (s_array .- s_array_min) ./ (s_array_max - s_array_min);
        end
        s_itp = interpolate(s_array_norm, BSpline(Quadratic(Line(OnGrid()))));
        s_itp = Interpolations.scale(s_itp, nRange);
        s_itp_ext = extrapolate(s_itp, Line());
        s_interp(n_val) = s_itp_ext(n_val);
        s_value = s_interp(n) * (s_array_max - s_array_min) + s_array_min; # interpolate savings policy function for a given n value

        b_array = policy.bPolicy[deltaInd, lambdaInd, :]; # bond policy function for a given state (delta, lambda)
        b_array_min, b_array_max = minimum(b_array), maximum(b_array);
        if b_array_min == b_array_max # if bPolicy is constant, return the constant value
            b_array_norm = b_array_min .* ones(length(b_array))
        else
            b_array_norm = (b_array .- b_array_min) ./ (b_array_max - b_array_min);
        end
        b_itp = interpolate(b_array_norm, BSpline(Quadratic(Line(OnGrid()))));
        b_itp = Interpolations.scale(b_itp, nRange);
        b_itp_ext = extrapolate(b_itp, Line());
        b_interp(n_val) = b_itp_ext(n_val);
        b_value = b_interp(n) * (b_array_max - b_array_min) + b_array_min; # interpolate bond policy function for a given n value

        return (l_value,s_value,b_value)
    end

    # 1. at (smallT, smallN, smallJ), calculate (delta, lambda, lambdaPrime) from shocks
    delta = params.deltaGrid[paths.deltaIndSim[smallT, smallJ, smallN]]; # delta state for bank J at time T
    lambda = params.lambdaGrid[paths.lambdaIndSim[smallT, smallJ, smallN]]; # lambda state for bank J at time T
    lambdaPrime = params.lambdaGrid[paths.lambdaIndSim[smallT+1, smallJ, smallN]]; # lambda state for bank J at time T+1

    paths.deltaSim[smallT, smallJ, smallN] = delta;
    paths.lambdaSim[smallT, smallJ, smallN] = lambda;

    # 2. at (smallT, smallN, smallJ), calculate choice variables (l,s,b) 
    # and the consequential market variables: dividend (div), asset, failure, next period asset (nPrime)
    if smallT == 1 # initial period 
        paths.nSim[smallT, smallJ, smallN] = params.nGrid[paths.nInitialIndSim[smallJ, smallN]]; ##### start with first grid n point
        paths.nIndSim[smallT, smallJ, smallN] = paths.nInitialIndSim[smallJ, smallN]; # index of nSim at (smallT, smallN, smallJ)
        paths.lSim[smallT, smallJ, smallN] = policy.lPolicy[paths.deltaIndSim[smallT, smallJ, smallN], paths.lambdaIndSim[smallT, smallJ, smallN], paths.nInitialIndSim[smallJ, smallN]]; # loan policy function
        paths.sSim[smallT, smallJ, smallN] = policy.sPolicy[paths.deltaIndSim[smallT, smallJ, smallN], paths.lambdaIndSim[smallT, smallJ, smallN], paths.nInitialIndSim[smallJ, smallN]]; # savings policy function
        paths.bSim[smallT, smallJ, smallN] = policy.bPolicy[paths.deltaIndSim[smallT, smallJ, smallN], paths.lambdaIndSim[smallT, smallJ, smallN], paths.nInitialIndSim[smallJ, smallN]]; # bond policy function
        paths.divSim[smallT, smallJ, smallN] = params.nGrid[paths.nInitialIndSim[smallJ, smallN]] + (1-params.cL) * params.β * delta + qBond_interpolated(paths.lSim[smallT, smallJ, smallN], paths.sSim[smallT, smallJ, smallN], paths.bSim[smallT, smallJ, smallN], paths.deltaIndSim[smallT, smallJ, smallN], paths.lambdaIndSim[smallT, smallJ, smallN], params, vFuncs) * paths.bSim[smallT, smallJ, smallN] - paths.lSim[smallT, smallJ, smallN] - params.g * delta - paths.sSim[smallT, smallJ, smallN] - params.cM * paths.lSim[smallT, smallJ, smallN]^2 - params.cO; # dividend (배당금); current period
        paths.assetSim[smallT, smallJ, smallN] = paths.lSim[smallT, smallJ, smallN] + paths.sSim[smallT, smallJ, smallN]; # total asset before lambdaPrime
        paths.assetPrimeSim[smallT, smallJ, smallN] = totalAsset(paths.lSim[smallT, smallJ, smallN], paths.sSim[smallT, smallJ, smallN], lambdaPrime, delta);  # total asset given lambdaPrime
        paths.failureSim[smallT, smallJ, smallN] = (capitalRequire(paths.lSim[smallT, smallJ, smallN], paths.sSim[smallT, smallJ, smallN],paths.bSim[smallT, smallJ, smallN],lambdaPrime, delta) < 0 ) ? 1 : 0; # failure or not, denoted as 1 or 0
        paths.liabilitySim[smallT, smallJ, smallN] = paths.bSim[smallT, smallJ, smallN] + delta; # total liability; current period
        paths.gSim[smallT, smallJ, smallN] = params.g * delta; # government guaranteed loan, current period
        paths.qBondSim[smallT, smallJ, smallN] = qBond_interpolated(paths.lSim[smallT, smallJ, smallN], paths.sSim[smallT, smallJ, smallN], paths.bSim[smallT, smallJ, smallN], paths.deltaIndSim[smallT, smallJ, smallN], paths.lambdaIndSim[smallT, smallJ, smallN], params, vFuncs); # bond price under the choice of (l,s,b) at state (delta, lambda), current period
        paths.leverageSim[smallT, smallJ, smallN] = paths.liabilitySim[smallT, smallJ, smallN] / paths.assetSim[smallT, smallJ, smallN]; # leverage ratio, current period
        paths.debtToLiability[smallT, smallJ, smallN] = paths.bSim[smallT, smallJ, smallN] / paths.liabilitySim[smallT, smallJ, smallN];
        paths.loanToAsset[smallT, smallJ, smallN] = paths.lSim[smallT, smallJ, smallN] / paths.assetSim[smallT, smallJ, smallN]; # 
        paths.loanToDeposit[smallT, smallJ, smallN] = paths.lSim[smallT, smallJ, smallN] / paths.deltaSim[smallT, smallJ, smallN];


        paths.failureSim[smallT, smallJ, smallN]
        if paths.failureSim[smallT, smallJ, smallN] == 1 # if failure
           paths.nPrimeSim[smallT, smallJ, smallN] = n_failure(paths.lSim[smallT, smallJ, smallN], lambdaPrime); # next period asset conditional on bank failure 
        else # if success
           paths.nPrimeSim[smallT, smallJ, smallN] = n_success(paths.lSim[smallT, smallJ, smallN], paths.sSim[smallT, smallJ, smallN], paths.bSim[smallT, smallJ, smallN], lambdaPrime, delta); # next period asset conditional on bank success 
        end
        paths.capitalToDeposit[smallT, smallJ, smallN] = (1-paths.failureSim[smallT, smallJ, smallN])*paths.nPrimeSim[smallT, smallJ, smallN]/delta;

    else # for smallT > 1,
        paths.nSim[smallT, smallJ, smallN] = paths.nPrimeSim[smallT-1, smallJ, smallN] # 잉여금 as a embodied state 
        paths.lSim[smallT, smallJ, smallN] = lsb_interpolated(params, policy, paths.deltaIndSim[smallT, smallJ, smallN], paths.lambdaIndSim[smallT, smallJ, smallN], paths.nSim[smallT, smallJ, smallN])[1]; # loan policy function
        paths.sSim[smallT, smallJ, smallN] = lsb_interpolated(params, policy, paths.deltaIndSim[smallT, smallJ, smallN], paths.lambdaIndSim[smallT, smallJ, smallN], paths.nSim[smallT, smallJ, smallN])[2]; # savings policy function
        paths.bSim[smallT, smallJ, smallN] = lsb_interpolated(params, policy, paths.deltaIndSim[smallT, smallJ, smallN], paths.lambdaIndSim[smallT, smallJ, smallN], paths.nSim[smallT, smallJ, smallN])[3]; # bond policy function
        paths.divSim[smallT, smallJ, smallN] = paths.nSim[smallT, smallJ, smallN] + (1-params.cL) * params.β * delta + qBond_interpolated(paths.lSim[smallT, smallJ, smallN], paths.sSim[smallT, smallJ, smallN], paths.bSim[smallT, smallJ, smallN], paths.deltaIndSim[smallT, smallJ, smallN], paths.lambdaIndSim[smallT, smallJ, smallN], params, vFuncs) * paths.bSim[smallT, smallJ, smallN] - paths.lSim[smallT, smallJ, smallN] - params.g * delta - paths.sSim[smallT, smallJ, smallN] - params.cM * paths.lSim[smallT, smallJ, smallN]^2 - params.cO; # dividend (배당금); current period
        paths.assetSim[smallT, smallJ, smallN] = paths.lSim[smallT, smallJ, smallN] + paths.sSim[smallT, smallJ, smallN]; # total asset before lambdaPrime
        paths.assetPrimeSim[smallT, smallJ, smallN] = totalAsset(paths.lSim[smallT, smallJ, smallN], paths.sSim[smallT, smallJ, smallN], lambdaPrime, delta);  # total asset given lambdaPrime
        paths.failureSim[smallT, smallJ, smallN] = (capitalRequire(paths.lSim[smallT, smallJ, smallN], paths.sSim[smallT, smallJ, smallN],paths.bSim[smallT, smallJ, smallN],lambdaPrime, delta) < 0 ) ? 1 : 0; # failure or not, denoted as 1 or 0
        paths.liabilitySim[smallT, smallJ, smallN] = paths.bSim[smallT, smallJ, smallN] + delta; # total liability; current period
        paths.gSim[smallT, smallJ, smallN] = params.g * delta; # government guaranteed loan, current period
        paths.qBondSim[smallT, smallJ, smallN] = qBond_interpolated(paths.lSim[smallT, smallJ, smallN], paths.sSim[smallT, smallJ, smallN], paths.bSim[smallT, smallJ, smallN], paths.deltaIndSim[smallT, smallJ, smallN], paths.lambdaIndSim[smallT, smallJ, smallN], params, vFuncs); # bond price under the choice of (l,s,b) at state (delta, lambda), current period
        paths.leverageSim[smallT, smallJ, smallN] = paths.liabilitySim[smallT, smallJ, smallN] / paths.assetSim[smallT, smallJ, smallN]; # leverage ratio, current period
        paths.debtToLiability[smallT, smallJ, smallN] = paths.bSim[smallT, smallJ, smallN] / paths.liabilitySim[smallT, smallJ, smallN];
        paths.loanToAsset[smallT, smallJ, smallN] = paths.lSim[smallT, smallJ, smallN] / paths.assetSim[smallT, smallJ, smallN];
        paths.loanToDeposit[smallT, smallJ, smallN] = paths.lSim[smallT, smallJ, smallN] / paths.deltaSim[smallT, smallJ, smallN];

        # 3. depending on failure or not, 
        if paths.failureSim[smallT, smallJ, smallN] == 1 # if failure
            paths.nPrimeSim[smallT, smallJ, smallN] = n_failure(paths.lSim[smallT, smallJ, smallN], lambdaPrime); # next period asset conditional on bank failure 
        else # if success
            paths.nPrimeSim[smallT, smallJ, smallN] = n_success(paths.lSim[smallT, smallJ, smallN], paths.sSim[smallT, smallJ, smallN], paths.bSim[smallT, smallJ, smallN], lambdaPrime, delta); # next period asset conditional on bank success 
        end
        paths.capitalToDeposit[smallT, smallJ, smallN] = (1-paths.failureSim[smallT, smallJ, smallN])*paths.nPrimeSim[smallT, smallJ, smallN]/delta;
    end

    if regime == "true"
        if paths.failureSim[smallT, smallJ, smallN] == 1
            paths.govSpendingSim[smallT, smallJ, smallN] = govGuarantedSpending(lambdaPrime, delta); # + bailoutCost_true; # government spending for bank bailout, next period with lambdaPrime
        else
            paths.govSpendingSim[smallT, smallJ, smallN] = govGuarantedSpending(lambdaPrime, delta); 
        end
    else # regime == "false"
        if paths.failureSim[smallT, smallJ, smallN] == 1
            paths.govSpendingSim[smallT, smallJ, smallN] = govGuarantedSpending(lambdaPrime, delta); #+ bailoutCost_false; # government spending for bank bailout, next period with lambdaPrime
        else
            paths.govSpendingSim[smallT, smallJ, smallN] = govGuarantedSpending(lambdaPrime, delta); # no government spending for bank bailout
        end
    end

end


# after obtaining simulation paths, trim some first parts and calculate moments
function calculate_moments(paths::SimPaths{T,S}, policy::PolicyFuncs{T,S}, vFuncs::VFuncs{T,S}, params::Params{T,S}, Rl::T, trim::S) where {T<:Real,S<:Integer}

    endT = Int(paths.bigT - 3); # trim the last 3 periods to avoid the terminal condition issue

    ############### target moments ############################
    debtToLiability = mean(paths.debtToLiability[trim:endT, :, :])
    loanToAsset = mean(paths.loanToAsset[trim:endT, :, :])
    # loanToAsset_afterLambdaPrime = mean(paths.loanToAsset[:, :, trim:paths.bigT])
    loanRate = Rl
    # capital to deposit conditional on no failure; by construction, capitalToDeposit = 0 if failure 
    A = paths.capitalToDeposit[trim:endT, :, :];
    nz = count(!iszero, A);
    if nz == 0
        capitalToDeposit = NaN
    else
        s = sum(x for x in A if x != 0; init =0.0)
        capitalToDeposit = s/nz;
    end 

    ################# non-targeted moments ######################
    # dividend to deposit 
    div = paths.divSim[trim:endT, :, :];
    deposit = paths.deltaSim[trim:endT, :, :];
    A = div./deposit;
    divToDeposit = mean(A)

    div_pos = [x for x in A if x > 0.0] 
    div_neg = [x for x in A if x < 0.0]

    divToDeposit_pos = mean(div_pos)
    divToDeposit_neg = mean(div_neg)

    failureRate = mean(paths.failureSim[trim:endT, :, :])

    if failureRate == 0.0
        nPrimeUnderFailure = NaN
    else
        nPrime_failure = paths.nPrimeSim[trim:endT, :, :] .* paths.failureSim[trim:endT, :, :]
        n = count(!iszero, nPrime_failure)
        s = sum(nPrime_failure)
        nPrimeUnderFailure = s/n
    end

    nPrime_noFailure = paths.nPrimeSim[trim:endT, :, :] .* (1 .- paths.failureSim[trim:endT, :, :])
    n = count(!iszero, nPrime_noFailure)
    s = sum(nPrime_noFailure)
    nPrimeUnderNoFailure = s/n
    loanToDeposit_temp = (paths.lSim[trim:endT, :, :] .+ params.g * paths.deltaSim[trim:endT, :, :]) ./ paths.deltaSim[trim:endT, :, :]
    loanToDeposit = mean(loanToDeposit_temp)

    mom = [debtToLiability, loanToAsset, capitalToDeposit, failureRate, loanRate]
    other_mom = [divToDeposit, divToDeposit_pos, divToDeposit_neg, nPrimeUnderFailure, nPrimeUnderNoFailure, loanToDeposit]
    return (moments = mom, other_moments = other_mom)
end

function simulate_and_moments(params::Params{T,S},params_cal::Params_cal{T,S},vFuncs::VFuncs{T,S},policy::PolicyFuncs{T,S},Rl::T,regime::F) where {T<:Real,S<:Integer,F<:Bool}

    bigT = params_cal.bigT;
    bigN = params_cal.bigN;
    bigJ = params_cal.bigJ;
    trim = params_cal.trim;

    shocks = makeShockSequence(params,vFuncs,bigT,bigN,bigJ);
    paths = Initiate_Paths(shocks, params);
    simulation_false = Simulate_paths(paths, policy, vFuncs, params, Rl, trim,regime);
    return (shocks = shocks, paths = paths, moments = simulation_false);
end

function solve_simulate_and_moments(params::Params{T,S},params_cal::Params_cal{T,S},regime::F) where {T<:Real,S<:Integer,F<:Bool}

    # 1. solve the model
    println("Solving the model...")
    sol = solve_model(params, regime, params_cal.a, params_cal.b);
    println("Solved the model.")
    policy = Get_PolicyFuncs(params, sol.eqq.Iterobj_is, sol.eqq.vFuncs, sol.Rl_star); # get policy functions from the solution
    println("Obtained policy functions.")
    simulation = simulate_and_moments(params,params_cal,sol.eqq.vFuncs,policy,sol.Rl_star,regime);
    @show simulation.moments.moments
    return (moments = simulation.moments, Rl_star = sol.Rl_star, sol = sol, policy = policy, paths = simulation.paths);
end



#### for more delicate calibration ###### 
function calibration(params::Params{T,S},params_cal::Params_cal{T,S},regime::F) where {T<:Real,S<:Integer,F<:Bool}

    target_moments = [
        params_cal.debt_to_liability,
        params_cal.loan_to_asset,
        params_cal.capital_to_deposit,
        params_cal.loan_rate]; #target moments from data

    calibrated_params = [
        params_cal.cM,
        params_cal.cO,
        params_cal.cL,
        params_cal.E,
        params_cal.dBar];
        # params_cal.σ]; #parameters to be calibrated  # params_cal.ϵ,
    
    function loss_function_given_params(x::Vector{T},params::Params{T,S},params_cal::Params_cal{T,S},Regime::F) where {T<:Real,S<:Integer,F<:Bool}
        
        # update params with calibrated parameters
        cM = x[1]; 
        cO = x[2];
        cL = x[3];
        E = x[4];
        dBar = x[5];
       #  σ = x[7];
       #   ϵ = x[4];

        params_candidate = Params{T,S}(
            params.qd,
            params.β,
            params.Rf,
            params.wr,
            params.α,
            params.ρ,
            params.g,
            params.ξ,
            params.cF,
            dBar,
            params_cal.σ,
            params.τC,
            params.z,
            params.δL,
            params.δM,
            params.δH,
            cM,
            cO,
            cL,
            params_cal.ϵ,
            E,
            params.H,
            params.F,
            params.M,
            params.λL,
            params.λM,
            params.λH,
            params.γ,
            params.ϕ,
            params.deltaGrid,
            params.lambdaGrid,
            params.nGrid,
            params.lGrid,
            params.sGrid,
            params.bGrid); 

        println("parameter candidates: cM = $cM , cO = $cO, cL = $cL, E = $E, dBar = $dBar ");
        println("1) solve_model block begins")
        @time sol = solve_model(params_candidate, Regime, params_cal.a, params_cal.b); # solve the model with candidate parameters
        println("2) Get_polciyFuncs block begins")
        @time policy = Get_PolicyFuncs(params_candidate, sol.eqq.Iterobj_is, sol.eqq.vFuncs, sol.Rl_star); # get policy functions from the solution
        println("3) simulate_and_moments block begins")
        @time sim = simulate_and_moments(params_candidate,params_cal,sol.eqq.vFuncs,policy,sol.Rl_star,Regime);
        @show msim.moments.moments
        @show sim.moments.other_moments

        println("4) Calculating the loss")
        @show sum((model_moments .- target_moments).^2)
        return sum((model_moments .- target_moments).^2)
    end

    # calibrated parameters: cM, cO, cL, E, dBar
    lb = zeros(length(calibrated_params)); # lower bound
    ub = zeros(length(calibrated_params)); # upper bound
    for i in [1]
        lb[i] = 0.25 * calibrated_params[i]; # 하한: 기준값의 25%
        ub[i] = 1.75 * calibrated_params[i]; # 상한: 기준값의 175%
    end

    for i in [2, 3, 4, 5]
        lb[i] = calibrated_params[i]; 
        ub[i] = calibrated_params[i]; 
    end

    @show lb
    @show ub


    n_sample = 5;
    s = SobolSeq(lb, ub);
    skip(s, 100_100); # skip the first 100,000 points
    x = next!(s);
    px = hcat([next!(s) for i in 1:n_sample]...)'; # generate matrix: n_sample x length(calibrated_params)
    PX = [px[i,:] for i in 1:size(px,1)]

    loss_temp(x) = loss_function_given_params(x, params, params_cal, regime)

    @time Parallel_test = pmap(loss_temp, PX)
    best_idx = argmin(Parallel_test);
    best_initial = vec(px[best_idx, :]);
    println("Best initial guess:", best_initial)

    # @time result = optimize(loss_function, lb, ub, best_initial, Fminbox(BFGS()))
    # println("Final estimated parameters:", Optim.minimizer(result))
    # println("Loss at optimum:", Optim.minimum(result))

    # return (best_loss = Parallel_test[best_idx], best_initial = best_initial, calibrated_params = result)
    return (result = Parallel_test, best_idx = best_idx, best_initial = best_initial)
end
