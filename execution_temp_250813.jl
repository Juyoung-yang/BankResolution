
pwd()

# include("C:\\Users\\master\\Desktop\\(2025) banking regulation\\BankResolution\\main.jl")
include("C:\\Users\\master\\Desktop\\(2025) banking regulation\\BankResolution\\parameters.jl")
include("/Users/juyoungyang_kdi/BankResolution/parameters.jl")
include("/Users/juyoungyang_kdi/BankResolution/main.jl")

using JuMP, Ipopt, Interpolations, MathOptInterface, LinearAlgebra, ForwardDiff, Plots, Ipopt, Roots, Statistics
using Profile, ProfileView, QuantEcon

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
        ϵ::T

        H  # transition matrix for delta: [delta, deltaPrime] so that sum(H[1,:]) = 1 
        F  # transition matrix for lambda: [lambdaPrime, lambda] so that sum(G[:,1]) = 1
        M::T # bank mass, policy tool 

        λL::T
        λM::T
        λH::T
        γ::T
        ϕ::T

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

function Initiate_Params(qd::T,β::T,Rf::T,wr::T,α::T,ρ::T,g::T,ξ::T,cF::T,dBar::T,σ::T,τC::T,z::T,α1::T,α2::T,α3::T,δL::T,δM::T,δH::T,cM::T,cO::T,cL::T,ϵ::T,H::Array{T,2},F::Array{T,2},M::T,λL::T,λM::T,λH::T,γ::T,ϕ::T,n_start::T,n_npts::S,n_stop::T,l_start::T,l_npts::S,l_stop::T,s_start::T,s_npts::S,s_stop::T,b_start::T,b_npts::S,b_stop::T) where {T<:Real,S<:Integer}
        deltaGrid = [δL,δM,δH] # Define a Tuple, immutable 
        lambdaGrid = [λL,λM,λH] # regular array, mutable
        nGrid = range(n_start,stop=n_stop,length=n_npts)
    
        lGrid = range(l_start,stop=l_stop,length=l_npts)
        sGrid = range(s_start,stop=s_stop,length=s_npts)
        bGrid = range(b_start,stop=b_stop,length=b_npts)
    
        pam = Params{T,S}(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,α1,α2,α3,δL,δM,δH,cM,cO,cL,ϵ,H,F,M,λL,λM,λH,γ,ϕ,deltaGrid,lambdaGrid,nGrid,lGrid,sGrid,bGrid)
        return pam
end

function Initiate_vFunc(params::Params{T,S}) where {T<:Real,S<:Integer}
        VF = ones(3, 3, length(params.nGrid)); # expected value function
        qBond = zeros(length(params.lGrid), length(params.nGrid), length(params.bGrid), 3, 3); # bank bond price schedule
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

function Initiate_IterObj_i(params::Params{T,S}) where {T<:Real,S<:Integer, F<:Bool}
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

function inter_v_temp(params::Params{T,S}, vFuncs::VFuncs{T,S}, delta::T,lambda::T,n::T)::T  where {T<:Real, S<:Integer, F<:Bool}
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
                @show iδ, iλ

                if failure 
                    println("failure at (iδ, iλ): ", iδ, iλ)
                        n_temp = n_failure(l,λprime)
                        V_temp[iδ, iλ] = regime ? zero(T) : (1 - params.ρ) * inter_v_temp(params, vFuncs, δprime, λprime, n_temp)
                else
                    println("success at (iδ, iλ): ", iδ, iλ)
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
        EV_conditional = dot(params.H[iDelta, :], V_temp * params.F[:, iLambda]) # the return is a scalar 
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

function NAV(l::T,s::T,b::T,lambda::T)::T where {T<:Real, S<:Integer, F<:Bool} 
        return Rl*( (1-lambda)*l + params.g*δ) + (1+params.Rf)*s - δ - b # net asset value of the bank, lambda here is lambda prime 
end

function tax(l::T,s::T,b::T,lambda::T)::T where {T<:Real, S<:Integer, F<:Bool} 
        return params.τC * max(0, (Rl-1)*((1-lambda)*l +params.g*δ) + params.Rf*s - params.Rf *(δ + b)) # tax on the bank's asset value
end

function n_failure(l::T, lambda::T)::T where {T<:Real, S<:Integer, F<:Bool} 
        return params.α * params.wr * Rl * (1-lambda)*l # next period asset conditional on bank failure 
end

function n_success(l::T,s::T,b::T,lambda::T)::T where {T<:Real, S<:Integer, F<:Bool} 
        return NAV(l,s,b,lambda) - tax(l,s,b,lambda) # next period asset conditional on bank success 
end

psi(params::Params, d) = d >= 0 ? (d + params.dBar)^params.σ  - params.dBar^params.σ : 1 - exp(-d)

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

function VFI_i(params::Params{T,S}, vFuncs::VFuncs{T,S}, vFuncsNew::VFuncsNew{T,S}, Rl::T, iterObj_i::IterObj_i{T,S}, iDelta::S, iLambda::S, regime::F) where {T<:Real, S<:Integer, F<:Bool}
        δ = params.deltaGrid[iDelta]; # get the state for Delta
        λ = params.lambdaGrid[iLambda]; # get the state for Lambda

        NAV(l::T,s::T,b::T,lambda::T)::T = Rl*( (1-lambda)*l + params.g*δ) + (1+params.Rf)*s - δ - b # net asset value of the bank, lambda here is lambda prime 
        tax(l::T,s::T,b::T,lambda::T)::T = params.τC * max(0, (Rl-1)*((1-lambda)*l +params.g*δ) + params.Rf*s - params.Rf *(δ + b)) # tax on the bank's asset value
        n_failure(l::T, lambda::T)::T = params.α * params.wr * Rl * (1-lambda)*l # next period asset conditional on bank failure 
        n_success(l::T,s::T,b::T,lambda::T)::T = NAV(l,s,b,lambda) - tax(l,s,b,lambda) # next period asset conditional on bank success 
        # div_func(l::T,s::T,b::T,n::T)::T = n + (1-params.cL)*params.β*δ + vFuncs.qBond[l,s,b,iDelta,iLambda]*b - l - params.g * δ -s-params.cM*l^2 -params.cO 

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
        function gen_EV_temp(l::T,s::T,b::T,params::Params{T,S},vFuncs::VFuncs{T,S},regime::F)::T # incorporate failure array and conditional asset array to get ex-post asset array 
            V_temp = gen_V_temp(l,s,b,params,vFuncs, regime) # [nDelta, nLambda], evaluated value functions at (delta prime, lambda prime) when choosing l,s,b
            EV_conditional = dot(params.H[iDelta, :], V_temp * params.F[:, iLambda]) # the return is a scalar 
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
        function solve_bank_problem(params::Params{T,S},vFuncs::VFuncs{T,S},iDelta::S,iLambda::S,iN::S,G::Array{T,3}) where {T<:Real, S<:Integer}

            lRange = range(first(params.lGrid), last(params.lGrid), length(params.lGrid))
            sRange = range(first(params.sGrid), last(params.sGrid), length(params.sGrid))
            bRange = range(first(params.bGrid), last(params.bGrid), length(params.bGrid))
    
            Gmin = minimum(G)
            Gmax = maximum(G)
            Gnorm = (G .- Gmin) ./ (Gmax - Gmin) # scaling 
            G_itp_rev = interpolate(Gnorm, BSpline(Quadratic(Line(OnGrid()))))
            G_itp_rev = Interpolations.scale(G_itp_rev, lRange, sRange, bRange)
            G_itp_ext_rev = extrapolate(G_itp_rev, Line())
            G_interp_rev(l, s, b) = G_itp_ext_rev(l, s, b)
    
            model = Model(Ipopt.Optimizer)
            set_optimizer_attribute(model, "tol", 1e-6)
            set_optimizer_attribute(model, "acceptable_tol", 1e-3)
            set_optimizer_attribute(model, "max_iter", 5000)
            set_optimizer_attribute(model, "acceptable_iter", 10)
            set_optimizer_attribute(model, "print_level", 0)
            set_optimizer_attribute(model, "hessian_approximation", "limited-memory")
         
            JuMP.register(model, :G_interp_rev, 3, G_interp_rev; autodiff = true)
            l_min, l_max = first(params.lGrid), last(params.lGrid) # (last(params.lGrid)+first(params.lGrid))/2, last(params.lGrid)
            s_min, s_max = first(params.sGrid), last(params.sGrid)
            b_min, b_max = first(params.bGrid), last(params.bGrid)
            l_start = (l_max+l_min)/2
            s_start = (s_max+s_min)/2
            b_start = (b_max+b_min)/2
            # println("c1 residual = ", (1-params.α*params.wr)*l_start + s_start - b_start - (1-params.g)*params.deltaGrid[iDelta])
            # println("c2 residual = ", l_start - (params.β - params.g)*params.deltaGrid[iDelta])
            @variable(model, l_min <= l <= l_max, start = l_start)
            @variable(model, s_min <= s <= s_max, start = s_start)
            @variable(model, b_min <= b <= b_max, start = b_start)
            # println("G_interp at start = ", G_interp_rev(l_start, s_start, b_start))
            # grad = ForwardDiff.gradient(u -> G_interp_rev(u[1], u[2], u[3]), [l_start, s_start, b_start])
            # println("gradient at start = ", grad)
    
            @NLobjective(model, Max, G_interp_rev(l, s, b))
            @constraint(model, (1-params.α*params.wr)*l + s - b >= (1-params.g)*params.deltaGrid[iDelta]) # constratint 1
            @constraint(model, l <= (params.β - params.g)*params.deltaGrid[iDelta]) # constratint 2
    
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
            else
                error("Optimization did not converge: status = $stat")
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
            iterObj_i.failure[iN] = [NAV(sol[1], sol[2], sol[3], lamPrime) < 0 for lamPrime in params.lambdaGrid] # if true, bank fails 
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
            println("iter=", iter, ", maxdiffs=", maxdiffs); # report the iteration progress 
            end

            # 3-3. if the difference is not small enough, do the iteration again
            iter += 1;
        end

        # 4. finish the iteration and return the value function 
        return (vFuncs = vFuncs, vFuncsNew = vFuncsNew, Iterobj_is = Iterobj_is);  
end

function Update_vFuncs_Diffs(vFuncs::VFuncs{T,S}, vFuncsNew::VFuncsNew{T,S}, params::Params{T,S}) where {T<:Real, S<:Integer}; # after the optimization, calculate the differene

        # 1. calculate the difference of VF_0 and VF_1 
        vFuncsNew.diffs .= vFuncsNew.VF - vFuncs.VF; # difference between VF and VF_new
        diffs = norm(vFuncsNew.diffs, Inf); # infinity norm of the difference of VF

        # 2. update vFuncsNew objects to vFuncs
        copyto!(vFuncs.VF, vFuncsNew.VF); # update EVF
        copyto!(vFuncs.qBond, vFuncsNew.qBond) # update qBond
        copyto!(vFuncs.X, vFuncsNew.X) # update qBond
        return diffs
end

function qBond_condiState(params::Params{T,S}, Rl::T, il::S ,is::S, ib::S, iDelta::S) where {T<:Real,S<:Integer} 
        
        l, s, b, Delta = params.lGrid[il], params.sGrid[is], params.bGrid[ib], params.deltaGrid[iDelta]
        lambdaStar(l,s,b,Delta) = (Rl*(l+params.g*Delta)+(1+params.Rf)*s-Delta-b)/(Rl*l)
    
        function return_temp_underFailure(l,s,b,Delta,Lambda)::T
            b == 0 && return zero(T)
            term = Rl*((1-Lambda)*l+params.g*Delta)+(1+params.Rf)*s 
            num = max(min(b,max(zero(T),params.cF*term-Delta)),term-Delta-params.α*params.wr*Rl*(1-Lambda)*l)
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
    tax(l::T,s::T,b::T,lambda::T)::T = params.τC * max(0, (Rl-1)*((1-lambda)*l +params.g*δ) + params.Rf*s - params.Rf *(δ + b)) # tax on the bank's asset value
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
        vFuncsNew.Γ[:, iLambdaPrime, 1] .+= params.H[iDelta, :] .* vFuncs.Γ[iDelta, iLambda, iN] .* params.F[iLambdaPrime,iLambda]
    elseif searchsortedlast(params.nGrid, nPrime) == length(params.nGrid) # if n is larger than the maximum nGrid
        vFuncs.Γ[:, iLambda, length(params.nGrid)] .+= params.H[iDelta, :] .* vFuncs.Γ[iDelta, iLambda, iN] .* params.F[iLambdaPrime,iLambda]
    else # set in between
        iNn, iNnPrime = searchsortedlast(params.nGrid, nPrime), searchsortedfirst(params.nGrid, nPrime)
        vFuncsNew.Γ[:,iLambdaPrime,iNnPrime] .+= params.H[iDelta, :] .* ((nPrime - params.nGrid[iNn]) / (params.nGrid[iNnPrime] - params.nGrid[iNn])) * vFuncs.Γ[iDelta, iLambda, iN] .* params.F[iLambdaPrime,iLambda]
        vFuncsNew.Γ[:,iLambdaPrime,iNn] .+= params.H[iDelta, :] .* ((params.nGrid[iNnPrime] - nPrime) / (params.nGrid[iNnPrime] - params.nGrid[iNn])) * vFuncs.Γ[iDelta, iLambda, iN] .* params.F[iLambdaPrime,iLambda]
    end

end

function stationary_distribution(params::Params{T,S}, Rl::T, vFuncs::VFuncs{T,S},vFuncsNew::VFuncsNew{T,S},Iterobj_is::Matrix{IterObj_i{T,S}}, maxiter::S, tol::T) where {T<:Real,S<:Integer}

    iter, maxdiffs, diffs = 1, one(T), zeros(T);

    while iter <= maxiter && maxdiffs > tol

        println("iter = ", iter, ", maxdiffs = ", maxdiffs); # report the iteration progress
        Threads.@threads for iDelta in 1:length(params.deltaGrid) # 3: number of states for Delta
         Threads.@threads for iLambda in 1:length(params.lambdaGrid)  # 3: number of states for Lambda
            Threads.@threads for iN in 1:length(params.nGrid)
                Threads.@threads for iLambdaPrime in 1:length(params.lambdaGrid)
                        Update_stationary_dist(params, vFuncs, vFuncsNew, Iterobj_is, iDelta, iLambda, iN, iLambdaPrime, Rl)
                    end
                end
            end
        end
        
        vFuncsNew.Γ ./= sum(vFuncsNew.Γ) # normalize the distribution
        diffs = Update_vFuncs_Diffs_stationary(vFuncs, vFuncsNew, params); # after updating the mass, calculate the difference
        maxdiffs = maximum(diffs);

        if mod(iter, 200) == 0
            println("iter=", iter, ", maxdiffs=", maxdiffs); # report the iteration progress 
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
    println("Regime is ", Regime) 
    eq = VFI(Params, Rl, Regime, 1000, 1e-4); # run the VFI algorithm with the given parameters and regime
    println("Regime is ", Regime) 
    eqq = stationary_distribution(Params, Rl, eq.vFuncs, eq.vFuncsNew, eq.Iterobj_is, 1000, 1e-4); # run the stationary distribution algorithm with the given parameters and regime
    excess_loan_supply = aggre_loan_supply(Params, eqq.vFuncs, eqq.Iterobj_is) - Rl^(Params.ϵ)
    return excess_loan_supply
end

# output: equilibrium loan rate
function solve_model(params::Params{T,S}, regime::F, a::T, b::T) where {T<:Real,S<:Integer,F<:Bool}
    
    solve_model_given_r_single = Rl -> solve_model_given_r(Rl; Params = params, Regime = regime)
    println("Regime is ", regime) 

    Rl_star = find_zero(solve_model_given_r_single, (a, b), Bisection(); tol=1e-8, maxevals=100);
    println("Rl_star", Rl_star) 

   # fa = solve_model_given_r_single(a)
   # fb = solve_model_given_r_single(b)
   # @assert fa * fb < 0 "The function must have opposite signs at the endpoints a and b."

    # obtain the equilibrium under equilibrium loan rate Rl_star
    eq =  VFI(params, Rl_star, regime, 1000, 1e-4);
    eqq = stationary_distribution(params, Rl_star, eq.vFuncs, eq.vFuncsNew, eq.Iterobj_is, 1000, 1e-4);
    return (Rl_star = Rl_star, eq = eq, eqq = eqq)
end

################################################################################################################################################
## code for simulation and calculating moments
struct PolicyFuncs{T,S} where {T<:Real,S<:Integer}
    lPolicy::Array{T,3} # loan policy function
    sPolicy::Array{T,3} # savings policy function
    bPolicy::Array{T,3} # bond policy function
    failure::Array{Bool,4} # failure or not
end

struct SimShocks{T,S} where {T<:Real,S<:Integer}
    deltaIndSim::Array{S,3}
    lambdaIndSim::Array{S,3}
end

struct SimPaths{T,S} where {T<:Real,S<:Integer}
    deltaIndSim::Array{S,3}
    lambdaIndSim::Array{S,3}
    bigT::S # number of periods (T, year)
    bigN::S # number of observations (N)
    bigJ::S # number of banks (J)
    deltaSim::Array{T,3} # delta states for each bank J
    lambdaSim::Array{T,3} # lambda states for each bank J
    lSim::Array{T,3} # loan 
    sSim::Array{T,3} # safe-asset saving
    bSim::Array{T,3} # bond issuance, bank borrowing 
    failureSim::Array{S,3} # failure or not, denoted as 1 or 0
    nSim::Array{T,3} # net asset value, current period 
    nIndSim::Array{S,3} # net asset value (index), current period 
    nPrimeSim::Array{T,3} # net asset value, next period with lambdaPrime
    nPrimeIndSim::Array{S,3} # net asset value (index), next period with lambdaPrime
    assetSim::Array{T,3} # asset, current period
    assetPrimeSim::Array{T,3} # asset, next period with lambdaPrime
    Γ::Array{T,3} # distribution of states, current period
    ΓPrime::Array{T,3} # distribution of states, next period with lambdaPrime
end

################################################################################################################################################

# CHECK if MarkovChaine and simulate are defined in the package: Need to use QuantEcon.jl
function makeShockSequence(params::Params{T,S}, bigT::S, bigN::S, bigJ::S) where {T<:Real,S<:Integer}

    deltaIndSim = zeros(Int, bigT, bigN, bigJ); # delta states for each bank J
    lambdaIndSim = zeros(Int, bigT, bigN, bigJ); # lambda states for each bank J

    mcDelta = QuantEcon.MarkovChain(params.H) # Create a Markov chain with transition probabilities
    mcLambda = QuantEcon.MarkovChain(params.F)

    for smallJ = 1:bigJ # for bank J
        for smallN = 1:bigN
            deltaIndSim[:, smallN, smallJ] .= QuantEcon.simulate(mcDelta, bigT); # Simulate delta states;
            lambdaIndSim[:, smallN, smallJ] .= QuantEcon.simulate(mcLambda, bigT); # Simulate lambda states;
        end
    end

    return SimShocks(deltaIndSim, lambdaIndSim)
end


function Initiate_Paths(shocks::SimShocks{T,S}, params::Params{T,S}) where {T<:Real,S<:Integer}

    deltaIndSim = shocks.deltaIndSim;
    lambdaIndSim = shocks.lambdaIndSim;

    bigT = size(deltaIndSim, 1); # number of periods (T, year)
    bigN = size(deltaIndSim, 2); # number of observations (N)
    bigJ = size(deltaIndSim, 3); # number of banks

    # below are variables needed to be filled
    deltaSim = zeros(T, bigT, bigN, bigJ);
    lambdaSim = zeros(T, bigT, bigN, bigJ);
    lSim = zeros(T, bigT, bigN, bigJ); # loan policy function
    sSim = zeros(T, bigT, bigN, bigJ); # savings policy function
    bSim = zeros(T, bigT, bigN, bigJ); # bond policy function
    failureSim = zeros(S, bigT, bigN, bigJ); # failure or not, denoted as 1 or 0
    nSim = zeros(T, bigT, bigN, bigJ); # next period asset
   # nIndSim = zeros(S, bigT, bigN, bigJ);
    nPrimeSim = zeros(T, bigT, bigN, bigJ);
   # nPrimeIndSim = zeros(S, bigT, bigN, bigJ);
    assetSim = zeros(T, bigT, bigN, bigJ); 
    assetPrimeSim = zeros(T, bigT, bigN, bigJ);

    Γ = zeros(3, 3, length(params.nGrid)); # stationary distribution
    ΓPrime = zeros(3, 3, length(params.nGrid)); # next period stationary distribution

    return SimPaths(deltaIndSim, lambdaIndSim, bigT, bigN, bigJ, deltaSim, lambdaSim, lSim, sSim, bSim, failureSim, nSim, nIndSim, nPrimeSim, nPrimeIndSim, assetSim, assetPrimeSim, Γ, ΓPrime)
end

################################################################################################################################################


# bank policy functions: (l,s,b), failure or not, 
function Initiate_PolicyFuncs(params::Params{T,S}) where {T<:Real,S<:Integer}

    return PolicyFuncs
end

function Get_PolicyFuncs()
    e
end

################################################################################################################################################

# after initiating paths, 
function Simulate_paths(paths::SimPaths{T,S}, policy::PolicyFuncs{T,S}, vFuncs::VFuncs{T,S}, params::Params{T,S}, Rl::T, trim::S) where {T<:Real,S<:Integer}

    bigT = size(paths.deltaIndSim, 1); # number of periods (T, year)
    bigN = size(paths.deltaIndSim, 2); # number of observations (N)
    bigJ = size(paths.deltaIndSim, 3); # number of banks


    # 1. for bank smallJ calculate the series of interest at time (smallT, smallN) or given state at (smallT, smallN)
    for smallJ = 1:bigJ
        for smallN = 1:bigN
            for smallT = 1:bigT-1
                calculate_series(paths, policy, vFuncs, params, Rl, smallT, smallN, smallJ)
            end
        end
    end

    # 2. after obtaining simulation paths, trim some first parts and calculate moments
    moments = calculate_moments(paths, policy, vFuncs, params, Rl, trim);

    # 3. return key moments of interest
    return moments 
end

# for bank smallJ, calculate the variable of interest at time (smallT, smallN); given state (delta, lambda, n), calculate the series of interest
function calculate_series(paths::SimPaths{T,S}, policy::PolicyFuncs{T,S}, vFuncs::VFuncs{T,S}, params::Params{T,S}, Rl::T, smallT::S, smallN::S, smallJ::S) where {T<:Real,S<:Integer}

    NAV(l::T,s::T,b::T,lambda::T,δ::T)::T = Rl*( (1-lambda)*l + params.g*δ) + (1+params.Rf)*s - δ - b; # 순자산, net asset value 
    totalAsset(l::T,s::T,lambda::T,δ::T)::T = Rl * ( (1-lambda)*l + params.g*δ) + (1+params.Rf)*s; # total asset of the bank, 총자산
    tax(l::T,s::T,b::T,lambda::T,δ::T)::T = params.τC * max(0, (Rl-1)*((1-lambda)*l +params.g*δ) + params.Rf*s - params.Rf *(δ + b)) # tax on the bank's asset value
    n_success(l::T,s::T,b::T,lambda::T,δ::T)::T = NAV(l,s,b,lambda,δ) - tax(l,s,b,lambda,δ) # next period asset conditional on bank success 
    n_failure(l::T, lambda::T)::T = params.α * params.wr * Rl * (1-lambda)*l # next period asset conditional on bank failure 

    delta = params.deltaGrid[paths.deltaIndSim[smallT, smallN, smallJ]]; # delta state for bank J at time T
    paths.deltaSim[smallT, smallN, smallJ] = delta;
    lambda = params.lambdaGrid[paths.lambdaIndSim[smallT, smallN, smallJ]]; # lambda state for bank J at time T
    paths.lambdaSim[smallT, smallN, smallJ] = lambda;
    lambdaPrime = params.lambdaGrid[paths.lambdaIndSim[smallT+1, smallN, smallJ]];
    
    if smallT = 1 # initial period 
        paths.nSim[smallT, smallN, smallJ] = params.nGrid[1]; ##### start with first grid n point
        paths.lSim[smallT, smallN, smallJ] = policy.lPolicy[paths.deltaIndSim[smallT, smallN, smallJ], paths.lambdaIndSim[smallT, smallN, smallJ], 1]; # loan policy function
        paths.sSim[smallT, smallN, smallJ] = policy.sPolicy[paths.deltaIndSim[smallT, smallN, smallJ], paths.lambdaIndSim[smallT, smallN, smallJ], 1]; # savings policy function
        paths.bSim[smallT, smallN, smallJ] = policy.bPolicy[paths.deltaIndSim[smallT, smallN, smallJ], paths.lambdaIndSim[smallT, smallN, smallJ], 1]; # bond policy function
        paths.assetSim[smallT, smallN, smallJ] = totalAsset(paths.lSim[smallT, smallN, smallJ], paths.sSim[smallT, smallN, smallJ], lambdaPrime, delta);  # total asset given lambdaPrime
        paths.failureSim[smallT, smallN, smallJ] = ( paths.assetSim[smallT, smallN, smallJ] - paths.bSim[smallT, smallN, smallJ] - delta < 0 ) ? 1 : 0; # failure or not, denoted as 1 or 0

        if paths.failureSim[smallT, smallN, smallJ] == 1 # if failure
            paths.nPrimeSim[smallT, smallN, smallJ] = n_failure(paths.lSim[smallT, smallN, smallJ], lambdaPrime); # next period asset conditional on bank failure 
        else # if success
            paths.nPrimeSim[smallT, smallN, smallJ] = n_success(paths.lSim[smallT, smallN, smallJ], paths.sSim[smallT, smallN, smallJ], paths.bSim[smallT, smallN, smallJ], lambdaPrime, delta); # next period asset conditional on bank success 
        end

    else
        paths.nSim[smallT, smallN, smallJ] = paths.nPrimeSim[smallT-1, smallN, smallJ] # 잉여금 as a embodied state 
        paths.lSim[smallT, smallN, smallJ] = policy.lPolicy[paths.deltaIndSim[smallT, smallN, smallJ], paths.lambdaIndSim[smallT, smallN, smallJ], ---]; # loan policy function
        paths.sSim[smallT, smallN, smallJ] = policy.sPolicy[paths.deltaIndSim[smallT, smallN, smallJ], paths.lambdaIndSim[smallT, smallN, smallJ], ---]; # savings policy function
        paths.bSim[smallT, smallN, smallJ] = policy.bPolicy[paths.deltaIndSim[smallT, smallN, smallJ], paths.lambdaIndSim[smallT, smallN, smallJ], ---]; # bond policy function
        paths.assetSim[smallT, smallN, smallJ] = totalAsset(paths.lSim[smallT, smallN, smallJ], paths.sSim[smallT, smallN, smallJ], lambdaPrime, delta);  # total asset given lambdaPrime
        paths.failureSim[smallT, smallN, smallJ] = ( paths.assetSim[smallT, smallN, smallJ] - paths.bSim[smallT, smallN, smallJ] - delta < 0 ) ? 1 : 0; # failure or not, denoted as 1 or 0

        if paths.failureSim[smallT, smallN, smallJ] == 1 # if failure
            paths.nPrimeSim[smallT, smallN, smallJ] = n_failure(paths.lSim[smallT, smallN, smallJ], lambdaPrime); # next period asset conditional on bank failure 
        else # if success
            paths.nPrimeSim[smallT, smallN, smallJ] = n_success(paths.lSim[smallT, smallN, smallJ], paths.sSim[smallT, smallN, smallJ], paths.bSim[smallT, smallN, smallJ], lambdaPrime, delta); # next period asset conditional on bank success 
        end
    end

    # calculate the bank distribution across states
    if smallT = 1 
        paths.Γ[smallT, smallN, smallJ] = 
        paths.ΓPrime[smallT, smallN, smallJ] = 
    else
        paths.Γ[smallT, smallN, smallJ] = paths.ΓPrime[smallT-1, smallN, smallJ]
        paths.ΓPrime[smallT, smallN, smallJ] = 
    end
end

# after obtaining simulation paths, trim some first parts and calculate moments
function calculate_moments(paths::SimPaths{T,S}, policy::PolicyFuncs{T,S}, vFuncs::VFuncs{T,S}, params::Params{T,S}, Rl::T, trim::S) where {T<:Real,S<:Integer}
    d
end

################################################################################################################################################


params = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,α1,α2,α3,δL,δM,δH,cM,cO,cL,ϵ,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
shocks = makeShockSequence(params, 100, 1000, 10); # make shock sequence for simulation
paths = Initiate_Paths(shocks, params); # initiate paths for simulation


@time sol = solve_model(params, false, 0.01, 1.9)
sol.eq.vFuncs == sol.eqq.vFuncs
mean(sol.eqq.vFuncs.qBond)
aggre_loan_supply2(params, sol.eqq.vFuncs, sol.eqq.Iterobj_is)
sol.eq.vFuncsNew.Γ
sol.eqq
sol.eq.Iterobj_is[1,1].solution
sol.eq.Iterobj_is[1,2].solution
sol.eq.Iterobj_is[1,3].solution
sol.eq.Iterobj_is[3,3].solution
sol.eq.Iterobj_is[2,3].solution

params.nGrid[5]
params.deltaGrid[1]






@time sol_true = solve_model(params, true, 0.01, 1.9)
sol_true.eq.vFuncs.qBond
mean(sol_true.eq.vFuncs.qBond)
aggre_loan_supply2(params, sol_true.eqq.vFuncs, sol_true.eqq.Iterobj_is)
sol_true.eq.vFuncsNew.Γ


Simulate_paths(paths, Get_PolicyFuncs(), sol.eq.vFuncs, params, ###); 

################################################################################################################################################


Rl_star = 0.028432711437424655
@time eq_example = VFI(params, Rl, regime, 1000, 1e-4)
@time eqq = stationary_distribution(params, Rl, eq_example.vFuncs, eq_example.vFuncsNew, eq_example.Iterobj_is, 1000, 1e-4)

eq_example.vFuncs.qBond
eq_example.vFuncsNew.qBond

eq_example.vFuncs.qBond .= params.β


# function solve_model(params::Params{T,S}, regime::F) where {T<:Real,S<:Integer,F<:Bool}
    Rl = 1.0; # initial guess for Rl
    maxiter = 1000; # maximum number of iterations
    tol = 1e-4; # tolerance for convergence

    for iter in 1:maxiter
        println("Iteration: ", iter, ", Rl: ", Rl);
        excess_loan_supply = solve_model_given_r(Rl; params=params, regime=regime); # solve the model given Rl
        if abs(excess_loan_supply) < tol # check convergence
            println("Converged with Rl: ", Rl);
            break;
        end
        Rl += excess_loan_supply * 0.01; # update Rl based on excess loan supply
    end

    return (Rl = Rl, excess_loan_supply = excess_loan_supply);
# end

aggre_loan_supply(params, eqq.vFuncs, eqq.Iterobj_is)

################################################################################################################################################
Profile.clear()
@profile VFI(params, Rl, regime, 1000, 1e-4) 
Profile.print(format = :flat) # print the profile results
Profile.print(format = :hier)
ProfileView.view()
Profile.print(sortedby = :time, format = :flat)


regime = false; 
Rl = 1.0
@time eq = VFI(params, Rl, true, 1000, 1e-4) # run the VFI algorithm with the given parameters and regime
eq.vFuncs.qBond
eq.vFuncsNew.qBond

@time eqq = stationary_distribution(params, Rl, eq.vFuncs, eq.vFuncsNew, eq.Iterobj_is, 1000, 1e-4); # run the stationary distribution algorithm with the given parameters and regime
# eq.vFuncs.VF[:,:,10]
# eq.vFuncs.qBond
# eq.vFuncs.X
example = [eqq.Iterobj_is[a,b].solution[c][1] for a in 1:3, b in 1:3, c in 1:11]
example
eqq.Iterobj_is[2,3].solution
public_given_state = [params.g * params.deltaGrid[a] for a in 1:3]
example .+ public_given_state
ex = (example .+ public_given_state) .* eqq.vFuncsNew.Γ
ex_ = sum(ex)
# eqq.vFuncsNew.Γ[:,:,:] 
# eqq.vFuncs.Γ[:,:,:] 

VFI_i(params::Params{T,S}, vFuncs::VFuncs{T,S}, vFuncsNew::VFuncsNew{T,S}, Rl::T, iterObj_i::IterObj_i{T,S}, iDelta::S, iLambda::S, regime::F)

eq.Iterobj_is[1,1].solution
eq.Iterobj_is[1,2].solution
eq.Iterobj_is[2,1].solution
sum(eq.Iterobj_is[2,2].failure)


@time eq_special = VFI(params, Rl, true, 1000, 1e-6) # run the VFI algorithm with the given parameters and regime
eq_special.vFuncs.qBond[:,5,5,3,:]
sum(eq_special.Iterobj_is[2,2].failure) 

###############################################################################################################################################

params = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,α1,α2,α3,δL,δM,δH,cM,cO,cL,ϵ,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
vFuncs = Initiate_vFunc(params)
vFuncsNew = Initiate_vFuncNew(params)
iterObj_i = Initiate_IterObj_i(params)
iDelta = 3; iLambda = 1;  
sol2 = VFI_i(params, vFuncs, vFuncsNew, Rl, iterObj_i, iDelta, iLambda, regime)

####################################





















