
include("C:\\Users\\master\\Desktop\\(2025) banking regulation\\BankResolution\\parameters.jl")
# include("C:\\Users\\master\\Desktop\\(2025) banking regulation\\BankResolution\\main.jl")


function Initiate_vFunc(params::Params{T,S}) where {T<:Real,S<:Integer}
    VF = ones(3, 3, length(params.nGrid)); # expected value function
    qBond = zeros(length(params.lGrid), length(params.nGrid), length(params.bGrid), 3, 3); # bank bond price schedule
    Rl = 0; # loan interest rate
    X = zeros(3, 3, length(params.nGrid), 3); # bank's failure decision 

    return VFuncs{T,S}(VF, qBond, Rl, X)
end

function Initiate_IterObj_i(params::Params{T,S}) where {T<:Real,S<:Integer}
    EV = zeros(length(params.lGrid), length(params.sGrid), length(params.bGrid)); # expected value function conditional on a choice of (l,s,b)
    G = zeros(length(params.lGrid), length(params.sGrid), length(params.bGrid), length(params.nGrid)); # bank bond price schedule
    solution = zeros(T, length(params.nGrid))
    solution_index = solution_index = zeros(T, length(params.nGrid))
    failure = zeros(T, length(params.nGrid))

    return IterObj_i{T,S}(EV, G, solution, solution_index, failure)
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
        EV_conditional = dot(params.H[iDelta, :], V_temp * params.Γ[:, iLambda]) # the return is a scalar 
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

function solve_bank_problem(params::Params{T,S},vFuncs::VFuncs{T,S},iDelta::S,iLambda::S,iN::S,G::Array{T,3}) where {T<:Real, S<:Integer}
        # construct objective function int_G: interpolated function[l,s,b] using G as inputs
        G_itp = interpolate(G, BSpline(Quadratic(Line(OnGrid())))) 
        lRange = range(first(params.lGrid), last(params.lGrid), length(params.lGrid))
        sRange = range(first(params.sGrid), last(params.sGrid), length(params.sGrid))
        bRange = range(first(params.bGrid), last(params.bGrid), length(params.bGrid))
        G_itp = Interpolations.scale(G_itp, lRange, sRange, bRange)
        G_itp_ext = extrapolate(G_itp, Line())     # enables linear extrapolation
        G_interp(l, s, b) = G_itp_ext(l, s, b)

        model = Model(Ipopt.Optimizer)
        set_optimizer_attribute(model, "tol", 1e-2)
        set_optimizer_attribute(model, "max_iter", 10)
        JuMP.register(model, :G_interp, 3, G_interp; autodiff = true)
        l_min, l_max = first(params.lGrid), last(params.lGrid)
        s_min, s_max = first(params.sGrid), last(params.sGrid)
        b_min, b_max = first(params.bGrid), last(params.bGrid)
        @variable(model, l_min <= l <= l_max)
        @variable(model, s_min <= s <= s_max)
        @variable(model, b_min <= b <= b_max)

        @NLobjective(model, Max, G_interp(l,s,b)) # maximizing interpolated G(l,s,b)
        @NLconstraint(model, (l + params.g*params.deltaGrid[iDelta] +s - params.deltaGrid[iDelta]-b)/(params.wr*l) >= params.α) # constratint 1
        @NLconstraint(model, (l + params.g*params.deltaGrid[iDelta]) <= params.β*params.deltaGrid[iDelta]) # constratint 2

        optimize!(model)

        stat = termination_status(model)
        if stat == MathOptInterface.OPTIMAL || stat == MathOptInterface.LOCALLY_SOLVED
            return [value(l), value(s), value(b)]
        else
            error("Optimization did not converge: status = $stat")
        end
end

psi(params::Params, d) = d >= 0 ? (d + params.dBar)^params.σ  - params.dBar^params.σ : 1 - exp(-d)

function VFI_i(params::Params{T,S}, vFuncs::VFuncs{T,S}, vFuncsNew::VFuncsNew{T,S}, Rl::T, iterObj_i::IterObj_i{T,S}, iDelta::S, iLambda::S, regime::F) where {T<:Real, S<:Integer, F<:Bool}
    δ = params.deltaGrid[iDelta]; # get the state for Delta
    λ = params.lambdaGrid[iLambda]; # get the state for Lambda

    NAV(l::T,s::T,b::T,lambda::T)::T = Rl*( (1-lambda)*l + params.g*δ) + (1+params.Rf)*s - δ - b # net asset value of the bank, lambda here is lambda prime 
    tax(l::T,s::T,b::T,lambda::T)::T = params.τC * max(0, (Rl-1)*((1-lambda)*l +params.g*δ) + params.Rf*s - params.Rf *(δ + b)) # tax on the bank's asset value
    n_failure(l::T, lambda::T)::T = params.α * params.wr * Rl * (1-lambda)*l # next period asset conditional on bank failure 
    n_success(l::T,s::T,b::T,lambda::T)::T = NAV(l,s,b,lambda) - tax(l,s,b,lambda) # next period asset conditional on bank success 

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
        EV_conditional = dot(params.H[iDelta, :], V_temp * params.Γ[:, iLambda]) # the return is a scalar 
        return EV_conditional 
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
        # construct objective function int_G: interpolated function[l,s,b] using G as inputs
        G_itp = interpolate(G, BSpline(Quadratic(Line(OnGrid())))) 
        lRange = range(first(params.lGrid), last(params.lGrid), length(params.lGrid))
        sRange = range(first(params.sGrid), last(params.sGrid), length(params.sGrid))
        bRange = range(first(params.bGrid), last(params.bGrid), length(params.bGrid))
        G_itp = Interpolations.scale(G_itp, lRange, sRange, bRange)
        G_itp_ext = extrapolate(G_itp, Line())     # enables linear extrapolation
        G_interp(l, s, b) = G_itp_ext(l, s, b)

        model = Model(Ipopt.Optimizer)
        set_optimizer_attribute(model, "tol", 1e-8)
        JuMP.register(model, :G_interp, 3, G_interp; autodiff = true)
        l_min, l_max = first(params.lGrid) + eps(), last(params.lGrid)
        s_min, s_max = first(params.sGrid), last(params.sGrid)
        b_min, b_max = first(params.bGrid), last(params.bGrid)
        @variable(model, l_min <= l <= l_max)
        @variable(model, s_min <= s <= s_max)
        @variable(model, b_min <= b <= b_max)

        @NLobjective(model, Max, G_interp(l,s,b)) # maximizing interpolated G(l,s,b)
        @NLconstraint(model, (l + params.g*params.deltaGrid[iDelta] +s - params.deltaGrid[iDelta]-b)/(params.wr*l) >= params.α) # constratint 1
        @NLconstraint(model, (l + params.g*params.deltaGrid[iDelta]) <= params.β*params.deltaGrid[iDelta]) # constratint 2

        optimize!(model)

        stat = termination_status(model)
        if stat == MathOptInterface.OPTIMAL || stat == MathOptInterface.LOCALLY_SOLVED
            return [value(l), value(s), value(b)]
        else
            error("Optimization did not converge: status = $stat")
        end
    end

    function update_VF_with_solution(sol,params::Params{T,S},vFuncs::VFuncs{T,S},regime)::T
        l,s,b = sol[1], sol[2], sol[3]
        ev = gen_EV_temp(l,s,b,params,vFuncs,regime) 
        div_val = div(l,s,b)
        return div_val + params.β* ev
    end

    # for each n in nGrid, find the optimal (l, s, b), store iterobj_i.solution, iterobj_i.failure, and update vFuncs.VF
    for (iN, n) in pairs(params.nGrid) # for each n / Threads.@threads 
        println(iN), println(n)

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
        iterObj_i.failure[iN] .= NAV(sol[1], sol[2], sol[3], params.lambdaGrid) .<= 0 ? 1 : 0; # store the failure decision in iterObj_i.failure
        vFuncsNew.VF[iDelta, iLambda, iN] = update_VF_with_solution(sol,params,vFuncs,regime)
    end
end   

####################################
params = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,α1,α2,α3,δL,δM,δH,cM,cO,cL,H,Γ,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)
vFuncs = Initiate_vFunc(params)
iterObj_i = Initiate_IterObj_i(params)

iDelta = 3; iLambda = 1; iN = 7; regime = false; 
n = params.nGrid[iN]; δ = params.deltaGrid[iDelta]; lambda = params.lambdaGrid[iLambda]; Rl = 1.03;
@show n, δ, lambda, Rl, params.deltaGrid[iDelta]
l, s, b = 0.0, 70.0, 10.0
L = reshape(params.lGrid, :, 1, 1)    # nL × 1 × 1       
SS = reshape(params.sGrid, 1, :, 1)   # 1  × nS × 1
B = reshape(params.bGrid, 1, 1, :)    # 1  × 1  × nB
QB = @view vFuncs.qBond[:, :, :, iDelta, iLambda] 
@show typeof(QB)
@show params.lambdaGrid, params.deltaGrid, params.nGrid

function done()
    # examine the evaluation works: Okay!
    @show tax(l,s,-b,lambda), tax(l+10.0,s,b,lambda) # OKAY
    @show NAV(l,s, b, lambda), NAV(10.0, 1000.0, 1.0, 0.001) # OKAY
    @show n_failure(10.0, 0.02), n_failure(10.0, 0.2) # OKAY
    @show n_success(10.0, 1000.0, 1.0, 0.001), NAV(10.0, 1000.0, 1.0, 0.001), tax(10.0, 1000.0, 1.0, 0.001) # OKAY

    # examine the interpotation works 
    @show inter_v_temp(params, vFuncs, δ, lambda, 1000.234) # OKAY, able to evaluate for any n, even a negative one 
    @show inter_v_temp.([params], [vFuncs], [δ], [lambda], collect(10.0:1.0:20.0)) # OKAY, using bracketing, able to evaluate for an array of n

    @show gen_V_temp(10.0,80.0,b,params, vFuncs, regime) # OKAY 
    gen_EV!(params,vFuncs,iterObj_i,regime) # OKAY
end

divv = n .+ (1 - params.cL) * params.β * δ .+ QB .* B .- L .- params.g * δ .- SS .- params.cM .* L.^2 .- params.cO
G = psi.([params], divv) .+ params.β * iterObj_i.EV

sol = solve_bank_problem(params, vFuncs, iDelta, iLambda, iN, G)

#################################### Dive into solve_bank_problem

G_itp = interpolate(G, BSpline(Quadratic(Line(OnGrid())))) 
lRange = range(first(params.lGrid), last(params.lGrid), length(params.lGrid))
sRange = range(first(params.sGrid), last(params.sGrid), length(params.sGrid))
bRange = range(first(params.bGrid), last(params.bGrid), length(params.bGrid))
G_itp = Interpolations.scale(G_itp, lRange, sRange, bRange)
G_itp_ext = extrapolate(G_itp, Line())     # enables linear extrapolation
G_interp(l, s, b) = G_itp_ext(l, s, b)

G_interp(1.2, -2.1, 10.005)

(l + params.g*params.deltaGrid[iDelta] +200.0 - params.deltaGrid[iDelta]-b)/(params.wr*l) >= params.α
(10.0 + params.g*params.deltaGrid[iDelta]) <= params.β*params.deltaGrid[iDelta]






























vFuncs = Initiate_vFunc(params)
iterObj_i = Initiate_IterObj_i(params)
QB = @view vFuncs.qBond[:, :, :, iDelta, iLambda] # nL × nS × nB

L = reshape(params.lGrid, :, 1, 1)    # nL × 1 × 1       
SS = reshape(params.sGrid, 1, :, 1)   # 1  × nS × 1
B = reshape(params.bGrid, 1, 1, :)    # 1  × 1  × nB

div = n .+ (1 - params.cL) * params.β * δ .+ QB .* B .- L .- params.g * δ .- SS .- params.cM .* L.^2 .- params.cO
        
# div(l,s,b) = n + (1-params.cL)*params.β*δ + vFuncs.qBond[l,s,b,iDelta,iLambda]*b - l - params.g * δ -s-params.cM*l^2 -params.cO 
# div = n .+ (1-params.cL)*params.β*δ .+ vFuncs.qBond[:, :, :, iDelta, iLambda]*b - L .- params.g*δ .- S - params.cM*L.^2 .- params.cO # (nL, nS, nB) matrix
gen_EV!(params, vFuncs, iterObj_i, regime)
G = psi.([params], div) .+ params.β * iterObj_i.EV # G(l,s,b,n) = flow utility[l,s,b,n] + β * EV[l,s,b], (nL, nS, nB) matrix with fixed n 



# G is well defined! 
G_itp = interpolate(G, BSpline(Quadratic(Line(OnGrid())))) 
lRange = range(first(params.lGrid), last(params.lGrid), length(params.lGrid))
sRange = range(first(params.sGrid), last(params.sGrid), length(params.sGrid))
bRange = range(first(params.bGrid), last(params.bGrid), length(params.bGrid))
G_itp = Interpolations.scale(G_itp, lRange, sRange, bRange)
G_itp_ext = extrapolate(G_itp, Line())     # enables linear extrapolation
G_interp(l, s, b) = G_itp_ext(l, s, b)


model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "tol", 1e-8)
set_optimizer_attribute(model, "max_iter", 10)
set_optimizer_attribute(model, "acceptable_tol", 1e-6)
JuMP.register(model, :G_interp, 3, G_interp; autodiff = true)
l_min, l_max = first(params.lGrid) + eps(), last(params.lGrid)
s_min, s_max = first(params.sGrid), last(params.sGrid)
b_min, b_max = first(params.bGrid), last(params.bGrid)
@variable(model, l_min <= l <= l_max, start = (l_min+l_max)/2)
@variable(model, s_min <= s <= s_max, start = (s_min+s_max)/2)
@variable(model, b_min <= b <= b_max, start = (b_min+b_max)/2)

@NLobjective(model, Max, G_interp(l,s,b)) # maximizing interpolated G(l,s,b)
@NLconstraint(model, (l + params.g*params.deltaGrid[iDelta] +s - params.deltaGrid[iDelta]-b)/(params.wr*l) >= params.α) # constratint 1
@NLconstraint(model, (l + params.g*params.deltaGrid[iDelta]) <= params.β*params.deltaGrid[iDelta]) # constratint 2

optimize!(model)




stat = termination_status(model)

if stat == MathOptInterface.OPTIMAL || stat == MathOptInterface.LOCALLY_SOLVED
    [value(l), value(s), value(b)]
else
   error("Optimization did not converge: status = $stat")
end