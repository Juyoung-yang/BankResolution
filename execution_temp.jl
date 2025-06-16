using JuMP
using Ipopt
using Interpolations
using MathOptInterface

include("C:\\Users\\master\\Desktop\\(2025) banking regulation\\BankResolution\\parameters.jl")
# include("C:\\Users\\master\\Desktop\\(2025) banking regulation\\BankResolution\\main.jl")

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
    Rl::T # loan interest rate
    X::Array{T,4} # bank's failure decision 
end

struct IterObj_i{T<:Real,S<:Integer} # objects given state (iDelta, iLambda), objects used in VFI_i
    EV::Array{T,3} # EV[l,s,b]
    G::Array{T,4} # G(l,s,b,n)
    solution::Array{T,1} # optimizer = [(l,s,b)]
    solution_index::Array{T,1} # optimizer = [(il,is,ib)]
    failure::Array{T,1} # failure decision = [(fail or not fail)]
end

function Initiate_Params(qd::T,β::T,Rf::T,wr::T,α::T,ρ::T,g::T,ξ::T,cF::T,dBar::T,σ::T,τC::T,z::T,α1::T,α2::T,α3::T,δL::T,δM::T,δH::T,cM::T,cO::T,cL::T,H::Array{T,2},Γ::Array{T,2},λL::T,λM::T,λH::T,γ::T,ϕ::T,n_start::T,n_npts::S,n_stop::T,l_start::T,l_npts::S,l_stop::T,s_start::T,s_npts::S,s_stop::T,b_start::T,b_npts::S,b_stop::T) where {T<:Real,S<:Integer}
    
    deltaGrid = [δL,δM,δH] # Define a Tuple, immutable 
    lambdaGrid = [λL,λM,λH] # regular array, mutable
    nGrid = range(n_start,stop=n_stop,length=n_npts)

    lGrid = range(l_start,stop=l_stop,length=l_npts)
    sGrid = range(s_start,stop=s_stop,length=s_npts)
    bGrid = range(b_start,stop=b_stop,length=b_npts)

    pam = Params{T,S}(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,α1,α2,α3,δL,δM,δH,cM,cO,cL,H,Γ,λL,λM,λH,γ,ϕ,deltaGrid,lambdaGrid,nGrid,lGrid,sGrid,bGrid)
    return pam
end

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
    solution_index = zeros(T, length(params.nGrid))
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

function solve_bank_problem2(params::Params{T,S},vFuncs::VFuncs{T,S},iDelta::S,iLambda::S,iN::S,G::Array{T,3}) where {T<:Real, S<:Integer}
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
        set_optimizer_attribute(model, "max_iter", 25000)
        set_optimizer_attribute(model, "print_level", 5)
        JuMP.register(model, :G_interp, 3, G_interp; autodiff = true)
        l_min, l_max = (last(params.lGrid)+first(params.lGrid))/2, last(params.lGrid)
        @show l_min, l_max
        s_min, s_max = first(params.sGrid), last(params.sGrid)
        b_min, b_max = first(params.bGrid), last(params.bGrid)
        l_start = l_max
        s_start = (s_max+s_min)/2
        b_start = (b_max+b_min)/2
        @variable(model, l_min <= l <= l_max, start = l_start)
        @variable(model, s_min <= s <= s_max, start = s_start)
        @variable(model, b_min <= b <= b_max, start = b_start)

        @NLobjective(model, Max, G_interp(l,s,b)) # maximizing interpolated G(l,s,b)
        @constraint(model, (1-params.α*params.wr)*l + s - b >= (1-params.g)*params.deltaGrid[iDelta])# constratint 1
        @constraint(model, l <= (params.β - params.g)*params.deltaGrid[iDelta]) # constratint 2

        optimize!(model)
        @show is_solved_and_feasible(model)
        @show stat = termination_status(model)

        if stat == MathOptInterface.OPTIMAL || stat == MathOptInterface.LOCALLY_SOLVED
            return [value(l), value(s), value(b)]
        else
            error("Optimization did not converge: status = $stat")
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
@show sum(G)
sol = solve_bank_problem2(params, vFuncs, iDelta, iLambda, iN, G)

#################################### Dive into solve_bank_problem

G_itp = interpolate(G, BSpline(Quadratic(Line(OnGrid())))) 
lRange = range(first(params.lGrid), last(params.lGrid), length(params.lGrid))
sRange = range(first(params.sGrid), last(params.sGrid), length(params.sGrid))
bRange = range(first(params.bGrid), last(params.bGrid), length(params.bGrid))
G_itp = Interpolations.scale(G_itp, lRange, sRange, bRange)
G_itp_ext = extrapolate(G_itp, Line())     # enables linear extrapolation
G_interp(l, s, b) = G_itp_ext(l, s, b)

G_interp(1.2, -2.1, 10.005)































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