pwd()
readdir()

include("C:\\Users\\master\\Desktop\\(2025) banking regulation\\BankResolution\\parameters.jl");
include("C:\\Users\\master\\Desktop\\(2025) banking regulation\\BankResolution\\main_v2.jl");

# 0. set the parameters
params = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,α1,α2,α3,δL,δM,δH,cM,cO,cL,ϵ,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)

# 1. solve the model under given params
@time sol = solve_model(params,true, 0.01, 1.9);
@time policy_ex = Get_PolicyFuncs(params, sol.eq.Iterobj_is, sol.Rl_star); # get policy functions from the solution

# 2. generate shock sequenc, simulate paths, and calculate moments
@time shocks = makeShockSequence(params, sol.eq.vFuncs, 100, 100, 10); # make shock sequence for simulation
paths = Initiate_Paths(shocks, params); # initiate paths for simulation

@time paths_ex = Simulate_paths(paths, policy_ex, sol.eq.vFuncs, params, sol.Rl_star, 10); # simulate paths and calculate moments
################################################################################################################################################

plot(paths_ex.failureSim[:,1,1]) # loan policy function for the first bank J at the first observation N
plot!(paths_ex.failureSim[:,1,3])
plot(paths_ex.lSim[:,1,1])
plot!(paths_ex.lSim[:,1,2])
paths_ex.failureSim[:,1,1]
plot(paths_ex.bSim[:,1,1])
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

struct PolicyFuncs{T<:Real,S<:Integer}
    lPolicy::Array{T,3} # loan policy function (delta, lambda, n)
    sPolicy::Array{T,3} # savings policy function
    bPolicy::Array{T,3} # bond policy function
    failure::Array{Bool,4} # failure or not (delta, lambda, n, lambdaPrime)
end

struct SimShocks{S<:Integer}
    deltaIndSim::Array{S,3}
    lambdaIndSim::Array{S,3}
    nInitialIndSim::Array{S,2} # initial net asset value (index) for each bank J
end

struct SimPaths{T<:Real,S<:Integer}
    deltaIndSim::Array{S,3}
    lambdaIndSim::Array{S,3}
    nInitialIndSim::Array{S,2} # initial net asset value (index) for each bank J
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
    divSim::Array{T,3} # dividend (배당금); current period
    nIndSim::Array{S,3} # net asset value (index), current period 
    nPrimeSim::Array{T,3} # net asset value, next period with lambdaPrime
    nPrimeIndSim::Array{S,3} # net asset value (index), next period with lambdaPrime
    assetSim::Array{T,3} # asset, current period
    assetPrimeSim::Array{T,3} # asset, next period with lambdaPrime
    Γ::Array{T,3} # distribution of states, current period
    ΓPrime::Array{T,3} # distribution of states, next period with lambdaPrime
end

################################################################################################################################################

# draw bank types from Γ and simulate (delta, lambda) from Markov Chain 
function makeShockSequence(params::Params{T,S}, vFuncs::VFuncs, bigT::S, bigN::S, bigJ::S) where {T<:Real,S<:Integer}

    deltaIndSim = zeros(Int, bigT, bigN, bigJ); # delta states for each bank J
    lambdaIndSim = zeros(Int, bigT, bigN, bigJ); # lambda states for each bank J
    nInitialIndSim = zeros(Int, bigN, bigJ); # initial net asset value (index) for each bank J

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
            deltaIndSim[1, smallN, smallJ] = idx[number][1]; 
            lambdaIndSim[1, smallN, smallJ] = idx[number][2]; 
            nInitialIndSim[smallN, smallJ] = idx[number][3]; 
        end
    end

    # 3. for each period T, draw (delta, lambda) from the markov process 
    H_Markov = [Categorical(params.H[i, :]) for i in 1:size(H,1)]
    F_Markov = [Categorical(params.F[i, :]) for i in 1:size(F,1)]

    for smallN = 1:bigN
        for smallJ = 1:bigJ # for bank J
            for smallT = 2:bigT # for each period T
                @show rand(H_Markov[Int(deltaIndSim[smallT-1, smallN, smallJ])])
                deltaIndSim[smallT, smallN, smallJ] = rand(H_Markov[Int(deltaIndSim[smallT-1, smallN, smallJ])]) # draw delta from the Markov chain
                lambdaIndSim[smallT, smallN, smallJ] = rand(F_Markov[Int(lambdaIndSim[smallT-1, smallN, smallJ])]) # draw lambda from the Markov chain                
            end
        end
    end

   return SimShocks{S}(deltaIndSim, lambdaIndSim, nInitialIndSim); # store the simulated shocks in SimShocks struct
end

function Initiate_PolicyFuncs(params::Params{T,S}) where {T<:Real,S<:Integer}

    lPolicy = zeros(T, length(params.deltaGrid), length(params.lambdaGrid), length(params.nGrid)); # loan policy function (delta, lambda, n)
    sPolicy = zeros(T, length(params.deltaGrid), length(params.lambdaGrid), length(params.nGrid)); # savings policy function
    bPolicy = zeros(T, length(params.deltaGrid), length(params.lambdaGrid), length(params.nGrid)); # bond policy function
    failure = falses(length(params.deltaGrid), length(params.lambdaGrid), length(params.nGrid), length(params.lambdaGrid)); # failure or not (delta, lambda, n)

    policyy = PolicyFuncs{T,S}(lPolicy, sPolicy, bPolicy, failure)
    return policyy;
end

function Initiate_Paths(shocks::SimShocks{S}, params::Params{T,S}) where {T<:Real,S<:Integer}

    deltaIndSim = shocks.deltaIndSim;
    lambdaIndSim = shocks.lambdaIndSim;
    nInitialIndSim = shocks.nInitialIndSim;

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
    nSim = zeros(T, bigT, bigN, bigJ); # retained earining (잉여금; 기업의 이익 중에서 배당금으로 배분하지 않고 재투자를 위해 남겨둔 금액); current period
    divSim = zeros(T, bigT, bigN, bigJ); # dividend (배당금); current period
    # 
    nPrimeSim = zeros(T, bigT, bigN, bigJ); # retained earining (잉여금; 기업의 이익 중에서 배당금으로 배분하지 않고 재투자를 위해 남겨둔 금액); next period with lambdaPrime
    # 
    assetPrimeSim = zeros(T, bigT, bigN, bigJ); # total asset (총자산); next period with lambdaPrime

    Γ = zeros(3, 3, length(params.nGrid)); # stationary distribution
    ΓPrime = zeros(3, 3, length(params.nGrid)); # next period stationary distribution

    nIndSim = zeros(S, bigT, bigN, bigJ);
    nPrimeIndSim = zeros(S, bigT, bigN, bigJ);
    assetSim = zeros(T, bigT, bigN, bigJ); # total asset (총자산); current period

    return SimPaths(deltaIndSim, lambdaIndSim, nInitialIndSim, bigT, bigN, bigJ, deltaSim, lambdaSim, lSim, sSim, bSim, failureSim, nSim, divSim, nIndSim, nPrimeSim, nPrimeIndSim, assetSim, assetPrimeSim, Γ, ΓPrime)
end

################################################################################################################################################

function Get_PolicyFuncs(params::Params{T,S}, Iterobj_is::Matrix{IterObj_i{T,S}}, Rl::T) where {T<:Real,S<:Integer,F<:Bool}

    NAV(l::T,s::T,b::T,lambda::T,δ::T)::T = Rl*( (1-lambda)*l + params.g*δ) + (1+params.Rf)*s - δ - b; # 순자산, net asset value 

    # 1. initiate policy functions
    policy = Initiate_PolicyFuncs(params);
    @show policy.failure
    @show policy.failure[1,1,1,:] # check the failure or not for the first state

    # 3. calculate the policy functions (l,s,b) and failure or not
    for iDelta in 1:length(params.deltaGrid)
        for iLambda in 1:length(params.lambdaGrid)
            iterobj = Iterobj_is[iDelta, iLambda]; 
            for iN in 1:length(params.nGrid)
                policy.lPolicy[iDelta, iLambda, iN] = iterobj.solution[iN][1]; # loan policy function
                policy.sPolicy[iDelta, iLambda, iN] = iterobj.solution[iN][2]; # savings policy function
                policy.bPolicy[iDelta, iLambda, iN] = iterobj.solution[iN][3]; # bond policy function
                policy.failure[iDelta, iLambda, iN, :] .= [NAV(iterobj.solution[iN][1], iterobj.solution[iN][2], iterobj.solution[iN][3], params.lambdaGrid[iLambdaPrime], params.deltaGrid[iDelta]) < 0 for iLambdaPrime in 1:length(params.lambdaGrid)]; # failure or not, denoted as 1 or 0
            end
        end
    end

    return policy;
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

    return paths
    # 2. after obtaining simulation paths, trim some first parts and calculate moments
   #  moments = calculate_moments(paths, policy, vFuncs, params, Rl, trim);

    # 3. return key moments of interest
   # return moments 
end

# for bank smallJ, calculate the variable of interest at time (smallT, smallN); given state (delta, lambda, n), calculate the series of interest
function calculate_series(paths::SimPaths{T,S}, policy::PolicyFuncs{T,S}, vFuncs::VFuncs{T,S}, params::Params{T,S}, Rl::T, smallT::S, smallN::S, smallJ::S) where {T<:Real,S<:Integer}

    NAV(l::T,s::T,b::T,lambda::T,δ::T)::T = Rl*( (1-lambda)*l + params.g*δ) + (1+params.Rf)*s - δ - b; # 순자산, net asset value 
    totalAsset(l::T,s::T,lambda::T,δ::T)::T = Rl * ( (1-lambda)*l + params.g*δ) + (1+params.Rf)*s; # total asset of the bank, 총자산
    tax(l::T,s::T,b::T,lambda::T,δ::T)::T = params.τC * max(0, (Rl-1)*((1-lambda)*l +params.g*δ) + params.Rf*s - params.Rf *(δ + b)); # tax on the bank's asset value
    n_success(l::T,s::T,b::T,lambda::T,δ::T)::T = NAV(l,s,b,lambda,δ) - tax(l,s,b,lambda,δ); # next period asset conditional on bank success 
    n_failure(l::T, lambda::T)::T = params.α * params.wr * Rl * (1-lambda)*l; # next period asset conditional on bank failure 

    # Interpolate qBond for a given value of choice (l,s,b) and state (delta, lambda)
    function qBond_interpolated(l::T, s::T, b::T, deltaInd::S, lambdaInd::S, params::Params{T,S}, vFuncs::VFuncs{T,S})::T

        lRange = range(first(params.lGrid), last(params.lGrid), length(params.lGrid));
        sRange = range(first(params.sGrid), last(params.sGrid), length(params.sGrid));
        bRange = range(first(params.bGrid), last(params.bGrid), length(params.bGrid));
        

        q = vFuncs.qBond[:, :, :, deltaInd, lambdaInd];
        qmin = minimum(q);
        qmax = maximum(q);
        qnorm = (q .- qmin) ./ (qmax - qmin);
        q_itp_rev = interpolate(qnorm, BSpline(Quadratic(Line(OnGrid()))));
        q_itp_rev = Interpolations.scale(q_itp_rev, lRange, sRange, bRange);
        q_itp_ext_rev = extrapolate(q_itp_rev, Line());
        q_interp_rev(ll, ss, bb) = q_itp_ext_rev(ll, ss, bb);

        q_value = q_interp_rev(l, s, b) * (qmax - qmin) + qmin; # interpolate qBond for a given choice (l,s,b)
        @show q_value, q_interp_rev(l, s, b), (l,s,b)
        return q_value
    end

    # Interpolate (l,s,b) for a given n value and a given state (delta, lambda)
    function lsb_interpolated(params::Params{T,S}, policy::PolicyFuncs{T,S}, deltaInd::S, lambdaInd::S, n::T)::NTuple{3,T}
        # Interpolate l, s, b for a given n value and a given state (delta, lambda)

        nRange = range(first(params.nGrid), last(params.nGrid), length(params.nGrid));

        l_array = policy.lPolicy[deltaInd, lambdaInd, :]; # loan policy function for a given state (delta, lambda)
        l_array_min, l_array_max = minimum(l_array), maximum(l_array);
        l_array_norm = (l_array .- l_array_min) ./ (l_array_max - l_array_min);
        l_itp = interpolate(l_array_norm, BSpline(Quadratic(Line(OnGrid()))));
        l_itp = Interpolations.scale(l_itp, nRange);    
        l_itp_ext = extrapolate(l_itp, Line());
        l_interp(n_val) = l_itp_ext(n_val);
        l_value = l_interp(n) * (l_array_max - l_array_min) + l_array_min; # interpolate loan policy function for a given n value


        s_array = policy.sPolicy[deltaInd, lambdaInd, :]; # savings policy function for a given state (delta, lambda)
        s_array_min, s_array_max = minimum(s_array), maximum(s_array);
        s_array_norm = (s_array .- s_array_min) ./ (s_array_max - s_array_min);
        s_itp = interpolate(s_array_norm, BSpline(Quadratic(Line(OnGrid()))));
        s_itp = Interpolations.scale(s_itp, nRange);
        s_itp_ext = extrapolate(s_itp, Line());
        s_interp(n_val) = s_itp_ext(n_val);
        s_value = s_interp(n) * (s_array_max - s_array_min) + s_array_min; # interpolate savings policy function for a given n value

        b_array = policy.bPolicy[deltaInd, lambdaInd, :]; # bond policy function for a given state (delta, lambda)
        b_array_min, b_array_max = minimum(b_array), maximum(b_array);
        b_array_norm = (b_array .- b_array_min) ./ (b_array_max - b_array_min);
        b_itp = interpolate(b_array_norm, BSpline(Quadratic(Line(OnGrid()))));
        b_itp = Interpolations.scale(b_itp, nRange);
        b_itp_ext = extrapolate(b_itp, Line());
        b_interp(n_val) = b_itp_ext(n_val);
        b_value = b_interp(n) * (b_array_max - b_array_min) + b_array_min; # interpolate bond policy function for a given n value

        return (l_value,s_value,b_value)
    end

    # 1. at (smallT, smallN, smallJ), calculate (delta, lambda, lambdaPrime)
    delta = params.deltaGrid[paths.deltaIndSim[smallT, smallN, smallJ]]; # delta state for bank J at time T
    lambda = params.lambdaGrid[paths.lambdaIndSim[smallT, smallN, smallJ]]; # lambda state for bank J at time T
    lambdaPrime = params.lambdaGrid[paths.lambdaIndSim[smallT+1, smallN, smallJ]];

    paths.deltaSim[smallT, smallN, smallJ] = delta;
    paths.lambdaSim[smallT, smallN, smallJ] = lambda;

    # 2. at (smallT, smallN, smallJ), calculate choice variables (l,s,b) and dividend (div), asset, failure, next period asset (nPrime)
    if smallT == 1 # initial period 
        paths.nSim[smallT, smallN, smallJ] = params.nGrid[paths.nInitialIndSim[smallN, smallJ]]; ##### start with first grid n point
        paths.lSim[smallT, smallN, smallJ] = policy.lPolicy[paths.deltaIndSim[smallT, smallN, smallJ], paths.lambdaIndSim[smallT, smallN, smallJ], paths.nInitialIndSim[smallN, smallJ]]; # loan policy function
        paths.sSim[smallT, smallN, smallJ] = policy.sPolicy[paths.deltaIndSim[smallT, smallN, smallJ], paths.lambdaIndSim[smallT, smallN, smallJ], paths.nInitialIndSim[smallN, smallJ]]; # savings policy function
        paths.bSim[smallT, smallN, smallJ] = policy.bPolicy[paths.deltaIndSim[smallT, smallN, smallJ], paths.lambdaIndSim[smallT, smallN, smallJ], paths.nInitialIndSim[smallN, smallJ]]; # bond policy function
        paths.divSim[smallT, smallN, smallJ] = params.nGrid[1] + (1-params.cL) * params.β * delta + qBond_interpolated(paths.lSim[smallT, smallN, smallJ], paths.sSim[smallT, smallN, smallJ], paths.bSim[smallT, smallN, smallJ], paths.deltaIndSim[smallT, smallN, smallJ], paths.lambdaIndSim[smallT, smallN, smallJ], params, vFuncs) * paths.bSim[smallT, smallN, smallJ] - paths.lSim[smallT, smallN, smallJ] - params.g * delta - paths.sSim[smallT, smallN, smallJ] - params.cM * paths.lSim[smallT, smallN, smallJ]^2 - params.cO; # dividend (배당금); current period
        paths.assetPrimeSim[smallT, smallN, smallJ] = totalAsset(paths.lSim[smallT, smallN, smallJ], paths.sSim[smallT, smallN, smallJ], lambdaPrime, delta);  # total asset given lambdaPrime
        paths.failureSim[smallT, smallN, smallJ] = 0; # no failure in the first period (i.e. period zero), denoted as 0
        paths.failureSim[smallT+1, smallN, smallJ] = ( paths.assetSim[smallT, smallN, smallJ] - paths.bSim[smallT, smallN, smallJ] - delta < 0 ) ? 1 : 0; # failure or not, denoted as 1 or 0

        if paths.failureSim[smallT+1, smallN, smallJ] == 1 # if failure
            paths.nPrimeSim[smallT, smallN, smallJ] = n_failure(paths.lSim[smallT, smallN, smallJ], lambdaPrime); # next period asset conditional on bank failure 
        else # if success
            paths.nPrimeSim[smallT, smallN, smallJ] = n_success(paths.lSim[smallT, smallN, smallJ], paths.sSim[smallT, smallN, smallJ], paths.bSim[smallT, smallN, smallJ], lambdaPrime, delta); # next period asset conditional on bank success 
        end

    else
        paths.nSim[smallT, smallN, smallJ] = paths.nPrimeSim[smallT-1, smallN, smallJ] # 잉여금 as a embodied state 
        paths.lSim[smallT, smallN, smallJ] = lsb_interpolated(params, policy, paths.deltaIndSim[smallT, smallN, smallJ], paths.lambdaIndSim[smallT, smallN, smallJ], paths.nSim[smallT, smallN, smallJ])[1]; # loan policy function
        paths.sSim[smallT, smallN, smallJ] = lsb_interpolated(params, policy, paths.deltaIndSim[smallT, smallN, smallJ], paths.lambdaIndSim[smallT, smallN, smallJ], paths.nSim[smallT, smallN, smallJ])[2]; # savings policy function
        paths.bSim[smallT, smallN, smallJ] = lsb_interpolated(params, policy, paths.deltaIndSim[smallT, smallN, smallJ], paths.lambdaIndSim[smallT, smallN, smallJ], paths.nSim[smallT, smallN, smallJ])[3]; # bond policy function
        paths.divSim[smallT, smallN, smallJ] = paths.nSim[smallT, smallN, smallJ] + (1-params.cL) * params.β * delta + qBond_interpolated(paths.lSim[smallT, smallN, smallJ], paths.sSim[smallT, smallN, smallJ], paths.bSim[smallT, smallN, smallJ], paths.deltaIndSim[smallT, smallN, smallJ], paths.lambdaIndSim[smallT, smallN, smallJ], params, vFuncs) * paths.bSim[smallT, smallN, smallJ] - paths.lSim[smallT, smallN, smallJ] - params.g * delta - paths.sSim[smallT, smallN, smallJ] - params.cM * paths.lSim[smallT, smallN, smallJ]^2 - params.cO; # dividend (배당금); current period
        paths.assetPrimeSim[smallT, smallN, smallJ] = totalAsset(paths.lSim[smallT, smallN, smallJ], paths.sSim[smallT, smallN, smallJ], lambdaPrime, delta);  # total asset given lambdaPrime
        paths.failureSim[smallT+1, smallN, smallJ] = ( paths.assetSim[smallT, smallN, smallJ] - paths.bSim[smallT, smallN, smallJ] - delta < 0 ) ? 1 : 0; # failure or not, denoted as 1 or 0

        if paths.failureSim[smallT+1, smallN, smallJ] == 1 # if failure
            paths.nPrimeSim[smallT, smallN, smallJ] = n_failure(paths.lSim[smallT, smallN, smallJ], lambdaPrime); # next period asset conditional on bank failure 
        else # if success
            paths.nPrimeSim[smallT, smallN, smallJ] = n_success(paths.lSim[smallT, smallN, smallJ], paths.sSim[smallT, smallN, smallJ], paths.bSim[smallT, smallN, smallJ], lambdaPrime, delta); # next period asset conditional on bank success 
        end
    end

    # return paths; # return the updated paths

    # 3. at (smallT, smallN, smallJ), calculate Γ and ΓPrime
   # if smallT == 1 
    #    paths.Γ[smallT, smallN, smallJ] = 
    #    paths.ΓPrime[smallT, smallN, smallJ] = 
   # else
   #     paths.Γ[smallT, smallN, smallJ] = paths.ΓPrime[smallT-1, smallN, smallJ]
   #     paths.ΓPrime[smallT, smallN, smallJ] = 
   # end
end