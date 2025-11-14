# main execution code that can reduce running time: 

function solve_simulate_and_moments_reuse(params::Params{T,S},params_cal::Params_cal{T,S},regime::F) where {T<:Real,S<:Integer,F<:Bool}

    println("<<<<<<Solving the model...>>>>>>")

    # --- 1️⃣ Preallocate once for reuse through root-finding and final solve ---
    vFuncs       = Initiate_vFunc(params)
    vFuncsNew    = Initiate_vFuncNew(params)
    Iterobj_is   = Initiate_MatrixIterObj_i(params)

    # --- 2️⃣ Define root function that reuses the same objects ---
    solve_model_given_r_single = Rl -> solve_model_given_r!(Rl;Params = params,Regime = regime,vFuncs = vFuncs,vFuncsNew = vFuncsNew,Iterobj_is = Iterobj_is)

    # --- 3️⃣ Root-finding for equilibrium Rl_star ---
    Rl_safe = clamp(params_cal.a, 1.0 + eps(), 2.0 - eps())
    Rl_star = find_zero(solve_model_given_r_single,Rl_safe,method = Roots.Secant();tol = 1e-2,maxevals = 100)

    println("<<<<<<Solved the model. Equilibrium Rl =" , Rl_star, ">>>>>>>")

    # --- 4️⃣ Final in-place equilibrium solve ---
    VFI!(params, vFuncs, vFuncsNew, Iterobj_is, Rl_star, regime, 1000, 1e-4)
    stationary_distribution!(params, Rl_star, vFuncs, vFuncsNew, Iterobj_is, 1000, 1e-4)

    # --- 5️⃣ Continue to policy and simulation using the already solved state ---
    policy = Get_PolicyFuncs(params, Iterobj_is, vFuncs, Rl_star, regime)
    println("<<<<<<Obtained policy functions>>>>>>")

    simulation = simulate_and_moments(params, params_cal, vFuncs, policy, Rl_star, regime)
    @show simulation.moments.moments

    return (moments = simulation.moments,
            Rl_star = Rl_star,
            vFuncs = vFuncs,
            vFuncsNew = vFuncsNew,
            Iterobj_is = Iterobj_is,
            policy = policy,
            paths = simulation.paths)
end

function solve_model_given_r!(Rl::T;Params::Params{T,S},Regime::F,vFuncs::VFuncs{T,S},vFuncsNew::VFuncsNew{T,S},Iterobj_is::Matrix{IterObj_i{T,S}}) where {T<:Real,S<:Integer,F<:Bool}

    println("(SUB BLOCK) -- VFI! (reuse enabled): Regime = ", Regime)
    VFI!(Params, vFuncs, vFuncsNew, Iterobj_is, Rl, Regime, 500, 1e-2)

    println("(SUB BLOCK) -- stationary_distribution! (reuse enabled)")
    stationary_distribution!(Params, Rl, vFuncs, vFuncsNew, Iterobj_is, 500, 1e-2)

    excess_loan_supply = aggre_loan_supply(Params, vFuncs, Iterobj_is) - Rl^(Params.ϵ) * Params.E
    return excess_loan_supply
end

function VFI!(params::Params{T,S},vFuncs::VFuncs{T,S},vFuncsNew::VFuncsNew{T,S},Iterobj_is::Matrix{IterObj_i{T,S}},Rl::T,regime::F, maxiter::S, tol::T) where {T<:Real, S<:Integer, F<:Bool}

        iter, maxdiffs, diffs = 1, one(T), zeros(T);

        βF = nothing
        fail_temp_buffer = nothing
        q_temp_buffer = nothing

        if regime == true
            βF = params.β .* params.F                     # 1회 계산
           #  nthreads = Threads.nthreads()
           # fail_temp_buffer = [Vector{T}(undef, 3) for _ in 1:nthreads]
          #  @show size(fail_temp_buffer)
          #  q_temp_buffer = [Vector{T}(undef, 3) for _ in 1:nthreads]
          fail_temp_buffer = [Vector(undef, 3) for _ in 1:5]
          q_temp_buffer = [Vector(undef, 3) for _ in 1:5]
        end

        println("VFI 시작 시 VF 평균값: ", mean(vFuncs.VF))

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
                qBond_specialRegime!(params,vFuncsNew,Rl,βF,fail_temp_buffer,q_temp_buffer) # set qBond to be 1 for all states 
            end

            # 감쇠(damping) 업데이트로 진동 억제
            @. vFuncs.VF = 0.7 * vFuncs.VF + 0.3 * vFuncsNew.VF
            @. vFuncs.qBond = 0.7 * vFuncs.qBond + 0.3 * vFuncsNew.qBond
            @. vFuncs.X = 0.7 * vFuncs.X + 0.3 * vFuncsNew.X

            # 수렴도 계산
            maxdiffs = maximum(abs.(vFuncs.VF .- vFuncsNew.VF))

            if mod(iter, 100) == 0
                println("VFI!: iter=$iter, maxdiffs=$maxdiffs")
            end

            iter += 1
        end

    println("VFI! 종료 시 VF 평균값: ", mean(vFuncs.VF))
    return nothing
end

# debt pricing inqBond_condiState! counterfactural regime
function qBond_condiState!(params::Params{T,S}, Rl::T, il::S ,is::S, ib::S, iDelta::S,fail_temp::Vector,q_temp::Vector,βF::AbstractMatrix{T}) where {T<:Real,S<:Integer} 
        
        l, s, b, Delta = params.lGrid[il], params.sGrid[is], params.bGrid[ib], params.deltaGrid[iDelta]
        # lambdaStar(l,s,b,Delta) = (Rl*(l+params.g*Delta)+(1+params.Rf)*s-Delta-b)/(Rl*l)
        lambdaStar_capital(l,s,b,Delta) = 1- (-s + b + (1-params.g)*Delta)/(1-params.α * params.wr)/l
        @show λStar = clamp( lambdaStar_capital(l,s,b,Delta), 0.0, 1.0) ; # lambda cutoff λStar in [0,1]
       
        
        # --- 3️⃣ 실패시 수익률 계산 (fail_temp에 in-place 저장) ---
        if b == eps()
            fill!(fail_temp, 0.0) # 대출이 없으면 fail과 상관없이 loan profit이 없음 
        else
            @inbounds for j in 1:3
                Λ = params.lambdaGrid[j]
                NetAsset = Rl * ((one(T) - Λ) * l + params.g * Delta) + (one(T) + params.Rf*0.8) * s - Delta 
                hi   = max(zero(T), NetAsset)
                VLB = min(b, hi)
                bHat = NetAsset - params.α * params.wr * Rl * (one(T) - Λ) * l
                num = max(VLB, bHat)
                fail_temp[j] = num / b
            end
        end
        @show fail_temp

        # --- 4️⃣ λ* 위치에 따라 q_temp 구성 (Nλ=3 기준), 기본적으로 lambda가 커야 은행 실패 ---
        @inbounds begin
            if λStar < params.λL # fail all the time  
               @. q_temp = fail_temp
            elseif λStar < params.λM
                q_temp[1] = one(T) # no fail when lambda is the lowest
                q_temp[2] = fail_temp[2]
                q_temp[3] = fail_temp[3]
            elseif λStar < params.λH
                q_temp[1] = one(T)
                q_temp[2] = one(T)
                q_temp[3] = fail_temp[3]
            else # failure at all time
                fill!(q_temp, one(T))
            end
        end

        # --- 5️⃣ βF * q_temp  (in-place 덮어쓰기) ---
        rhs = copy(q_temp)
        mul!(q_temp, βF, rhs)
end

function qBond_specialRegime!(params::Params{T,S},vFuncsNew::VFuncsNew{T,S},Rl::T,βF::AbstractMatrix{T},fail_temp_buffer::Vector,q_temp_buffer::Vector) where {T<:Real,S<:Integer}

  Threads.@threads for il in 1:length(params.lGrid)
   #  for il in length(params.lGrid)
           #  tid = Threads.threadid()
             fail_temp = fail_temp_buffer[il]
             q_temp = q_temp_buffer[il]

            for is in 1:length(params.sGrid)
                for ib in 1:length(params.bGrid)
                    for iDelta in 1:length(params.deltaGrid)
                        qBond_condiState!(params,Rl,il,is,ib,iDelta,fail_temp,q_temp,βF)
                        @. vFuncsNew.qBond[il,is,ib,iDelta,:] = 0.7 * vFuncsNew.qBond[il,is,ib,iDelta,:] + 0.3 * q_temp
                    end
                end
            end
        end
end


function stationary_distribution!(params::Params{T,S},Rl::T,vFuncs::VFuncs{T,S},vFuncsNew::VFuncsNew{T,S},Iterobj_is::Matrix{IterObj_i{T,S}},maxiter::S,tol::T) where {T<:Real,S<:Integer}

    iter, maxdiffs, diffs = 1, one(T), zeros(T);

    # --- 상태공간 총 조합 수 계산 (flat loop용) ---
    nΔ  = length(params.deltaGrid)
    nΛ  = length(params.lambdaGrid)
    nN  = length(params.nGrid)

    if iter == 1
        println("Γ 합계 초기값: ", sum(vFuncs.Γ))
    end

    while iter ≤ maxiter && maxdiffs > tol

        # 먼저 vFuncsNew.Γ를 초기화
        fill!(vFuncsNew.Γ, 0.0)

        # (1) 한 층만 thread 병렬화
        Threads.@threads for iΛp in 1:nΛ 
      #   for iΔ in 1:nΔ
            for iΛ in 1:nΛ, iN in 1:nN, iΔ in 1:nΔ
                Update_stationary_dist(params, vFuncs, vFuncsNew,Iterobj_is, iΔ, iΛ, iN, iΛp, Rl)
            end
        end


        # (2) 감쇠(damping) 업데이트
        damping = 0.1;
        @. vFuncs.Γ = (1 - damping) * vFuncs.Γ + damping * vFuncsNew.Γ

        # (3) 정규화, 수렴도 계산
        vFuncs.Γ ./= sum(vFuncs.Γ)
        vFuncsNew.Γ ./= sum(vFuncsNew.Γ)   # 추가 normalize

        diffs = abs.(vFuncsNew.Γ .- vFuncs.Γ)
        maxdiffs = maximum(diffs)

        if mod(iter, 100) == 0
        println("iter=$iter, sumΓ=", sum(vFuncs.Γ), " maxdiffs=", maxdiffs)
            println("stationary_distribution!: iter=$iter, maxdiffs=$(round(maxdiffs, sigdigits=6))")
        end

        iter += 1
    end

    println("✅ stationary_distribution! 종료: iter=$(iter-1), maxdiffs=$(round(maxdiffs, sigdigits=6))")
    

    return nothing
end
