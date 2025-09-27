pwd()

include("parameters.jl");
include("main_v2_250922.jl");

Random.seed!(5877);
params = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,α1,α2,α3,δL,δM,δH,cM,cO,cL,ϵ,H,F,M,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)

startt = 1.0+eps();
endd = 2.0-eps();
# @time sol_false = solve_model(params, false, startt, endd);

function boundary_check()

    solve_model_given_r_single = Rl -> solve_model_given_r(Rl; Params = params, Regime = false)
    
    @show a = 1 + eps()
    @show b = 2 - eps()
    A = solve_model_given_r_single(a)
    B = solve_model_given_r_single(b)
    x = A * B
    return (A, B, x)
end

@show boundary_check()

# l_min, l_max = first(params.lGrid), 322.0 # (last(params.lGrid)+first(params.lGrid))/2, last(params.lGrid)
# s_min, s_max = first(params.sGrid), last(params.sGrid)
# b_min, b_max = first(params.bGrid), last(params.bGrid)
# @show l_cap = min.(l_max, (params.β - params.g) * params.deltaGrid) 
# @show b_cap = min.(b_max, (1 - params.α * params.wr) * l_cap[1] .+ s_max .- (1 - params.g) * params.deltaGrid[1]) # b_start = (b_max+b_min)/2
# @show b_cap = min.(b_max, (1 - params.α * params.wr) * l_cap[2] .+ s_max .- (1 - params.g) * params.deltaGrid[2]) # b_start = (b_max+b_min)/2
# @show b_cap = min.(b_max, (1 - params.α * params.wr) * l_cap[3] .+ s_max .- (1 - params.g) * params.deltaGrid[3]) # b_start = (b_max+b_min)/2



# @show l0 = (3 .* l_min .+ l_cap)/4
# @show s0 = (3 * s_min + s_max)/4
# @show b0 = (3 * b_min + b_max)/4

# @show (1-params.α*params.wr)*l0 .+ s0 .- b0 .>= (1-params.g)*params.deltaGrid
# see if end points satisfy root finding condition 
function experiment()
    regime = false 

    function solve_bank_problem(params::Params{T,S},vFuncs::VFuncs{T,S},iDelta::S,iLambda::S,iN::S,G::Array{T,3}; warm_start::Union{Nothing,NTuple{3,Float64}}=nothing) where {T<:Real, S<:Integer}

            lRange = range(first(params.lGrid), last(params.lGrid), length(params.lGrid))
            sRange = range(first(params.sGrid), last(params.sGrid), length(params.sGrid))
            bRange = range(first(params.bGrid), last(params.bGrid), length(params.bGrid))
    
            Gmin = minimum(G)
            Gmax = maximum(G)
            Gnorm = (G .- Gmin) ./ max(Gmax - Gmin, eps(T)) # scaling 
            G_itp_rev = interpolate(Gnorm, BSpline(Cubic(Line(OnGrid()))))
            G_itp_rev = Interpolations.scale(G_itp_rev, lRange, sRange, bRange)
            G_itp_ext_rev = extrapolate(G_itp_rev, Line())
            G_interp_rev(l, s, b) = G_itp_ext_rev(l, s, b)
    
            model = Model(Ipopt.Optimizer)
            set_optimizer_attribute(model, "tol", 1e-6)
            set_optimizer_attribute(model, "acceptable_tol", 1e-3)
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
            set_optimizer_attribute(model, "max_wall_time", 30.0) # 시간 방어막: iteration_limit 대신 acceptable로 떨어지게 유도
            set_optimizer_attribute(model, "print_options_documentation", "yes")

            JuMP.register(model, :G_interp_rev, 3, G_interp_rev; autodiff = true)
            l_min, l_max = first(params.lGrid), last(params.lGrid) # (last(params.lGrid)+first(params.lGrid))/2, last(params.lGrid)
            s_min, s_max = first(params.sGrid), last(params.sGrid)
            b_min, b_max = first(params.bGrid), last(params.bGrid)
            l_cap = min(l_max, (params.β - params.g) * params.deltaGrid[iDelta]) # l_start = (l_max+l_min)/2

              # start: warm-start 있으면 사용, 없으면 중앙값
            @show l0 = warm_start === nothing ? (l_min + l_cap)/2 : warm_start[1]
            @show l0
            @show s0 = warm_start === nothing ? (s_min + s_max)/2 : warm_start[2]
            @show b0 = warm_start === nothing ? (b_min + b_max)/2 : warm_start[3]

            # println("c1 residual = ", (1-params.α*params.wr)*l_start + s_start - b_start - (1-params.g)*params.deltaGrid[iDelta])
            # println("c2 residual = ", l_start - (params.β - params.g)*params.deltaGrid[iDelta])
            @variable(model, l_min <= l <= l_cap, start = clamp(l0, l_min, l_cap))
            @variable(model, s_min <= s <= s_max, start = clamp(s0, s_min, s_max))
            @variable(model, b_min <= b <= b_max, start = clamp(b0, b_min, b_max))
            # println("G_interp at start = ", G_interp_rev(l_start, s_start, b_start))
            # grad = ForwardDiff.gradient(u -> G_interp_rev(u[1], u[2], u[3]), [l_start, s_start, b_start])
            # println("gradient at start = ", grad)
    
            @NLobjective(model, Max, G_interp_rev(l, s, b))
            @constraint(model, (1-params.α*params.wr)*l + s - b >= (1-params.g)*params.deltaGrid[iDelta]) # constratint 1
           # @constraint(model, l <= (params.β - params.g)*params.deltaGrid[iDelta]) # constratint 2

           # turn on logging
            set_silent(model, false)
            # set_optimizer_attribute(model, "print_level", 12)

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
            else
                error("Optimization did not converge: status=$stat, primal=$(primal_status(model))")
            end
    end
    
    function boundary_check()

        solve_model_given_r_single = Rl -> solve_model_given_r(Rl; Params = params, Regime = regime)
    
        @show a = 1 + eps()
        @show b = 2 - eps()
        A = solve_model_given_r_single(a)
        B = solve_model_given_r_single(b)
        x = A * B
        return (A, B, x)
    end

    return boundary_check()
 end