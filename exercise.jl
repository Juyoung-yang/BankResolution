
pwd()
using JuMP
using Ipopt

## example 
# 1. create a JuMP model that uses Ipopt as the solver
model = Model(Ipopt.Optimizer) 
set_optimizer_attribute(model, "tol", 1e-8)

# 2. Decision variables 
@variable(model, x >= 0)
@variable(model, y >= 0)

# 3. nonlinear objective function
@NLobjective(model, Min, (1-x)^2 + (y-2)^2)

# 4. Nonlinear constraints
@NLconstraint(model, x*y >= 1)
@NLconstraint(model, x^2 + y^2 <= 5) 

# 5. solve
print(model)
optimize!(model)

# 6. extract results
is_solved_and_feasible(model)
set_optimizer(model, Ipopt.Optimizer)

println(termination_status(model))
println("objective_value: ", objective_value(model))
println("x: ", value(x))
println("y: ", value(y))

function returning_optimizer(a, b)
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "tol", 1e-8)
    @variable(model, x >= 0)
    @variable(model, y >= 0)
    f(x,y) = (a-x)^2 + (y-b)^2
    @NLobjective(model, Min, f(x,y))
    @NLconstraint(model, x*y >= 1)
    # @NLconstraint(model, x^2 + y^2 <= 5) 
    optimize!(model)
    x = value(x)
    y = value(y)
    return (x,y) 
end

returning_optimizer(1, 2)
returning_optimizer(20, 300)
returning_optimizer(0.5, 0.2)



##########################################################################################################
using Optim
f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2

x0 = [0.0, 0.0]
optimize(f, x0)


function example_rosenbrock()
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, x)
    @variable(model, y)
    @objective(model, Min, (1 - x)^2 + 100 * (y - x^2)^2)
    @show optimize!(model)
    assert_is_solved_and_feasible(model)
    Test.@test objective_value(model) ≈ 0.0 atol = 1e-10
    Test.@test value(x) ≈ 1.0
    Test.@test value(y) ≈ 1.0
    return
end

@show example_rosenbrock()

function example_clnlbeam()
    N = 1000
    h = 1 / N
    alpha = 350
    model = Model(Ipopt.Optimizer)
    @variables(model, begin
        -1 <= t[1:(N+1)] <= 1
        -0.05 <= x[1:(N+1)] <= 0.05
        u[1:(N+1)]
    end)
    @objective(
        model,
        Min,
        sum(
            0.5 * h * (u[i+1]^2 + u[i]^2) +
            0.5 * alpha * h * (cos(t[i+1]) + cos(t[i])) for i in 1:N
        ),
    )
    @constraint(
        model,
        [i = 1:N],
        x[i+1] - x[i] - 0.5 * h * (sin(t[i+1]) + sin(t[i])) == 0,
    )
    @constraint(
        model,
        [i = 1:N],
        t[i+1] - t[i] - 0.5 * h * u[i+1] - 0.5 * h * u[i] == 0,
    )
    optimize!(model)
    println("""
    termination_status = $(termination_status(model))
    primal_status      = $(primal_status(model))
    objective_value    = $(objective_value(model))
    """)
    assert_is_solved_and_feasible(model)
    return
end

example_clnlbeam()