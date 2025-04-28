
pwd()
using JuMP
using Ipopt
using Interpolations
using ForwardDiff

bgrid_start = 0.0
bgrid_stop = 1.0
bgrid_npts=5
typeof(range(bgrid_start,stop=bgrid_stop,length=bgrid_npts))
AA = [1, 1.2, 30]
A = (1, 1.2, 20)
typeof(A)
AA[1] = 1.2
AA



## example for nonlinear programming with JuMP and Ipopt
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
    return x,y
end

returning_optimizer(1, 2)
returning_optimizer(1, 2)[1] + returning_optimizer(1, 2)[2] 
returning_optimizer(20, 300)
returning_optimizer(0.5, 0.2)

## example for interpolation
x = 0:0.5:3
x[1]
typeof(x), size(x)
y = @. sin(x)

itp = interpolate( (x,), y, Gridded(Linear()))
x_star = 0.73
itp(0.73)

## example for interpolating a 3-D function

xx = 0:0.5:2
size(xx)
yy = 0:0.25:1
size(yy)
zz = 0:π/6:π/2
size(zz)

f(x,y,z) = sin(x)*cos(y) + z^2
fvals = [sin(xi)*cos(yi) + zi^2 for zi in zz, yi in yy, xi in xx]

itp_3d = interpolate((zz,yy,xx), fvals, Gridded(Linear()))
itp_3d(1.5, 0.18, 0.73) ## the evaluating point should be in the grid of the interpolation
itp_3d.(0:π/4:π/2, 0:0.5:1, 0:1:2)
xn = 0:1:2
yn = 0:0.5:1
zn = 0:π/4:π/2
F = [itp_3d(z, y, x) for z in zn, y in yn, x in xn]

## multidemensional uniformly spaced grid BSpline interpolation
A_x1 = 1:.1:10
A_x2 = 1:.5:20
f(x1, x2) = log(x1+x2)
A = [f(x1,x2) for x1 in A_x1, x2 in A_x2]
itp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))
sitp = Interpolations.scale(itp, A_x1, A_x2)
sitp(5., 10.) # exactly log(5 + 10)
log(5+10)
sitp(5.6, 7.1) # approximately log(5.6 + 7.1)
log(5.6 + 7.1) # exactly log(5.6 + 7.1)

eps()


## example of optimizing interpolated function
# given xx, yy, zz and fvals 
itp_3d_Bspline = interpolate((zz,yy,xx), fvals, Gridded(Linear()))
f_interp(z, y, x) = itp_3d_Bspline(z, y, x) # a function that evalutes at u using interpolation
f_interp(z::T, y::T, x::T) where {T<:Real} = itp_3d_Bspline(z, y, x)


function optimize_example()
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "tol", 1e-8)
    register(model, :foo, 3, f_interp; autodiff = true)
    (zlb, zub) = extrema(zz)
    (ylb, yub) = extrema(yy)
    (xlb, xub) = extrema(xx)
    @variable(model, zlb <= z <= zub, start = 0.5)
    @variable(model, ylb <= y <= yub, start = 0.5)
    @variable(model, xlb <= x <= xub, start = 0.5)
    @NLobjective(model, Min, foo(z, y, x))
    optimize!(model)
    return termination_status(model), objective_value(model), value.(z, y, x)
end

optimize_example()

extrema(zz), extrema(yy), extrema(xx)
##########################################################################################################

function minmax_broadcast_example()
    a = [1,2,3]
    b = [10,20]
    c = [1,2,3,4,5]
    d = [1,2,3]
    e = [1,2,3,4]

    f(x,y,z,w,u) = max(min(x/y, y/z), (x+y)/z -w + u )
    fvals = [f.(ia, ib, ic, id, ie) for ia in a, ib in b, ic in c, id in d, ie in e]
    return fvals
end

minmax_broadcast_example()[1,1,1,1, :]

x = ones(2,3,4)
x .= x * 10
x

f_example(x,y,z) = x + y + z
f_example.(1,1,[1,2,3])


v = [10, 20, 30]                 # Vector
for i in eachindex(v)
    @show i, v[i]                # i is 1,2,3            (Int)
end

a = 10:20
for (iλ, λprime) in pairs(a)
    @show iλ, λprime             # iλ is 1,2,3...10      (Int)
end
x =  1 < 2 ? "Fail" : "yes" 

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