pwd()
readdir()

include("C:\\Users\\master\\Desktop\\(2025) banking regulation\\BankResolution\\main.jl")
include("C:\\Users\\master\\Desktop\\(2025) banking regulation\\BankResolution\\parameters.jl")

params = Initiate_Params(qd,β,Rf,wr,α,ρ,g,ξ,cF,dBar,σ,τC,z,α1,α2,α3,δL,δM,δH,cM,cO,cL,H,Γ,λL,λM,λH,γ,ϕ,n_start,n_npts,n_stop,l_start,l_npts,l_stop,s_start,s_npts,s_stop,b_start,b_npts,b_stop)

Rl = 1.03;
regime = false
VFI(params, Rl, regime, 10, 0.001)



params.nGrid
params.lGrid
Iterobj_is = Initiate_MatrixIterObj_i(params); 
Iterobj_is[1,1]

Threads.@threads for (iN, n) in pairs(params.nGrid) # for each n 
    println(iN), println(n)
end