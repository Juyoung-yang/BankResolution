


function graph2()
    theme(:bright) # theme(:ggplot2)
    plot(params.bGrid[7:15]./(params.lGrid[10]+params.sGrid[10]), eq.vFuncs.qBond[1,19,7:15,2,1], label ="low loan-to-safe asset", legend =:bottomleft, legendfontsize = 12, guidefontsize = 12, lw = 3, xlabel = "Bank debt/ total asset", ylabel = "Debt price")
    plot!(params.bGrid[7:15]./(params.lGrid[15]+params.sGrid[5]), eq.vFuncs.qBond[19,1,7:15,2,1], label ="high loan-to-safe asset", lw = 3)
    plot!(params.bGrid[7:15]./(params.lGrid[10]+params.sGrid[10]), eq_false.vFuncs.qBond[1,1,1,1,1] .*ones(size(eq.vFuncs.qBond[10,10,7:15,2,3])), label ="benchmark", ls =:dash, lw = 3)
    savefig("plot_highres2.pdf")
end

function graph1()
    theme(:bright) # theme(:ggplot2)
    plot(params.bGrid[7:15]./(params.lGrid[10]+params.sGrid[10]),eq.vFuncs.qBond[10,10,7:15,2,1], label ="low λ", lw = 3)
    plot!(params.bGrid[7:15]./(params.lGrid[10]+params.sGrid[10]),eq.vFuncs.qBond[10,10,7:15,2,3], label ="high λ", lw = 3)
    plot!(params.bGrid[7:15]./(params.lGrid[10]+params.sGrid[10]), eq_false.vFuncs.qBond[1,1,1,1,1] .*ones(size(eq.vFuncs.qBond[10,10,7:15,2,3])), label ="benchmark", ls =:dash, lw = 3)
    savefig("plot_highres.pdf")
end
