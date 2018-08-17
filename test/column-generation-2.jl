module ColumnGeneration2

    using TransitNetworks, NetworkDesignModels, Gurobi, JuMP
    using DataFrames, JLD, JLD2
    using Base.Test

    # read network
    pathfiles = ["data/tidy/toy/0-stationpaths/stationpath$(k).csv" for k in 1:8]
    segfiles = ["data/tidy/toy/0-stationpaths/segment$(k).csv" for k in 1:8]
    stns, paths = TransitNetworks.readstationpaths(pathfiles)
    segs, seglist = TransitNetworks.readsegmentpaths(segfiles)
    latlngs = readtable("data/tidy/toy/1-stops.csv")
    nstns = length(stns)
    demand = 100*(ones(nstns,nstns) - eye(nstns))

    # construct network
    np = TransitNetworks.TransitNetworkProblem(paths, segs, latlngs, demand);
    gridtype = :latlong

    @testset "Loading Network" begin
        @test np.nstations == 16
        @test length(np.lines) == 8
        @test sum(np.odmatrix) == 24000
    end 
rmp = nothing
primal_objs = nothing
dual_objs = nothing
    @testset "Generating Columns" begin
        budget = sum(NetworkDesignModels.linecost(np,line,gridtype) for line in np.lines)
        rmp = NetworkDesignModels.MasterProblem(
            np, 
            linelist = [[1,6]],
            initialbudget = budget,
            nlegs = 2,
            transferparam = -0.25,
            gridtype = gridtype
        );
        NetworkDesignModels.optimize(rmp, budget)
        @test getobjectivevalue(rmp.model) == 200

        primal_objs = Float64[]
        dual_objs = Float64[]
        for direction in [[1.0,0.0],[1.0,1.0],[0.0,1.0],[-1.0,1.0]]
            for i in 1:10
                sp = NetworkDesignModels.SubProblem(rmp,
                    direction = direction, 
                    delta = 0.8);

                # solve master
                NetworkDesignModels.optimize(rmp, budget)
                push!(primal_objs, getobjectivevalue(rmp.model))
                
                # generate new column
                p = getdual(rmp.choseline);
                q = getdual(rmp.bcon);
                s = getdual(rmp.choseub)
                path = NetworkDesignModels.generatecolumn(sp, p, q, s)
                if length(path) == 0
                    break
                end
                
                push!(dual_objs, getobjectivevalue(sp.model))     
                NetworkDesignModels.addcolumn!(rmp, path)
            end
        end
        primal_vals = [200.0,7600.0,14800.0,18723.813395198402,20689.8510808549,
                       21590.84841034586,22566.451740632667,22750.64905711167,
                       23061.637598663154,23231.961790427915,23231.961790427915,
                       23466.951717820542,23560.065140604653,23560.065140604653,
                       23603.07854663059,23628.652205095023,23630.450957003137,
                       23630.450957003137,23633.210666631596,23633.210666631596,
                       23908.28961173809,23938.58392475759,23974.69252174918,
                       23974.69252174918,24000.0,24000.0]
        dual_vals = [7300.0, 5400.0, 5400.0, 1340.82, 1317.99,
                     1524.59, 827.837, 567.095, 429.145, 383.678, 
                     518.713, 347.696, 143.485, 119.274, 96.5043,
                     73.6098, 61.2911, 25.9272, 544.472, 429.654, 
                     238.853, 226.258, 303.004]
        @test isapprox(primal_vals, primal_objs, atol = 0.01)                     
        @test isapprox(dual_vals, dual_objs, atol = 0.01)

    end 

end
