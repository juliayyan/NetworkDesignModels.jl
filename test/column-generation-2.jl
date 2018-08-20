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
        primal_vals = [200.0, 7600.0, 15200.0, 20317.3, 20606.9, 
                       21142.0, 22445.1, 22445.1, 22706.7, 23239.6, 
                       23563.6, 23902.2, 23902.2, 23902.2, 23902.2, 
                       23902.2, 23902.2, 23909.6, 24000.0, 24000.0]
        dual_vals = [7300.0, 5600.0, 4900.0, 3118.26, 1136.21, 
                     1076.77, 573.208, 552.95, 494.655, 427.498, 
                     155.546, 119.092, 26.3887, 131.487, 32.6297, 
                     150.499, 214.456]
        @test isapprox(primal_vals, primal_objs, atol = 0.1)                     
        @test isapprox(dual_vals, dual_objs, atol = 0.1)

    end 

end
