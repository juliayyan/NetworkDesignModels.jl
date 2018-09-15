module ColumnGeneration2

    using TransitNetworks, NetworkDesignModels, Gurobi, JuMP
    using DataFrames, JLD, JLD2
    using Base.Test

    # read network
    pathfiles = [
        "data/tidy/toy/0-stationpaths/stationpath$(k).csv" for k in 1:8
    ]
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
        budget = sum(NetworkDesignModels.linecost(np,line,gridtype)
                     for line in np.lines)
        rmp = NetworkDesignModels.MasterProblem(
            np, 
            linelist = [[1,6]],
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
                # solve master
                NetworkDesignModels.optimize(rmp, budget)
                push!(primal_objs, getobjectivevalue(rmp.model))
                
                sp = NetworkDesignModels.SubProblem(rmp,
                    direction = direction, 
                    delta = 0.8);

                # generate new column
                p = getdual(rmp.choseline);
                q = getdual(rmp.bcon);
                s = getdual(rmp.choseub)
                path = NetworkDesignModels.generatecolumn(sp, p, q, 
                    coeffs = NetworkDesignModels.spcoeffs(rmp,sp))
                if length(path) == 0
                    break
                end

                push!(dual_objs, getobjectivevalue(sp.model))     
                NetworkDesignModels.addcolumn!(rmp, path)
            end
        end
        @test all(rmp.linelist[2] .== [3,1,2,5,13,11,9,10,12])
        @test isapprox(dual_objs[1], 7400.0)
        @test isapprox(dual_objs[2], 7400.0)
        @test isapprox(dual_objs[3], 5400.0)
        @test isapprox(dual_objs[4], 2600.0)
        for i in 1:3
            @test isapprox(primal_objs[i+1], primal_objs[i] + dual_objs[i])
        end
        @test isapprox(primal_objs[end], sum(np.odmatrix))

    end 

end
