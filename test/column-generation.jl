module ColumnGeneration

    using TransitNetworks, NetworkDesignModels, Gurobi, JuMP
    using DataFrames, JLD, JLD2
    using Base.Test

    # read network
    np = load("data/processed/networks/9b-transitnetwork.jld2", "keynetwork")
    stns = np.stnnames
    nstns = length(stns)
    stnindex = Dict(zip(stns, 1:nstns));
    
    # read o-d demand
    demand = DataFrames.readtable(
        "data/processed/odmatrices/odmatrices-output-bus-linelevel.csv"
    );
    demand[:origin]      = string.(demand[:origin])
    demand[:destination] = string.(demand[:destination])

    # read stopclusters
    stopclusters = DataFrames.readtable("data/tidy/1-stopclusters.csv");
    stopclusters = Dict(zip(stopclusters[:stop_id], stopclusters[:cluster_stop]))
    for u in setdiff(unique(vcat(demand[:origin], demand[:destination])), collect(keys(stopclusters)))
        stopclusters[u] = u
    end 
            
    @testset "Loading Network" begin
        @test nstns == 410
        @test length(np.lines) == 46
        @test sum(demand[:demand]) == 364189
        @test nrow(demand) == 116798
    end 

    # cluster stops of odmatrix
    demand[:origin]      = [stopclusters[u] for u in demand[:origin]]
    demand[:destination] = [stopclusters[u] for u in demand[:destination]]
    demand = by(demand, [:origin, :destination, :hour], d -> DataFrame(demand = sum(d[:demand])))
    demand = demand[demand[:origin] .!= demand[:destination],:]
    @testset "Clustering Stations" begin
        @test sum(demand[:demand]) == 364124
        @test nrow(demand) == 101758
    end 

    # create odmatrix
    demand = demand[findin(demand[:origin], stns),:];
    demand = demand[findin(demand[:destination], stns),:];
    odmatrix = zeros(nstns, nstns);
    TransitNetworks.aggregateODdemand!(demand, stnindex, 0:23, odmatrix);
    odmatrix[odmatrix .< 10] = 0;
    np.odmatrix = odmatrix;
    @testset "Create O-D Matrix" begin
        @test sum(np.odmatrix) == 124285
    end 
    
    @testset "Generating Columns" begin
        rmp = NetworkDesignModels.MasterProblem(
            np, 
            linelist = copy(np.lines), 
            initialbudget = 100) # corresponds to about 1/3 of the key budget (296)
        solve(rmp.model)
        @test JuMP.getobjectivevalue(rmp.model) == 93587.42418280317

        sp = NetworkDesignModels.SubProblem(np);
        p = getdual(rmp.choseline);
        q = getdual(rmp.bcon);
        path = NetworkDesignModels.generatecolumn(sp, p, q)
        @test path == [8, 27, 53, 175, 177, 179, 182, 255]
        @test JuMP.getobjectivevalue(sp.model) == 294.650731361031

        lines = copy(np.lines)
        push!(lines, path)
        push!(lines, reverse(path))
        rmp = NetworkDesignModels.MasterProblem(np, 
            linelist = lines, 
            initialbudget = 100)
        solve(rmp.model)
        @test JuMP.getobjectivevalue(rmp.model) == 93882.0749141642    
    end 

end
