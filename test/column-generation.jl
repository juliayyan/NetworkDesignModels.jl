module ColumnGeneration

    using TransitNetworks, NetworkDesignModels, Gurobi, JuMP
    using CSV, FileIO
    using Test

    # read network
    np = load("data/processed/networks/9b-transitnetwork.jld2", "keynetwork")
    stns = np.stnnames
    nstns = length(stns)
    stnindex = Dict(zip(stns, 1:nstns));
    
    # read o-d demand
    demand = CSV.read(
        "data/processed/odmatrices/odmatrices-output-bus-linelevel.csv"
    );
    demand[:origin]      = string.(demand[:origin])
    demand[:destination] = string.(demand[:destination])

    # read stopclusters
    stopclusters = CSV.read("data/tidy/1-stopclusters.csv");
    stopclusters = Dict(zip(
        string.(stopclusters[:stop_id]), string.(stopclusters[:cluster_stop])
    ))
    for u in setdiff(unique(vcat(demand[:origin], demand[:destination])),
                     collect(keys(stopclusters)))
        stopclusters[u] = u
    end 

    @testset "Loading Network" begin
        @test nstns == 410
        @test length(np.lines) == 46
        @test sum(demand[:demand]) == 364189
    end 

    # cluster stops of odmatrix
    demand[:origin]      = [stopclusters[u] for u in demand[:origin]]
    demand[:destination] = [stopclusters[u] for u in demand[:destination]]
    demand = by(demand, [:origin, :destination, :hour],
                d -> DataFrame(demand = sum(d[:demand])))
    demand = demand[demand[:origin] .!= demand[:destination],:]
    @testset "Clustering Stations" begin
        @test sum(demand[:demand]) == 364124
    end 

    # create odmatrix
    demand = demand[findall((in)(stns), demand[:origin]),:];
    demand = demand[findall((in)(stns), demand[:destination]),:];
    odmatrix = zeros(nstns, nstns);
    TransitNetworks.aggregateODdemand!(demand, stnindex, 0:23, odmatrix);
    odmatrix[odmatrix .< 10] = 0;
    np.odmatrix = odmatrix;
    @testset "Create O-D Matrix" begin
        @test sum(np.odmatrix) == 124285
    end 
    
    @testset "Generating Columns" begin
        budget = 100.0 # corresponds to about 1/3 of the key budget (296)
        rmp = NetworkDesignModels.MasterProblem(np)
        NetworkDesignModels.optimize(rmp, budget)
        @test length(rmp.linelist) == 23
        @test isapprox(JuMP.getobjectivevalue(rmp.model),93587.42418280317)

        sp = NetworkDesignModels.SubProblem(rmp);
        p = getdual(rmp.choseline);
        q = getdual(rmp.bcon);
        s = getdual(rmp.choseub);
        path = NetworkDesignModels.generatecolumn(sp, p, q)
        @test path == [227, 166, 264, 342, 98, 99, 22, 
                       36, 19, 31, 24, 16, 253, 201, 159, 
                       349, 141, 385, 25, 46, 364, 317, 
                       139, 359, 363, 371, 160, 377, 305, 331]
        @test isapprox(JuMP.getobjectivevalue(sp.model), 2316.560035938597)
        @test sum(getvalue(sp.edg)) == length(path)-1

        NetworkDesignModels.addcolumn!(rmp, path)
        NetworkDesignModels.optimize(rmp, budget)
        @test isapprox(JuMP.getobjectivevalue(rmp.model), 95674.68627657647)
        @test length(rmp.linelist) == 24
    end 

    @testset "Two-Leg Commute Model" begin
        rmp = NetworkDesignModels.MasterProblem(np, 
            linelist = [np.lines[1]], 
            nlegs = 1);
        NetworkDesignModels.optimize(rmp, 5.0)
        @test JuMP.getobjectivevalue(rmp.model) == 10764.0

        rmp = NetworkDesignModels.MasterProblem(np, 
            linelist = [np.lines[1][1:10],
                        np.lines[1][10:end]], 
            nlegs = 2);
        NetworkDesignModels.optimize(rmp, 5.0)
        @test JuMP.getobjectivevalue(rmp.model) == 10764.0
    end

end
