# NetworkDesignModels

Code supplement for the paper "Data-driven transit network design at scale" by Dimitris Bertsimas, Yee Sian Ng, and Julia Yan.

## Introduction

Given a set of stations with known demand between pairs of stations (u, v), this code is intended to design a transit network to service the demand subject to a budget constraint.

## Installation

This code was written in the Julia programming language.  You can install a recent version of Julia from the [downloads page](http://julialang.org/downloads/). The most recent version of Julia at the time this code was last tested was Julia 1.1.1.

Several packages must be installed in Julia before the code can be run. These packages include JuMP, Gurobi, JLD2, DataFrames, CSV, ProgressMeter, and others listed in `src/NetworkDesignModels.jl`.

These packages can be added with the following commands:

```pkg> add PACKAGE_NAME```

This code does not use the most recent version of JuMP and requires JuMP v0.18.6 to run.  For JuMP, the installation command is:

```pkg> add JuMP@0.18.6```

## Scripts

Each script in the ```scripts/``` directory is labeled with the section number of the paper that it corresponds to.  To run a script, navigate to the directory where this package is installed, and run the following:

```julia> include("scripts/SCRIPT_NAME.jl")```

Data for the scripts are in the ```scripts/data``` directory, saved in the JLD2 format.

The general format for solving a problem is as follows:

Load a transit network data structure from the data called `np`:

```
JLD2.@load "data/transit-network-boston.jld2" np 
```

Build and solve the master problem:

```
rmp = NDM.MasterProblem(np); NDM.optimize(rmp, budget);
```

Build and solve the subproblem:

```
sp = NDM.SubProblem(rmp); path = NDM.generatecolumn(sp, rmp);
```

and repeat.  Each function comes with various parameters that are explained in the comments.
