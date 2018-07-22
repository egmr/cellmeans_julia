#!/bin/julia
# Type for containing marker and genotype data
type Markers
    markerInfo::Array{String,2} # m x 5 columns
    IID::Array{String,1} # n x 1
    genotypes::Array{Float64,2} # n x m
end

# Type for containing phenotype data
type Phenotypes
    IID::Array{String,1}
    trait::Array{Float64,2} # n x k
end


# Type for output summary statistics
type Summary
    markerInfo::Array{String,2} # 5 columns
    f::Array{Float64,1}
    stats::Array{Float64,2} # 8 columns
end

# Type for Alpha data
type Alphas
    markerInfo::Array{String,2} # m x 5
    values::Array{Float64,2} # k x m
end

## Perhaps MarkerData should be a subtype
