#!/bin/julia
# read_assist.jl

## For functions which assist the functions in read.jl
## For example: sorting by IIDs, finding headers, calculating genotypes
function sort_cellmeans(cellmeans, y_iid)
    new_order = get_new_order(y_iid, cellmeans.IID)
    cellmeans.genotypes = cellmeans.genotypes[new_order,:]
    cellmeans.IID = cellmeans.IID[new_order]
    return cellmeans
end
#
function sort_phenotypes(phenotypes, m_iid)
    new_order = get_new_order(m_iid, phenotypes.IID)
    phenotypes.trait = phenotypes.trait[new_order,:]
    phenotypes.IID = phenotypes.IID[new_order]
    return phenotypes
end
# find which line the header of the VCF is # REQUIRES EDIT
function find_header(file)
    # edit if header is line 1!!!
    header = 1
    x = CSV.Source(file)
    while(CSV.readline(x)[1:6] != "#CHROM")
        header += 1
    end
    return header + 1
end
# calculate genotypes from
function calc_genotypes(vcf, nIds, nMarkers)
    genotypes = zeros(Int64, nIds, nMarkers)
    for i in 1:nMarkers
        for j in 10:(9+nIds)
            genotypes[(j-9),i] = (vcf[i,j][1] == '1') + (vcf[i,j][3] == '1')
        end
    end
    return genotypes
end
#
function get_new_order(ordered, x)
    ordered_n = size(ordered,1)
    x_n = size(x,1)
    new_order = zeros(Int64,ordered_n)
    for i in 1:x_n
        for j in 1:ordered_n
            if x[i] == ordered[j]
                new_order[j] = i
            end
        end
    end
    return new_order
end

# calculate genotype means for matrix
function colMeans(x)
    colMeans = mean(x,1)
end
# center genotype matrix
function center(x, colMeans)
    return (x .- colMeans)
end
