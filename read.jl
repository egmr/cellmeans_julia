#!/bin/julia
# read.jl
## For functions used to read in datasets
# Near double up of some functions
using CSV, DataFrames, GeneticVariation
# for gz files
## using CodecZlib
## reader = VCF.Reader(GzipDecompressionStream(open("example.vcf.gz"))
## Genotypes (vcf), Phenotypes, and Alphas
#
# Reads using CSV.read
function read_genotypes(file, centered)
    # Fails if first line is header :O
    header = find_header(file)
    fstream = open(file)
    for i in 1:(header-1)
        readline(fstream)
    end
    ncol = length(split(readline(fstream), "\t"))
    close(fstream)
    filetypes = Array{DataType}(ncol)
    fill!(filetypes, String)
    vcf = CSV.read(file; delim = "\t", datarow =header+1, header=header, types= filetypes)
    ids = convert(Array{String,1}, names(vcf)[10:end])
    marker_info = vcf[:,1:5]
    genos = calc_genotypes(vcf, length(ids), size(vcf,1))

    if centered
        geno_means = colMeans(genos)
        genox = genos
        genos = center(genox, geno_means)
    end
    return Markers(marker_info, ids, genos)
end
#
function read_alphas(afile, info)
    alphas = Array{Float64,2}(CSV.read(afile, delim = "\t", header = true))
    return Alphas(info, alphas)
end
#
function read_phenotypes(pfile, y_types)
    y = CSV.read(pfile, delim = " ", header = true, types = y_types)
    return Phenotypes(y[:,2], y[:,3:3])
end
# reads in residual variance values. returns a transposed k x 1 [ => 1 x k]
function read_residuals(rvfile)
    rv = CSV.read(rvfile, delim = "\t", header = true, types = [Float64])
    return (Array{Float64, 2}(rv))'
end
#
#
#
## Read VCF using GeneticVariation.jl
# Looks to be faster than above
function read_genotypes_vcfreader(gfile, centered)
    line_count = countlines(gfile)
    header_count = find_vcfheader(gfile)
    m = line_count - header_count # number of markers in file
    #
    reader = VCF.Reader(open(gfile, "r"))
    ids = header(reader).sampleID
    n = size(ids, 1) # number of individuals
    genotypes = zeros(Int64, n, m)
    markerInfo = Array{String, 2}(m,5)
    #
    i = 1
    for record in reader
        genotypes[:,i] = calc_geno(n, VCF.genotype(record, 1:n, "GT"))
        markerInfo[i,:] = get_markerInfo(record)
        i += 1
        if i >= m
            break
        end
    end
    close(reader)
    if centered
        geno_means = colMeans(genotypes)
        genotypes = center(genotypes, geno_means)
    end
    println("END")
    return Markers(markerInfo, ids, genotypes)
end

function find_vcfheader(filename)
    header = 1
    file = open(filename, "r")
    while( readline(file)[1:6] != "#CHROM")
        header += 1
    end
    return header
end
# Converts 0/0,0/1,1/1 genotypes to 0,1,2 allele counts
function calc_geno(n, genotypes)
    record_genos = zeros(Int64, n)
    for i in 1:n
        record_genos[i] = (genotypes[i][1] == '1') + (genotypes[i][3] == '1')
    end
    return record_genos
end

function get_markerInfo(record)
    c = VCF.chrom(record)
    p = string(VCF.pos(record))
    i = VCF.id(record)[1]
    r = VCF.ref(record)
    a = "."
    try
        a = VCF.alt(record)[1]  # Where alternate genotype is '.', returns error
    catch
        a = "."
    end
    #println([c,p,i,r,a])
    return [c,p,i,r,a]
end
