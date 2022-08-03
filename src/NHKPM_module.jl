__precompile__()

module NHKPM

using ITensors
using LinearAlgebra
using DelimitedFiles
using Statistics
using Dates

using Random
using SparseArrays
using Arpack
using BenchmarkTools

using PyCall
const plt=PyNULL()
const colormap=PyNULL()
function __init__()
    copy!(plt,pyimport_conda("matplotlib.pyplot","matplotlib"))
    copy!(colormap,pyimport_conda("matplotlib.colors","matplotlib"))
end

include("NHKPM.jl")
include("KrylovSchur.jl")
include("NH_hamiltonian.jl")
include("spin_hamiltonian.jl")

end
