# use Krylov algorithm to get both left and right eigenvectors of non-hermitian H, with MPO representation

include("../../src/spin_hamiltonian.jl")
include("../../src/NHKPM.jl")
include("../../src/KrylovSchur.jl")

using LinearAlgebra

using PyCall
plt = pyimport("matplotlib.pyplot")

N=8
J=[1,0,0,1,0,0]
B=[-1*im,0,0,1*im,0,0]

#mode="MPO"
mode="sparse_tensor"

H=get_generalized_SSH_spin_chain(J,B,N=N,mode=mode)
H2=get_H_dag_MPO(H)

sites=firstsiteinds(H)  
psi0 = randomMPS(sites,1)

energy,psi,err=get_eigs_krylovMPO(H,psi0,which="SR",k=1,m=5,mode=mode,tol=1e-4,maxerror=1e-4,max_iter=10,maxdim=maxdim,factor=100)
H2=get_H_dag_MPO(H)
energy2,psi2,err2=get_eigs_krylovMPO(H2,psi0,which="SR",k=1,m=5,mode=mode,tol=1e-4,maxerror=1e-4,max_iter=10,maxdim=maxdim,factor=100)

println(energy)
println(energy2)
println(dot(psi,psi2))  # in general this is not 1 since psi and psi2 are nomalized by themselves


