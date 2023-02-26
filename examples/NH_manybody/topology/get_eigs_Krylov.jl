include("../../../src/spin_hamiltonian.jl")
include("../../../src/NHKPM.jl")
include("../../../src/KrylovSchur.jl")

using LinearAlgebra
using Serialization

struct eigen_structure
  H
  psi_R
  psi_L
  E_R
  E_L
end

mode="MPO"
#mode="sparse_tensor"

maxdim=20    

mu=1     # on-site dissipation
gamma=0   # skin term
Bz=0      # magnetic field

n=2
H=get_four_partite_spin_chain([0.5,1+gamma,1-gamma],mu,N=n,B=Bz,mode=mode)
N=n*4

psi0=0
if mode=="MPO"
  sites=firstsiteinds(H)
  psi0 = randomMPS(sites,1)
else
  dim=2^N
  psi0=sprand(dim,0.5)+im*sprand(dim,0.5)
end

K=careful_converge_krylov_MPO(H,psi0,which="SR",k=4,m=10,mode=mode,tol=1e-4,maxerror=1e-6,maxdim=maxdim,maxiter=50,dim_increase=10)
eigvecs=K.U
eigenvals=zeros(Complex,4)
psi=eigvecs[1]
B=K.B
for i=1:4
  eigenvals[i]=B[i,i]
end
println(eigenvals)

H2=get_H_dag_MPO(H)
K=careful_converge_krylov_MPO(H2,psi0,which="SR",k=4,m=10,mode=mode,tol=1e-4,maxerror=1e-6,maxdim=maxdim,maxiter=50,dim_increase=10)
psi_L=K.U[1]
E_L=K.B[1,1]

F=eigen_structure(H,psi,psi_L,eigenvals[1],E_L)
open(f->serialize(f,F),string("Eigen_L=",string(N),"_mu=",string(mu),"_gamma=",string(gamma),".jls"),"w")

