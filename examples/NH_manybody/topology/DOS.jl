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

maxdim=16 # maximum bond dimension in MPS   

mu=0      # on-site dissipation
gamma=0   # skin term
Bz=0      # magnetic field
N=8       # number of sites

# load eigenvectors computed with Krylov-Schur
F=open(deserialize,string("Eigen_L=",string(N),"_mu=",string(mu),"_gamma=",string(gamma),".jls"))
H=F.H
psi=F.psi_R
psi2=F.psi_L
psi2=psi2/dot(psi,psi2)  # normalize
E_R,E_L=F.E_R,F.E_L

energy=real(E_R) #spectrum is real
sites_H=firstsiteinds(H)
N=size(sites_H)[1]

# energy range and number of points, for simplicity only spectrum on real axis is shown
nER=51
ER_max=1
Ey=0

sites=zeros(N*nER)
Exs=zeros(N*nER)
rhos=zeros(N*nER)

Delta=5 # scaling factor

for site_index=1:N
  println(site_index)
  Sx=get_Sx_MPO(sites_H,site_index)   
  v_R=contract(Sx,psi,cutoff=1e-10,maxdim=maxdim)
  v_L=contract(Sx,psi2,cutoff=1e-10,maxdim=maxdim)
    
  for i=1:nER
    E=ER_max*(i-1)/(nER-1)+im*Ey

    rho=get_spec_kpm_NH(H,E+energy,v_L,v_R,E_max=Delta,N=100,mode=mode,maxdim=maxdim,kernel="Jackson")   # since we have x-y rotational symmetry, Sx correlator = Sy correlator

    Exs[(site_index-1)*nER+i]=E
    sites[(site_index-1)*nER+i]=site_index
    rhos[(site_index-1)*nER+i]=abs(rho)
  end
end

open(string("DOS_L=",string(N),"_mu=",string(mu),"_gamma=",string(gamma),".txt"),"w") do io
  writedlm(io, [Exs sites rhos])
end






