include("../../../src/special_spin_hamiltonian.jl")
include("../../../src/spin_hamiltonian.jl")
include("../../../src/NHKPM.jl")
include("../../../src/KrylovSchur.jl")

using LinearAlgebra
using Arpack
using Random

using PyCall
plt = pyimport("matplotlib.pyplot")


function get_Liouvillian_DOS(i)
  N=4
  L=2*N

  Jx=1
  Jy=0.5
  
  # energy mesh
  nEx=11
  nEy=11
  NN=nEx*nEy

  Exs=zeros(NN)
  Eys=zeros(NN)
  DOSs_KPM=zeros(NN)

  gamma=2*i/100
  gamma=round(gamma,digits=2)
  H=get_Liouvillian_spin_chain(Jx,Jy,gamma,N=N,Bx=0.0001,mode="MPO")  # Bx to let maximally mixed state to be GS
  
  sites=firstsiteinds(H)
    
  # construct state for the correlator
  psi_1=get_Liouvillian_GS(N)
  Sz1=get_Sz_MPO(sites,1)
  psi_1=Sz1*psi_1
  psi_2=psi_1    
    
  for ii=1:nEx
    Ex=0.5*(-1+(ii-1)/(nEx-1))
    for jj=1:nEy
      Ey=-0.5+0.5*(jj-1)/(nEy-1)
      z=Ex+im*Ey  
   
      DOS_KPM=get_spec_kpm_NH(H,z,psi_2,psi_1,E_max=10,N=300,kernel="Jackson",mode="MPO",maxdim=4)
      
      Exs[(ii-1)*nEy+jj]=Ex
      Eys[(ii-1)*nEy+jj]=Ey
      DOSs_KPM[(ii-1)*nEy+jj]+=abs(DOS_KPM)
    end
  end
  
  open(string("spec_Liouvillian_gamma=",string(gamma),".OUT"),"w") do io
    writedlm(io, [Exs Eys DOSs_KPM])
  end
end

for i=25:25
  get_Liouvillian_DOS(i)
  println(i)
end
