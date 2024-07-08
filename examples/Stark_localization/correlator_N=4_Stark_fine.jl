include("../NHKPM_source/special_spin_hamiltonian.jl")
include("../NHKPM_source/NHKPM.jl")
include("../NHKPM_source/KrylovSchur.jl")

using LinearAlgebra
using Arpack
using Random

using PyCall
plt = pyimport("matplotlib.pyplot")

using Serialization

struct H_and_random_vec
  H
  psi
end

function get_Liouvillian(Jx,Jy,B,gamma)
  N=4
  L=2*N
 
  H=get_Liouvillian_spin_chain_Stark(Jx,Jy,gamma,B,N=N,Bx=0,mode="MPO")
  sites=firstsiteinds(H)

  x=[ITensor() for i=1:2*N]
  M=[1 0;0 1]

  bonds=[Index(2,"link") for i=1:(2*N-1)]
  for i=1:(N-1)
    bonds[2*i]=Index(1,"link")
  end

  x[1]=ITensor(M,sites[1],bonds[1])
  x[2*N]=ITensor(M,bonds[2*N-1],sites[2*N])
  for i=2:2*N-1
    x[i]=ITensor(M,bonds[i-1],sites[i],bonds[i])
  end
  psi_1=MPS(x)

  Sz7=get_Sz_MPO(sites,7)
  psi_1=Sz7*psi_1

  return H, psi_1
end

function get_Liouvillian_DOS(i,H,psi_1,Jx,Jy,B,gamma)
  nEx=121
  nEy=821

  ii=i%nEx
  jj=Int((i-ii)/nEx)

  Ex=-0.2+0.6*(-1+(ii)/(nEx-1))
  Ey=0+4.1*(jj)/(nEy-1)
  z=Ex+im*Ey
  println(z)

  #println(dot(psi_1,psi_1))
  psi_2=psi_1 

  Delta=10
  npol=100*Delta
    
  DOS=get_spec_kpm_NH(H,z,psi_2,psi_1,E_max=Delta,N=npol,kernel="Jackson",mode="MPO",maxdim=16)     # Jx=0.25->E_max=8, Jx=0.75->E_max=9 
  DOS_re=real(DOS)
  DOS_im=imag(DOS)  

  open(string("correlator_fine_Stark_N=4_Sz7_Jx=",string(Jx),"_Jy=",string(Jy),"_B=",string(B),"_gamma=",string(gamma),"_index=",string(i),".OUT"),"w") do io
    writedlm(io, [Ex Ey DOS_re DOS_im])
  end
end


Jx=0.75
Jy=0.5
B=0.25
gamma=0.2
H, psi_1=get_Liouvillian(Jx,Jy,B,gamma)
i=parse(Int,ARGS[1])
for j=i:i+99
  get_Liouvillian_DOS(j,H,psi_1,Jx,Jy,B,gamma)
end


