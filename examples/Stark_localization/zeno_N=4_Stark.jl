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

function get_random_pauli_chain(l,L,psi_0)
  Sp=Array{Any}(undef,L)
  Sm=Array{Any}(undef,L)
  Sz=Array{Any}(undef,L)
  
  for ii=1:L
    Sp[ii],Sm[ii],Sz[ii]=Si(ii,L)
  end  
      
  rand_index=rand(0:3,l)
  println(rand_index)

  psi_1=psi_0  
  for jj=1:l
    if rand_index[jj]==1
      psi_1=Sp[jj]*psi_1
    elseif rand_index[jj]==2
      psi_1=Sm[jj]*psi_1
    elseif rand_index[jj]==3
      psi_1=Sz[jj]*psi_1
    end
  end 
  return psi_1
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

  Sz4=get_Sz_MPO(sites,7)
  psi_1=Sz4*psi_1

  return H, psi_1
end

function get_Liouvillian_DOS(i,H,psi_1,gamma,Ex)
  Ey=2.25*(i-1)/(151-1)
  z=Ex+im*Ey
  println(z)

  #println(dot(psi_1,psi_1))
  psi_2=psi_1 

  Delta=25
  npol=30*Delta
    
  DOS=get_spec_kpm_NH(H,z,psi_2,psi_1,E_max=Delta,N=npol,kernel="Jackson",mode="MPO",maxdim=16)     # Jx=0.25->E_max=20, Jx=0.75->E_max=25
 
  return real(DOS), imag(DOS)
end

Jx=0.75
Jy=0.5
B=0.25

nEx=101
n_gamma=51
i=parse(Int,ARGS[1])
ii=i%nEx
jj=Int((i-ii)/nEx)

Ex=1.5*(-1+(ii)/(nEx-1))
gamma=1*jj/(n_gamma-1)

Ex=round(Ex,digits=3)
gamma=round(gamma,digits=2)

H, psi_1=get_Liouvillian(Jx,Jy,B,gamma)
DOSs=zeros(151)
DOSs_im=zeros(151)
DOSs=zeros(151)
for j=1:151
  dos_re, dos_im=get_Liouvillian_DOS(j,H,psi_1,gamma,Ex)
  DOSs[j]=dos_re
  DOSs_im[j]=dos_im
  DOSs[j]=sqrt(dos_re^2+dos_im^2)
end

DOS=sum(DOSs)
open(string("Zeno_N=4_Stark_Jx=",string(Jx),"_Jy=",string(Jy),"_gamma=",string(gamma),"_Ex=",string(Ex),"B=",string(B),".OUT"),"w") do io
  writedlm(io, [gamma Ex DOS])
end

open(string("Zeno_N=4_Stark_origin_Jx=",string(Jx),"_Jy=",string(Jy),"_gamma=",string(gamma),"_Ex=",string(Ex),"B=",string(B),".OUT"),"w") do io
  writedlm(io, [DOSs DOSs_im])
end

