include("../../../src/spin_hamiltonian.jl")
include("../../../src/KPM.jl")
include("../../../src/KrylovSchur.jl")
include("../../../src/NH_hamiltonian.jl")

using LinearAlgebra

using PyCall
plt = pyimport("matplotlib.pyplot")

symmetry="Hatano"

if symmetry=="chiral_inv"
  t=pi/6
  v1=exp(im*pi/5)*sin(t)
  v2=exp(im*pi/5)*sin(t)
  w1=cos(t)
  w2=cos(t)
  u=0
elseif symmetry=="chiral"
  t=pi/6
  v1=sin(t)
  v2=exp(im*pi/10)*sin(t)
  w1=cos(t)
  w2=exp(im*pi/5)*cos(t)
  u=0  
elseif symmetry=="PT"
  t=pi/6
  v1=sin(t)
  v2=sin(t)
  w1=cos(t)
  w2=cos(t)
  u=0.3
elseif symmetry=="Hatano"
  gamma=0.2
  v1=1-gamma
  v2=1+gamma
  w1=1-gamma
  w2=1+gamma
  u=0
elseif symmetry=="YSY"
  v1=0.5
  v2=1.5
  w1=1
  w2=1
  u=0
end

mode="matrix"

dim=8

t2=0.2
H=get_NH_SSH(v1,v2,w1,w2,u,dim)
H=add_t2_1D(H,t2)
H=H.intra

N=dim
vs=Array{Any}(undef,N)
for i=1:N
  vs[i]=zeros(N)
  vs[i][i]=1
end

nE=201
NN=nE*N

sites=zeros(NN)
final_Es=zeros(NN)
final_rhos=zeros(NN)

E_max=4
Es=zeros(nE)
for i=1:nE
  Es[i]=E_max*(i-1)/(nE-1)
end

E0=6
H=H/E0
Es=Es/E0
energy=0
n_pol=500
gn=Jackson_kernel(n_pol)

for i=1:N
  mu_n=get_mu_n(H,vs[i],vs[i],n_pol,mode=mode)
  for j=1:nE
    ii=j+nE*(i-1)
    final_Es[ii]=Es[j]*E0
    E=Es[j]+energy
    Tn=Cheb_poly_iter(E,n_pol)
    final_rhos[ii]=get_spec_kpm(E,mu_n,Tn,gn)   # never forget GS energy
    #final_rhos[ii]=ldos_ED(H,E,i,1e-2)
    sites[ii]=i
    println(ii)
  end
end

open("LDOS_fermion_chain.txt","w") do io
  writedlm(io, [sites final_Es final_rhos])
end
