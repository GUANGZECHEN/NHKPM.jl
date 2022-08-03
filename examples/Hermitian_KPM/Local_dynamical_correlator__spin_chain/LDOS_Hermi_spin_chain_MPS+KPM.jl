include("../../../src/spin_hamiltonian.jl")
include("../../../src/KPM.jl")
include("../../../src/KrylovSchur.jl")

using LinearAlgebra

using PyCall
plt = pyimport("matplotlib.pyplot")

N=8
gamma=0
Bz=0
hx=1
J=[1,0-gamma,0+gamma,1,0-gamma,0+gamma]    
B=[Bz,hx/2,hx/2,Bz,hx/2,hx/2]

mode="sparse_tensor"  # use "MPO" for MPS formalism
H=get_generalized_SSH_spin_chain(J,B,N=N,mode=mode)

psi0=0
if mode=="MPO"
  sites=firstsiteinds(H)
  psi0 = randomMPS(sites,10)
else
  dim=2^N
  psi0=sprand(dim,0.5)+im*sprand(dim,0.5)
end

#energy,psi,err=get_eigs_krylovMPO(H,psi0,which="SR",k=1,m=5,mode=mode,tol=1e-4,maxerror=1e-4,max_iter=10,maxdim=50)  # uncomment for MPS formalism
#energy=energy[1]
#psi=psi[1]

psi2=0
if mode=="sparse_tensor"
  H2=Matrix(H)
  F=eigen(H2)
  eigenvals,eigenvecs=real(F.values),F.vectors
  local ii=argmin(eigenvals)
  energy=eigenvals[ii]
  psi=eigenvecs[:,ii]
  println(norm(eigenvecs[:,ii]))
  L_eigvecs=inv(eigenvecs)
  psi2=conj(L_eigvecs[ii,:])
end

vs=Array{Any}(undef,N)
vs_2=Array{Any}(undef,N)
vs_3=Array{Any}(undef,N)
vs_4=Array{Any}(undef,N)
vs_5=Array{Any}(undef,N)
vs_6=Array{Any}(undef,N)
for i=1:N
  if mode=="MPO"
    Sp=get_Sp_MPO(sites,i)
    vs[i]=contract(Sp,psi,cutoff=1e-10,maxdim=50)
    #vs[i]=vs[i]/norm(vs[i])
  else
    Sp,Sm,Sz=Si(i,N,S=1/2,mode=mode)
    vs[i]=Sp*psi
    #vs[i]=multiply_JW_string(vs[i],i,mode=mode,N=N)      # apply JW string for fermion correlator
    vs_2[i]=Sp*psi2
    #vs_2[i]=multiply_JW_string(vs_2[i],i,mode=mode,N=N)
    vs_3[i]=Sm*psi
    vs_4[i]=Sm*psi2    
    vs_5[i]=Sz*psi
    vs_6[i]=Sz*psi2  
  end
end

nE=201
NN=nE*N

sites=zeros(NN)
final_Es=zeros(NN)
final_rhos=zeros(NN)

E_max=2
Es=zeros(nE)
for i=1:nE
  Es[i]=E_max*(i-1)/(nE-1)
end

E0=6
H=H/E0
Es=Es/E0
energy=energy/E0
println(energy)
n_pol=500
gn=Jackson_kernel(n_pol)

for i=1:N
  mu_n=get_mu_n(H,vs_2[i],vs[i],n_pol,mode=mode)
  mu_n_2=get_mu_n(H,vs_4[i],vs_3[i],n_pol,mode=mode)
  mu_n_3=get_mu_n(H,vs_6[i],vs_5[i],n_pol,mode=mode)
  for j=1:nE
    local ii=j+nE*(i-1)
    final_Es[ii]=Es[j]*E0
    E=Es[j]+energy
    Tn=Cheb_poly_iter(E,n_pol)
    final_rhos[ii]=get_spec_kpm(E,mu_n,Tn,gn)+get_spec_kpm(E,mu_n_2,Tn,gn)+get_spec_kpm(E,mu_n_3,Tn,gn)  # never forget GS energy
    sites[ii]=i
    println(ii)
  end
end

open("LDOS_Hermi_MPS_KPM.txt","w") do io
  writedlm(io, [sites final_Es final_rhos])
end
