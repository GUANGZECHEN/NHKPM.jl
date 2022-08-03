include("../src/spin_hamiltonian.jl")
include("../src/NHKPM.jl")
include("../src/KrylovSchur.jl")

using LinearAlgebra

using PyCall
plt = pyimport("matplotlib.pyplot")

# LDOS on every site for complex E

N=8
#gamma=0
#Bz=-0.5
#delta=-0
#hx=0
#J=[1/2-gamma,1-gamma,1-gamma,1/2+gamma,1+gamma,1+gamma]
#B=[Bz+delta*im,hx/2,hx/2,Bz-delta*im,hx/2,hx/2]

mode="sparse_tensor"
#H=get_generalized_SSH_spin_chain(J,B,N=N,mode=mode)

Bz=0
mu=0
gamma=0.5
#H=get_four_partite_imag_onsite(H,mu,mode=mode)
H=get_four_partite_spin_chain([0.5,1+gamma,1-gamma],mu,B=Bz,N=2,mode=mode,BC="OBC")

psi0=0
if mode=="MPO"
  sites=firstsiteinds(H)
  psi0 = randomMPS(sites,10)
else
  dim=2^N
  psi0=sprand(dim,0.5)+im*sprand(dim,0.5)
end

psi=0
psi2=0
if mode=="MPO"    # compute also 2nd smallest, check orthogonality
  energy,psi,err=get_eigs_krylovMPO(H,psi0,which="SR",k=1,m=5,mode=mode,tol=1e-4,maxerror=1e-4,max_iter=10,maxdim=maxdim)
  energy=energy[1]
  psi=psi[1]
  H2=get_H_dag_MPO(H)
  energy_2,psi2,err_2=get_eigs_krylovMPO(H2,psi0,which="SR",k=1,m=5,mode=mode,tol=1e-4,maxerror=1e-4,max_iter=10,maxdim=maxdim)
  energy_2=energy_2[1]
  println(energy,energy_2)
  psi2=psi2[1]
else
  H2=Matrix(H)
  F=eigen(H2)
  eigvalues,eigenvecs=F.values,F.vectors
  println(F.values)
  local ii=argmin(real(eigvalues))
  println(ii)
  println(eigvalues[ii])
  psi=eigenvecs[:,ii]
  L_eigvecs=inv(eigenvecs)
  psi2=conj(L_eigvecs[ii,:])
  energy=eigvalues[ii]
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
    vs[i]=contract(Sp,psi,cutoff=1e-10,maxdim=maxdim)
    vs_2[i]=contract(Sp,psi2,cutoff=1e-10,maxdim=maxdim)
    vs_3[i]=contract(Sm,psi,cutoff=1e-10,maxdim=maxdim)
    vs_4[i]=contract(Sm,psi2,cutoff=1e-10,maxdim=maxdim)
    vs_5[i]=contract(Sz,psi,cutoff=1e-10,maxdim=maxdim)
    vs_6[i]=contract(Sz,psi2,cutoff=1e-10,maxdim=maxdim)
  else
    Sp,Sm,Sz=Si(i,N,S=1/2,mode=mode)
    factor=sqrt((1+gamma)/(1-gamma))
    Sx=(factor^(-i+1)*Sp+factor^(i-1)*Sm)/2
    Sx2=(factor^(-i+1)*Sm+factor^(i-1)*Sp)/2
    #Sy=(Sp-Sm)/(2*im)
    vs[i]=Sx*psi
    vs_2[i]=Sx2*psi2
    vs_3[i]=Sm*psi
    vs_4[i]=Sm*psi2
    vs_5[i]=Sz*psi
    vs_6[i]=Sz*psi2
    #println(dot(psi2,vs_5[i])*dot(vs_6[i],psi))
    #println(dot(vs_6[i],psi))
    #println(dot(psi2,vs[i])^2)
    #println(dot(vs_2[i],psi))
  end
end

nER=101
nEI=11
NN=nER*nEI*N

# all DOS
sites=zeros(NN)
Exs=zeros(NN)
Eys=zeros(NN)
rhos=zeros(NN)

# tDOS
tDOSs=zeros(nER*nEI)
tDOSs_imag=zeros(nER*nEI)
Exs_tDOS=zeros(nER*nEI)
Eys_tDOS=zeros(nER*nEI)

# LDOS_real
LDOS_real=zeros(nER*N)
LDOS_imag=zeros(nER*N)
Es_LDOS_real=zeros(nER*N)
sites_LDOS_real=zeros(nER*N)

ER_max=1
EI_max=0.05
for i=0:nER-1
  println(i)
  for j=0:nEI-1
#for i=50:50
#  for j=25:25
    Ex=ER_max*i/(nER-1)
    Ey=-EI_max+2*EI_max*j/(nEI-1)
    E=Ex+im*Ey
    local tDOS=0
    local tDOS2=0
    #for site_index=1:Int(N/2)
    for site_index=1:Int(N)
      ii=site_index+N*(j+i*nEI)
     # println(ii)
     # rho=abs(get_rho_NH_ED(eigvalues,eigenvecs,L_eigvecs,E+energy,vs[site_index],vs_2[site_index],eta=2e-2))
     # rho+=abs(get_rho_NH_ED(eigvalues,eigenvecs,L_eigvecs,E+energy,vs_3[site_index],vs_4[site_index],eta=2e-2))
     # rho+=abs(get_rho_NH_ED(eigvalues,eigenvecs,L_eigvecs,E+energy,vs_5[site_index],vs_6[site_index],eta=1e-2))
     # if real(rho)<0 || abs(imag(rho))>1e-5
     #   println(rho)
     # end
      rho=get_spec_kpm_NH(H,E+energy,vs_2[site_index],vs[site_index],E_max=5,N=200,mode=mode,kernel="Jackson")
     # rho+=get_spec_kpm_NH(H,E+energy,vs_4[site_index],vs_3[site_index],E_max=5,N=400,mode=mode,kernel="Jackson")
     # rho+=abs(get_spec_kpm_NH(H,E+energy,vs_6[site_index],vs_5[site_index],E_max=8,N=1200,mode=mode,kernel="Jackson"))
      Exs[ii]=Ex
      Eys[ii]=Ey
      sites[ii]=site_index
      rhos[ii]=real(rho)
      #println(rho)
      tDOS+=rho
      #tDOS2+=abs(rho)
      #println(ii)

      Es_LDOS_real[site_index+i*N]=Ex
      sites_LDOS_real[site_index+i*N]=site_index
      LDOS_real[site_index+i*N]+=real(rho)
      LDOS_imag[site_index+i*N]+=imag(rho)
      #rho=get_spec_kpm_NH(H,E+energy,vs_2[site_index],vs[site_index],E_max=10,N=300,mode=mode)
    end
    Exs_tDOS[j+i*nEI+1]=Ex
    Eys_tDOS[j+i*nEI+1]=Ey
    #println(tDOS)
    #println(tDOS2)
    tDOSs[j+i*nEI+1]=real(tDOS)
    tDOSs_imag[j+i*nEI+1]=imag(tDOS)
  end
end

open("all_DOS_NH_spin_chain.txt","w") do io
  writedlm(io, [Exs Eys sites rhos])
end

open("tDOS_NH_spin_chain.txt","w") do io
  writedlm(io, [Exs_tDOS Eys_tDOS tDOSs tDOSs_imag])
end

open("LDOS_real_projected_NH_spin_chain.txt","w") do io
  writedlm(io, [Es_LDOS_real sites_LDOS_real LDOS_real LDOS_imag])
end
