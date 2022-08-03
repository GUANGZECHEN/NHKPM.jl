include("../../../src/spin_hamiltonian.jl")
include("../../../src/NHKPM.jl")
include("../../../src/KrylovSchur.jl")

using LinearAlgebra

N=8

mode="sparse_tensor"   

dof_sites=ones(N)*3
dof_sites[1]=3
dof_sites[N]=3  

H=get_Haldane_chain(1,0.5,N)

H2=Matrix(H)
F=eigen(H2)
eigvalues,eigenvecs=F.values,F.vectors

println(eigvalues)

Exs=real(eigvalues)
Eys=imag(eigvalues)

plt.scatter(Exs,Eys)
plt.show()

energy=eigvalues[1]
psi=eigenvecs[:,1]
L_eigvecs=inv(eigenvecs)
psi2=conj(L_eigvecs[1,:])

vs=Array{Any}(undef,N)
vs_2=Array{Any}(undef,N)
for i=1:N
  Sp,Sm,Sz=get_op_S(i,dof_sites)
  vs[i]=Sz*psi
  vs_2[i]=Sz*psi2
end

nE=51
NN=nE*N

sites=zeros(NN)
Es=zeros(NN)
rhos=zeros(NN)

E_max=1
for i=0:nE-1
  E=real(E_max*i/(nE-1))
  for j=1:N
    ii=j+N*i
    println(ii)
    rho=get_rho_NH_ED(eigvalues,eigenvecs,L_eigvecs,E+energy,vs[j],vs_2[j],eta=5e-2)
    Es[ii]=E
    sites[ii]=j
    rhos[ii]=real(rho)
  end
end

open("LDOS_AKLT_ED.txt","w") do io
  writedlm(io, [sites Es rhos])
end


