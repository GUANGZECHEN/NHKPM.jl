include("KPM.jl")
include("General_algorithms.jl")
include("spin_hamiltonian.jl")

function get_vn_NH_old(H,v,N::Int;mode="matrix",maxdim=50)  # Chebyshev moments for matrix H(=H-z), vector v, of order 0 to 2*N-2, vn[n]=v_(n-1)   # maybe reduce number of operations to make better for MPS 
  vn=Array{Any}(undef,2*N)                # vn is dzdz*(H)|v>
  beta=Array{Any}(undef,2*N)             # beta is H|v>
  gamma=Array{Any}(undef,2*N)          # gamma is dz(H)|v>
  gamma2=Array{Any}(undef,2*N)         # gamma2 is dz*(H)|v>  

  H_dag=0
  zero_vec=0
  if mode=="MPO"
    H_dag=get_H_dag_MPO(H)
    sites=firstsiteinds(H)
    zero_vec=get_zero_MPS(sites)
  else
    dim=size(v)[1]
    H_dag=conj(transpose(H))
    zero_vec=zeros(Complex,dim)
    if mode=="sparse_tensor"
      zero_vec=spzeros(Complex,dim)
    end
  end

  beta[1]=v             # beta_0=v                                                # this part only computes the HH_dag contribution. For DOS it's fine since the H_dagH part has same contribution, but for LDOS,this is not fine
  beta[2]=general_product(H_dag,v,mode=mode,maxdim=maxdim)

  gamma[1]=zero_vec
  gamma[2]=zero_vec

  gamma2[1]=zero_vec
  gamma2[2]=-v

  vn[1]=zero_vec
  vn[2]=zero_vec 

  for n=1:N-1    
    x=2*n+1
    y=2*n+2                            
    beta[x]=general_vec_sum(2*general_product(H,beta[x-1],mode=mode,maxdim=maxdim),-beta[x-2],mode=mode,maxdim=maxdim)
    beta[y]=general_vec_sum(2*general_product(H_dag,beta[y-1],mode=mode,maxdim=maxdim),-beta[y-2],mode=mode,maxdim=maxdim)
 
    gamma[x]=general_vec_sum_2([-2*beta[x-1],2*general_product(H,gamma[x-1],mode=mode,maxdim=maxdim),-gamma[x-2]],mode=mode,maxdim=maxdim)
    gamma[y]=general_vec_sum(2*general_product(H_dag,gamma[y-1],mode=mode,maxdim=maxdim),-gamma[y-2],mode=mode,maxdim=maxdim)

    gamma2[x]=general_vec_sum(2*general_product(H,gamma2[x-1],mode=mode,maxdim=maxdim),-gamma2[x-2],mode=mode,maxdim=maxdim)
    gamma2[y]=general_vec_sum_2([-2*beta[y-1],2*general_product(H_dag,gamma2[y-1],mode=mode,maxdim=maxdim),-gamma2[y-2]],mode=mode,maxdim=maxdim)

    vn[x]=general_vec_sum_2([-2*gamma2[x-1],2*general_product(H,vn[x-1],mode=mode,maxdim=maxdim),-vn[x-2]],mode=mode,maxdim=maxdim)
    vn[y]=general_vec_sum_2([-2*gamma[y-1],2*general_product(H_dag,vn[y-1],mode=mode,maxdim=maxdim),-vn[y-2]],mode=mode,maxdim=maxdim)
  end

  #vn_2=Array{Any}(undef,2*N) 
  #beta[1]=v             # beta_0=v                                                # this part computes the H_dagH contribution.
  #beta[2]=general_product(H,v,mode=mode)

  #gamma[1]=zero_vec
  #gamma[2]=-v

  #gamma2[1]=zero_vec
  #gamma2[2]=zero_vec

  #vn_2[1]=zero_vec
  #vn_2[2]=zero_vec 

  #for n=1:N-1    
    #x=2*n+2
    #y=2*n+1
    #beta[y]=general_vec_sum(2*general_product(H_dag,beta[y-1],mode=mode),-beta[y-2],mode=mode)                            
    #beta[x]=general_vec_sum(2*general_product(H,beta[x-1],mode=mode),-beta[x-2],mode=mode)

   # gamma[y]=general_vec_sum(2*general_product(H_dag,gamma[y-1],mode=mode),-gamma[y-2],mode=mode)
   # gamma[x]=general_vec_sum_2([-2*beta[x-1],2*general_product(H,gamma[x-1],mode=mode),-gamma[x-2]],mode=mode)

   # gamma2[y]=general_vec_sum_2([-2*beta[y-1],2*general_product(H_dag,gamma2[y-1],mode=mode),-gamma2[y-2]],mode=mode)
   # gamma2[x]=general_vec_sum(2*general_product(H,gamma2[x-1],mode=mode),-gamma2[x-2],mode=mode)

  #  vn_2[y]=general_vec_sum_2([-2*gamma[y-1],2*general_product(H_dag,vn[y-1],mode=mode),-vn[y-2]],mode=mode)
  #  vn_2[x]=general_vec_sum_2([-2*gamma2[x-1],2*general_product(H,vn[x-1],mode=mode),-vn[x-2]],mode=mode)
  #end

  return vn#,vn_2       # v includes both v_u and v_d
end

function get_vn_NH_MPO_old(H,v,N::Int;maxdim=50)  # Chebyshev moments for matrix H(=H-z), vector v, of order 0 to 2*N-2, vn[n]=v_(n-1)   # maybe reduce number of operations to make better for MPS 
  vn=Array{Any}(undef,2*N)                # vn is dzdz*(H)|v>
  beta=Array{Any}(undef,2*N)             # beta is H|v>
  gamma=Array{Any}(undef,2*N)          # gamma is dz(H)|v>
  gamma2=Array{Any}(undef,2*N)         # gamma2 is dz*(H)|v>  

  H_dag=get_H_dag_MPO(H)
  sites=firstsiteinds(H)
  zero_vec=get_zero_MPS(sites)

  beta[1]=v             # beta_0=v                                                # this part only computes the HH_dag contribution. For DOS it's fine since the H_dagH part has same contribution, but for LDOS,this is not fine
  beta[2]=contract(H_dag,v,cutoff=1e-10,maxdim=maxdim)
  beta[3]=2*contract(H,beta[2],cutoff=1e-10,maxdim=maxdim)-beta[1]
  beta[4]=2*contract(H_dag,beta[3],cutoff=1e-10,maxdim=maxdim)-beta[2]

  gamma[1]=zero_vec
  gamma[2]=zero_vec
  gamma[3]=-2*beta[2]
  gamma[4]=2*contract(H_dag,gamma[3],cutoff=1e-10,maxdim=maxdim)

  gamma2[1]=zero_vec
  gamma2[2]=-v
  gamma2[3]=2*contract(H,gamma2[2],cutoff=1e-10,maxdim=maxdim)
  gamma2[4]=-2*beta[3]+2*contract(H_dag,gamma2[3],cutoff=1e-10,maxdim=maxdim)-gamma2[2]

  vn[1]=zero_vec
  vn[2]=zero_vec 
  vn[3]=-2*gamma2[2]
  vn[4]=-2*gamma[3]+2*contract(H_dag,vn[3],cutoff=1e-10,maxdim=maxdim)

  for n=2:N-1    
    x=2*n+1
    y=2*n+2                            
    beta[x]=2*contract(H,beta[x-1],cutoff=1e-10,maxdim=maxdim)-beta[x-2]  # beta_odd=no H
    beta[y]=2*contract(H_dag,beta[y-1],cutoff=1e-10,maxdim=maxdim)-beta[y-2] # beta_even= one H
 
    gamma[x]=-2*beta[x-1]+2*contract(H,gamma[x-1],cutoff=1e-10,maxdim=maxdim)-gamma[x-2]  # H
    gamma[y]=2*contract(H_dag,gamma[y-1],cutoff=1e-10,maxdim=maxdim)-gamma[y-2]           # no H

    gamma2[x]=2*contract(H,gamma2[x-1],cutoff=1e-10,maxdim=maxdim)-gamma2[x-2]  # H
    gamma2[y]=-2*beta[y-1]+2*contract(H_dag,gamma2[y-1],cutoff=1e-10,maxdim=maxdim)-gamma2[y-2] # no H

    vn[x]=-2*gamma2[x-1]+2*contract(H,vn[x-1],cutoff=1e-10,maxdim=maxdim)-vn[x-2]  # no H
    vn[y]=-2*gamma[y-1]+2*contract(H_dag,vn[y-1],cutoff=1e-10,maxdim=maxdim)-vn[y-2] # H
  end

  return vn
end

function get_mu_n_NH_old(H,v1,v2,N::Int; mode="matrix",maxdim=50)  # mu_n[n]=mu_2n
  vn=Array{Any}(undef,2*N) 
  mu_n=zeros(Complex,N-1)
  if mode=="matrix" || mode=="sparse_tensor"
    vn=get_vn_NH_old(H,v2,N,mode=mode,maxdim=maxdim)
    for n=1:N-1
      mu_n[n]=general_inner_product(v1,vn[2*n+1],mode=mode)#+0*general_inner_product(v1,vn_2[2*n+1],mode=mode)    # vn[2n+1] is v_u
    end
  else
    mu_n=get_mu_n_NH_MPO_old(H,v1,v2,N,maxdim=maxdim)
  end

  return mu_n
end

function get_mu_n_NH_MPO_old(H,v_L,v,N::Int;maxdim=50)  # saves computational space
  mu_n=zeros(Complex,N-1)
  vn=Array{Any}(undef,4)                # vn is dzdz*(H)|v>
  beta=Array{Any}(undef,4)             # beta is H|v>
  gamma=Array{Any}(undef,4)          # gamma is dz(H)|v>
  gamma2=Array{Any}(undef,4)         # gamma2 is dz*(H)|v>

  H_dag=get_H_dag_MPO(H)
  sites=firstsiteinds(H)
  zero_vec=get_zero_MPS(sites)

  beta[1]=v             # beta_0=v                                                # this part only computes the HH_dag contribution. For DOS it's fine since the H_dagH part has same contribution, but for LDOS,this is not fine
  beta[2]=contract(H_dag,v,cutoff=1e-10,maxdim=maxdim)
  beta[3]=2*contract(H,beta[2],cutoff=1e-10,maxdim=maxdim)-beta[1]
  beta[4]=2*contract(H_dag,beta[3],cutoff=1e-10,maxdim=maxdim)-beta[2]

  gamma[1]=zero_vec
  gamma[2]=zero_vec
  gamma[3]=-2*beta[2]
  gamma[4]=2*contract(H_dag,gamma[3],cutoff=1e-10,maxdim=maxdim)

  gamma2[1]=zero_vec
  gamma2[2]=-v
  gamma2[3]=2*contract(H,gamma2[2],cutoff=1e-10,maxdim=maxdim)
  gamma2[4]=-2*beta[3]+2*contract(H_dag,gamma2[3],cutoff=1e-10,maxdim=maxdim)-gamma2[2]

  vn[1]=zero_vec
  vn[2]=zero_vec
  vn[3]=-2*gamma2[2]
  vn[4]=-2*gamma[3]+2*contract(H_dag,vn[3],cutoff=1e-10,maxdim=maxdim)

  mu_n[1]=dot(v_L,vn[3])

  for n=2:N-1
    x=2*n+1
    y=2*n+2
    beta[1]=2*contract(H,beta[4],cutoff=1e-10,maxdim=maxdim)-beta[3]  # beta_odd=no H
    beta[2]=2*contract(H_dag,beta[1],cutoff=1e-10,maxdim=maxdim)-beta[4] # beta_even= one H

    gamma[1]=-2*beta[4]+2*contract(H,gamma[4],cutoff=1e-10,maxdim=maxdim)-gamma[3]  # H
    gamma[2]=2*contract(H_dag,gamma[1],cutoff=1e-10,maxdim=maxdim)-gamma[4]           # no H

    gamma2[1]=2*contract(H,gamma2[4],cutoff=1e-10,maxdim=maxdim)-gamma2[3]  # H
    gamma2[2]=-2*beta[1]+2*contract(H_dag,gamma2[1],cutoff=1e-10,maxdim=maxdim)-gamma2[4] # no H

    vn[1]=-2*gamma2[4]+2*contract(H,vn[4],cutoff=1e-10,maxdim=maxdim)-vn[3]  # no H
    vn[2]=-2*gamma[1]+2*contract(H_dag,vn[1],cutoff=1e-10,maxdim=maxdim)-vn[4] # H

    mu_n[n]=dot(v_L,vn[1])

    beta[3]=beta[1]
    beta[4]=beta[2]
    gamma[3]=gamma[1]
    gamma[4]=gamma[2]
    gamma2[3]=gamma2[1]
    gamma2[4]=gamma2[2]
    vn[3]=vn[1]
    vn[4]=vn[2]
  end

  return mu_n
end


function get_spec_kpm_NH_old(H,z,v1,v2;E_max=10,N=100,kernel="Jackson",mode="matrix",maxdim=50) # compute <v1|delta(H-E)|v2>
  if mode=="MPO"
    sites=firstsiteinds(H)
    Identity=get_Id_MPO(sites)
    H=+(H,-z*Identity,maxdim=maxdim)/E_max
    H=permute_inds_MPO(H)
  else
    H=(H-z*I)/E_max
  end

  mu_n=get_mu_n_NH_old(H,v1,v2,N,mode=mode,maxdim=maxdim) 

  if kernel=="Jackson"
    gn=Jackson_kernel(2*N-1)
  else
    gn=Lorentz_kernel(2*N-1)
  end

  rho=0
  for n=1:N-1
    rho+=((-1)^n)*gn[2*n+1]*mu_n[n]/n     
  end
  rho=-rho*1/(pi)
  return abs(rho)
end

function ldos_kpm_NH_old(H,z,ii;E_max=10,N=100,kernel="Jackson",mode="matrix")   # rho_ii(E), 2*N is the number of moments   # although tDOS is the same for top and bottom blocks (trace), lDOS is not, thus need to sum
  dim=size(H)[1]
  v=zeros(dim)
  v[ii]=1

  rho=get_spec_kpm_NH_old(H,z,v,v,E_max=E_max,N=N,kernel=kernel,mode=mode)
  return rho
end

function dos_kpm_NH_old(H,z;E_max=10,N=100,kernel="Jackson",mode="matrix")
  dim=size(H)[1]
  rho=0+0*im
  for i=1:dim
    rho+=ldos_kpm_NH(H,z,i,E_max=E_max,N=N,kernel=kernel,mode=mode)
  end
  return rho
end






function get_spec_kpm_NH_new(H,z,v1,v2;E_max=10,N=100,kernel="Jackson",mode="matrix",maxdim=50) # compute <v1|delta(H-E)|v2>
  if mode=="MPO"
    sites=firstsiteinds(H)
    Identity=get_Id_MPO(sites)
    H=+(-H,z*Identity,maxdim=maxdim)/E_max
    H=permute_inds_MPO(H)
  else
    H=(-H+z*I)/E_max
  end

  mu_n=get_mu_n_NH_new(H,v1,v2,N,mode=mode,maxdim=maxdim) 

  if kernel=="Jackson"
    gn=Jackson_kernel(2*N-1)
  else
    gn=Lorentz_kernel(2*N-1)
  end

  rho=0
  for n=1:N
    rho+=gn[2*n-1]*mu_n[n]*sin((2*n-1)/2*pi)     
  end
  #rho=-rho*1/(pi)
  return real(rho)
end

function get_mu_n_NH_MPO_new(H,v_L,v,N::Int;maxdim=50)  # saves computational space
  mu_n=zeros(Complex,N)
  vn=Array{Any}(undef,2)                # vn is dz*(H)|v>
  alpha=Array{Any}(undef,2)             # alpha is H|v>

  H_dag=get_H_dag_MPO(H)
  sites=firstsiteinds(H)
  #zero_vec=get_zero_MPS(sites)
  
  alpha[1]=contract(H_dag,v,cutoff=1e-10,maxdim=maxdim)                                                 
  alpha[2]=2*contract(H,alpha[1],cutoff=1e-10,maxdim=maxdim)-v

  vn[1]=v
  vn[2]=2*contract(H,vn[1],cutoff=1e-10,maxdim=maxdim)

  mu_n[1]=dot(v_L,vn[1])

  for n=2:N-1
    x=2*n-1
    y=2*n 
    alpha[1]=2*contract(H_dag,alpha[2],cutoff=1e-10,maxdim=maxdim)-alpha[1]
    vn[1]=2*alpha[2]+2*contract(H_dag,vn[2],cutoff=1e-10,maxdim=maxdim)-vn[1]
        
    alpha[2]=2*contract(H,alpha[1],mode=mode,maxdim=maxdim)-alpha[2]
    vn[2]=2*contract(H,vn[1],mode=mode,maxdim=maxdim)-vn[2]

    mu_n[n]=dot(v_L,vn[1])
  end

  return mu_n
end

function get_mu_n_NH_new(H,v1,v2,N::Int; mode="matrix",maxdim=50)  # mu_n[n]=v1*vn[2n-1], vn[n]=vn_n
  vn=Array{Any}(undef,2*N) 
  mu_n=zeros(Complex,N)
  if mode=="matrix" || mode=="sparse_tensor"
    vn=get_vn_NH_new(H,v2,N,mode=mode,maxdim=maxdim)
    for n=1:N
      mu_n[n]=general_inner_product(v1,vn[2*n-1],mode=mode)
    end
  else
    mu_n=get_mu_n_NH_MPO_new(H,v1,v2,N,maxdim=maxdim)
  end

  return mu_n
end

function get_vn_NH_new(H,v,N::Int;mode="matrix",maxdim=50)  # Chebyshev moments for matrix H(=z-H), vector v, of order 1 to 2*N, vn[n]=v_n   
  vn=Array{Any}(undef,2*N)                # vn is dz*(H)|v>
  alpha=Array{Any}(undef,2*N)             # alpha is H|v>

  H_dag=0
  zero_vec=0
  if mode=="MPO"
    H_dag=get_H_dag_MPO(H)
    sites=firstsiteinds(H)
    zero_vec=get_zero_MPS(sites)
  else
    dim=size(v)[1]
    H_dag=conj(transpose(H))
    zero_vec=zeros(Complex,dim)
    if mode=="sparse_tensor"
      zero_vec=spzeros(Complex,dim)
    end
  end

  alpha[1]=general_product(H_dag,v,mode=mode,maxdim=maxdim)                                                 
  alpha[2]=2*general_vec_sum(2*general_product(H,alpha[1],mode=mode,maxdim=maxdim),-v,mode=mode,maxdim=maxdim)

  vn[1]=v
  vn[2]=2*general_product(H,vn[1],mode=mode,maxdim=maxdim)

  for n=2:N    
    x=2*n-1
    y=2*n                            
    
    alpha[x]=general_vec_sum(2*general_product(H_dag,alpha[x-1],mode=mode,maxdim=maxdim),-alpha[x-2],mode=mode,maxdim=maxdim)
    alpha[y]=general_vec_sum(2*general_product(H,alpha[y-1],mode=mode,maxdim=maxdim),-alpha[y-2],mode=mode,maxdim=maxdim)

    vn[x]=general_vec_sum_2([2*alpha[x-1],2*general_product(H_dag,vn[x-1],mode=mode,maxdim=maxdim),-vn[x-2]],mode=mode,maxdim=maxdim)
    vn[y]=general_vec_sum(2*general_product(H,vn[y-1],mode=mode,maxdim=maxdim),-vn[y-2],mode=mode,maxdim=maxdim)
  end

  return vn
end

function dos_kpm_NH_new(H,z;E_max=10,N=100,kernel="Jackson",mode="matrix")
  dim=size(H)[1]
  rho=0+0*im
  for i=1:dim
    rho+=ldos_kpm_NH_new(H,z,i,E_max=E_max,N=N,kernel=kernel,mode=mode)
  end
  return rho
end

function ldos_kpm_NH_new(H,z,ii;E_max=10,N=100,kernel="Jackson",mode="matrix")   # rho_ii(E), 2*N is the number of moments   # although tDOS is the same for top and bottom blocks (trace), lDOS is not, thus need to sum
  dim=size(H)[1]
  v=zeros(dim)
  v[ii]=1

  rho=get_spec_kpm_NH_new(H,z,v,v,E_max=E_max,N=N,kernel=kernel,mode=mode)
  return rho
end

#function get_dos_kpm_NH(H::hamiltonian,E;E_max=10,nk=50,N=100,kernel="Jackson")   # in NH case there is no Bloch wave
#  a1=H.g.inter_vector.x
#  a2=H.g.inter_vector.y
#  b1=2*pi/(a1[2]*a2[1]-a1[1]*a2[2])*[-a2[2],a2[1]]
#  b2=2*pi/(a2[2]*a1[1]-a2[1]*a1[2])*[-a1[2],a1[1]]
#  rho=0
#  for i=1:nk
#    for j=1:nk
#      k=i/nk*b1+j/nk*b2
#      Hk=get_Hk(H,k)
#      rho+=dos_kpm(Hk,E,E_max=E_max,N=N,kernel=kernel)
#    end
#  end
#  return rho
#end

#function get_ldos_kpm_NH(H::hamiltonian,E,ii;E_max=10,nk=50,N=100,kernel="Jackson")
#  a1=H.g.inter_vector.x
#  a2=H.g.inter_vector.y
#  b1=2*pi/(a1[2]*a2[1]-a1[1]*a2[2])*[-a2[2],a2[1]]
#  b2=2*pi/(a2[2]*a1[1]-a2[1]*a1[2])*[-a1[2],a1[1]]
#  rho=0
#  for i=1:nk
#    for j=1:nk
#      k=i/nk*b1+j/nk*b2
#      Hk=get_Hk(H,k)
#      rho+=ldos_kpm(H,E,ii,E_max=E_max,N=N,kernel=kernel)
#    end
#  end
#  return rho
#end

function dos_ED_NH(H,z,eta)
  F=eigen(H)
  eigvals=F.values
  N=size(eigvals)[1]
  rho=0
  for i=1:N
    rho+=1/pi*eta/((real(z-eigvals[i]))^2+eta^2)*(1/pi*eta/((imag(z-eigvals[i]))^2+eta^2))
  end
  return real(rho)
end

function ldos_ED_NH(H,z,ii,eta)
  F=eigen(H)
  eigvals=F.values
  eigvecs_r=F.vectors    # F.vectors[:,i]=|v[i]>, right eigvec
  eigvecs_l=inv(eigvecs_r)  # eigvecs_l[i,:]=<u[i]|, left eigvec
  N=size(eigvals)[1]
  rho=0
  for i=1:N
    v=eigvecs_r[:,i]
    u=eigvecs_l[i,:]    
    rho+=v[ii]*u[ii]/pi*eta/((abs(z-eigvals[i]))^2+eta^2)
  end
  return real(rho)
end
  
function get_E_max(H)
  return abs(eigs(H, which=:LM,nev=1)[1][1])
end 

function dos_Green_NH(H,z,eta)     # this is wrong           
  G=inv(H-z*I)     # Green's function
  dim=size(G,1)
  rho=0
  for i=1:dim
    rho+=1/pi*eta/((1/abs(G[i,i]))^2+eta^2)
  end
  return real(rho)
end

function get_rho_NH_ED(eigvalues,eigenvecs,L_eigvecs,omega,v_R,v_L;eta=1e-2)
  N=size(eigvalues)[1]
  rho=0
  for i=1:N
    x=real(omega-eigvalues[i])
    y=imag(omega-eigvalues[i])
    f_L=dot(v_L,eigenvecs[:,i])
    f_R=dot(v_R,L_eigvecs[i,:])
    rho+=eta/(x^2+eta^2)*(eta/(y^2+eta^2))*f_L*f_R
  end
  rho=rho/pi/pi
end


  





