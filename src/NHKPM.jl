include("KPM.jl")
include("General_algorithms.jl")
include("spin_hamiltonian.jl")

function get_S_ent(psi,b)  # get entanglement entropy on bond b
  psi=psi/norm(psi)
  orthogonalize!(psi, b)
  U,S,V = svd(psi[b], (linkind(psi, b-1), siteind(psi,b)))
  SvN = 0.0
  for n=1:dim(S, 1)
    p = S[n,n]^2
    SvN -= p * log(p)
  end
  return SvN
end


function get_mu_n_NH_MPO(H,v_L,v,N::Int;maxdim=50,return_S=false)  # saves computational space, return_S returns the entanglement entropy of moments
  mu_n=zeros(Complex,N)
  vn=Array{Any}(undef,2)                # vn is dz*(H)|v>
  alpha=Array{Any}(undef,2)             # alpha is H|v>

  H_dag=get_H_dag_MPO(H)
  sites=firstsiteinds(H)
  println(size(sites)[1])
  #zero_vec=get_zero_MPS(sites)
  
  alpha[1]=contract(H_dag,v,cutoff=1e-10,maxdim=maxdim)                                                 
  alpha[2]=2*contract(H,alpha[1],cutoff=1e-10,maxdim=maxdim)-v

  vn[1]=v
  vn[2]=2*contract(H,vn[1],cutoff=1e-10,maxdim=maxdim)

  mu_n[1]=dot(v_L,vn[1])
  #println(get_S_ent(vn[1],Int(size(sites)[1]/2)))

  for n=2:(N-1)
    x=2*n-1
    y=2*n 
    alpha[1]=2*contract(H_dag,alpha[2],cutoff=1e-10,maxdim=maxdim)-alpha[1]
    vn[1]=2*alpha[2]+2*contract(H_dag,vn[2],cutoff=1e-10,maxdim=maxdim)-vn[1]
        
    alpha[2]=2*contract(H,alpha[1],cutoff=1e-10,maxdim=maxdim)-alpha[2]
    vn[2]=2*contract(H,vn[1],cutoff=1e-10,maxdim=maxdim)-vn[2]

    mu_n[n]=dot(v_L,vn[1])
    #println(get_S_ent(alpha[1],Int(size(sites)[1]/2)))
    #println(get_S_ent(alpha[2],Int(size(sites)[1]/2)))
    #println(get_S_ent(vn[1],Int(size(sites)[1]/2)))
    #println(get_S_ent(vn[2],Int(size(sites)[1]/2)))
  end

  return mu_n
end

function get_mu_n_NH(H,v1,v2,N::Int; mode="matrix",maxdim=50)  # mu_n[n]=v1*vn[2n-1], vn[n]=vn_n
  vn=Array{Any}(undef,2*N) 
  mu_n=zeros(Complex,N)
  if mode=="matrix" || mode=="sparse_tensor"
    vn=get_vn_NH(H,v2,N,mode=mode,maxdim=maxdim)
    for n=1:N
      mu_n[n]=general_inner_product(v1,vn[2*n-1],mode=mode)
    end
  else
    mu_n=get_mu_n_NH_MPO(H,v1,v2,N,maxdim=maxdim)
  end

  return mu_n
end

function get_vn_NH(H,v,N::Int;mode="matrix",maxdim=50)  # Chebyshev moments for matrix H(=z-H), vector v, of order 1 to 2*N, vn[n]=v_n   
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
  alpha[2]=general_vec_sum(2*general_product(H,alpha[1],mode=mode,maxdim=maxdim),-v,mode=mode,maxdim=maxdim)

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

function dos_kpm_NH(H,z;E_max=10,N=100,kernel="Jackson",mode="matrix")
  dim=size(H)[1]
  rho=0+0*im
  for i=1:dim
    rho+=ldos_kpm_NH(H,z,i,E_max=E_max,N=N,kernel=kernel,mode=mode)
  end
  
  return rho
end

function ldos_kpm_NH(H,z,ii;E_max=10,N=100,kernel="Jackson",mode="matrix")   # rho_ii(E), 2*N is the number of moments   # although tDOS is the same for top and bottom blocks (trace), lDOS is not, thus need to sum
  dim=size(H)[1]
  v=zeros(dim)
  v[ii]=1

  rho=get_spec_kpm_NH(H,z,v,v,E_max=E_max,N=N,kernel=kernel,mode=mode)
  return rho
end

function get_spec_kpm_NH(H,z,v1,v2;E_max=10,N=100,kernel="Jackson",mode="matrix",maxdim=50) # compute <v1|delta(H-E)|v2>
  if mode=="MPO"
    sites=firstsiteinds(H)
    Identity=get_Id_MPO(sites)
    H=+(-H,z*Identity,maxdim=maxdim)/E_max
    H=permute_inds_MPO(H)
  else
    H=(-H+z*I)/E_max
  end

  mu_n=get_mu_n_NH(H,v1,v2,N,mode=mode,maxdim=maxdim) 

  if kernel=="Jackson"
    gn=Jackson_kernel(2*N-1)
  elseif kernel=="Dirichlet"
    gn=ones(2*N-1)
  else
    gn=Lorentz_kernel(2*N-1)
  end

  rho=0
  for n=1:N
    rho+=gn[2*n-1]*mu_n[n]*sin((2*n-1)/2*pi)     
  end
  #rho=-rho*1/(pi)
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
  #println(dot(v_L,eigenvecs[:,1]))
  #println(dot(conj(L_eigvecs[1,:]),v_R))
  for i=1:N
    x=real(omega-eigvalues[i])
    y=imag(omega-eigvalues[i])
    f_L=dot(v_L,eigenvecs[:,i])
    f_R=dot(conj(L_eigvecs[i,:]),v_R)
    #rho+=eta/(x^2+eta^2)*(eta/(y^2+eta^2))*f_L*f_R
    rho+=eta/(x^2+y^2+eta^2)*f_L*f_R
  end
  rho=rho/pi
end


  





