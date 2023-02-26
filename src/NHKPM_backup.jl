include("KPM.jl")

function get_vn_NH(H::Array,v::Array,N::Int)  # Chebyshev moments for matrix H(=H-z), vector v, of order 0 to 2*N-2, vn[n]=v_(n-1)
  #print(H)
  #print(N,"\n")
  dim=size(v)[1]
  vn=Array{Array,1}(undef,2*N)                # vn is dzdz*(H)|v>
  beta=Array{Array,1}(undef,2*N)             # beta is H|v>
  gamma=Array{Array,1}(undef,2*N)          # gamma is dz(H)|v>
  gamma2=Array{Array,1}(undef,2*N)         # gamma2 is dz*(H)|v>  

  beta[1]=v             # beta_0=v
  beta[2]=conj(transpose(H))*v

  gamma[1]=zeros(Complex,dim)
  gamma[2]=zeros(Complex,dim)

  gamma2[1]=zeros(Complex,dim)
  gamma2[2]=-v

  vn[1]=zeros(Complex,dim)
  vn[2]=zeros(Complex,dim) 

  for n=1:N-1                                
    beta[2*n+1]=2*H*beta[2*n]-beta[2*n-1]
    beta[2*n+2]=2*conj(transpose(H))*beta[2*n+1]-beta[2*n]

    gamma[2*n+1]=-2*beta[2*n]+2*H*gamma[2*n]-gamma[2*n-1]
    gamma[2*n+2]=2*conj(transpose(H))*gamma[2*n+1]-gamma[2*n]

    gamma2[2*n+1]=2*H*gamma2[2*n]-gamma2[2*n-1]
    gamma2[2*n+2]=-2*beta[2*n+1]+2*conj(transpose(H))*gamma2[2*n+1]-gamma2[2*n]

    vn[2*n+1]=-2*gamma2[2*n]+2*H*vn[2*n]-vn[2*n-1]
    vn[2*n+2]=-2*gamma[2*n+1]+2*conj(transpose(H))*vn[2*n+1]-vn[2*n]
  end
  return vn       # v includes both v_u and v_d
end

function get_mu_n_NH(H::Array,v::Array,N::Int)  # mu_n[n]=mu_2n
  vn=get_vn_NH(H,v,N)
  #print(size(vn)[1],"\n")
  #print(vn,"\n")
  vT=conj(transpose(v))
  mu_n=zeros(Complex,N-1)
  for n=1:N-1
    mu_n[n]=(vT*vn[2*n+1])     # vn[2n+1] is v_u
  end
  return mu_n
end

function ldos_kpm_NH(H,z,ii;E_max=10,N=100,kernel="Jackson")   # rho_ii(E), 2*N is the number of moments   # although tDOS is the same for top and bottom blocks (trace), lDOS is not, thus need to sum
  if N<3
    N=3
    print("N need to be at least 3")
  end
  #E_max=abs(eigs(H, which=:LM,nev=1)[1][1])
  H=(H-z*I)/E_max
  dim=size(H)[1]
  v=zeros(dim)
  v[ii]=1

  mu_n=get_mu_n_NH(H,v,N)
  #mu_n=get_mu_n_2(H,v,N)
  
  if kernel=="Jackson"
    gn=Jackson_kernel(2*N-1)
  else
    gn=Lorentz_kernel(2*N-1)
  end

  rho=0
  for n=1:N-1
    rho+=((-1)^n)*gn[2*n+1]*mu_n[n]/n     
  end
  rho=-rho*1/(pi*dim)
  return real(rho)
end

function dos_kpm_NH(H,z;E_max=10,N=100,kernel="Jackson")
  dim=size(H)[1]
  rho=0+0*im
  for i=1:dim
    rho+=ldos_kpm_NH(H,z,i,E_max=E_max,N=N,kernel=kernel)
  end
  return rho
end

#function get_dos_kpm_NH(H::hamiltonian,E;E_max=10,nk=50,N=100,kernel="Jackson")
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
    rho+=1/pi*eta/((abs(z-eigvals[i]))^2+eta^2)
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


#function get_dos_ED(H::hamiltonian,E;nk=50,eta=1e-2)
  #a1=H.g.inter_vector.x
  #a2=H.g.inter_vector.y
  #b1=2*pi/(a1[2]*a2[1]-a1[1]*a2[2])*[-a2[2],a2[1]]
  #b2=2*pi/(a2[2]*a1[1]-a2[1]*a1[2])*[-a1[2],a1[1]]
  #rho=0
  #for i=1:nk
    #for j=1:nk
      #k=i/nk*b1+j/nk*b2
     # Hk=get_Hk(H,k)
    #  rho+=dos_ED(Hk,E,eta)
   # end
  #end
 # return rho
#end

#function get_ldos_ED(H::hamiltonian,E,ii,nk=50,eta=1e-2)
 # a1=H.g.inter_vector.x
 # a2=H.g.inter_vector.y
 # b1=2*pi/(a1[2]*a2[1]-a1[1]*a2[2])*[-a2[2],a2[1]]
 ## b2=2*pi/(a2[2]*a1[1]-a2[1]*a1[2])*[-a1[2],a1[1]]
 # rho=0
 # for i=1:nk
 #   for j=1:nk
#      k=i/nk*b1+j/nk*b2
#      Hk=get_Hk(H,k)
#      rho+=ldos_ED(Hk,E,eta)
#    end
#  end
#  return rho
#end

  





