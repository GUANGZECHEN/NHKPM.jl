include("hamiltonian.jl")
include("spin_hamiltonian.jl")
include("General_algorithms.jl")
using Arpack

function Cheb_poly(x,n::Int)
  return cos(n*acos(x))
end

function Cheb_poly_iter(x,N::Int)   # Chebyshev polynomial of order 0 to N-1, Tn[n]=T_(n-1)
  Tn=Array{Any,1}(undef,N)
  Tn[1]=1
  Tn[2]=x
  for n=3:N
    Tn[n]=2*x*Tn[n-1]-Tn[n-2]
  end  
  return Tn
end

function get_vn(H,v,N::Int;mode="matrix")  
  #print(H)
  vn=Array{Any}(undef,N)
  vn[1]=v
  vn[2]=general_product(H,v,mode=mode)
  for n=3:N
    f=general_product(H,vn[n-1],mode=mode)
    vn[n]=general_vec_sum(2*f,-vn[n-2],mode=mode)
  end
  return vn
end

function get_mu_n(H,v1,v2,N::Int;mode="matrix")  # Chebyshev moments for matrix H, of order 0 to N-1, mu_n[n]=mu_(n-1)
  vn=get_vn(H,v2,N,mode=mode)
  #print(vn,"\n")
  mu_n=zeros(Complex,N)
  for n=1:N
    mu_n[n]=general_inner_product(v1,vn[n])     
  end
  return mu_n
end

function get_mu_n_2(H,v,N::Int;mode="matrix")  # mu_n[n]=mu_(n-1), when v1=v2, can use this formula
  vn=get_vn(H,v,Int(N/2+1),mode=mode)
  #print(vn,"\n")
  #vT=conj(transpose(v))
  mu_n=zeros(N)
  mu_n[1]=general_inner_product(v,vn[1])
  mu_n[2]=general_inner_product(v,vn[2])
  for i=1:Int(N/2)-1
    f1=general_inner_product(vn[i+1],vn[i+1])
    f2=general_inner_product(vn[i+2],vn[i+2])
    mu_n[2*i+1]=general_vec_sum(2*f1,-mu_n[1],mode=mode)
    mu_n[2*i+2]=general_vec_sum(2*f2,-mu_n[2],mode=mode)      
  end
  return mu_n
end

function Jackson_kernel(N::Int)    # Jackson kernel of order 0 to N-1, gn[n]=g_(n-1)
  gn=zeros(N)
  q=pi/(N+1)
  for n=0:N-1
    gn[n+1]=((N-n+1)*cos(n*q)+sin(n*q)*cot(q))/(N+1)
  end
  return gn
end

function Lorentz_kernel(N::Int)   # Lorentz kernel of order 0 to N-1, gn[n]=g_(n-1)
  lambda=3
  gn=zeros(N)
  for n=0:N-1
    gn[n+1]=sinh(lambda*(1-n/N))/sinh(lambda)
  end
  return gn
end

function get_spec_kpm(E,mu_n,Tn,gn) # compute <v1|delta(H-E)|v2>  # mu_n are Chebyshev moments.
  N=size(mu_n)[1]

  rho=gn[1]*mu_n[1]*Tn[1]
  for n=2:N
    rho+=2*gn[n]*mu_n[n]*Tn[n]     
  end
  rho=rho*1/(pi*sqrt(1-E^2))
  #if real(rho)<0
  #  println(E)
  #end
  return real(rho)
end

function ldos_kpm(H,Es,ii;E_max=10,N=100,kernel="Jackson",mode="matrix")   # rho_ii(E), N is the number of moments
  #E_max=abs(eigs(H, which=:LM,nev=1)[1][1])
  H=H/E_max
  Es=Es/E_max
  dim=size(H)[1]
  v=zeros(dim)
  v[ii]=1
  
  mu_n=get_mu_n(H,v,v,N,mode=mode)
  
  gn=zeros(N)
  if kernel=="Jackson"
    gn=Jackson_kernel(N)
  else
    gn=Lorentz_kernel(N)
  end

  nE=size(Es)[1]
  rhos=zeros(Complex,nE)
  for i=1:nE
    Tn=Cheb_poly_iter(Es[i],N)
    rhos[i]=get_spec_kpm(Es[i],mu_n,Tn,gn)
  end
  return rhos
end

function dos_kpm(H,Es;E_max=10,N=100,kernel="Jackson",mode="matrix")
  dim=size(H)[1]

  nE=size(Es)[1]
  rhos=zeros(Complex,nE)
  for i=1:dim
    rhos+=ldos_kpm(H,Es,i,E_max=E_max,N=N,kernel=kernel,mode=mode)
  end
  return rhos
end

function get_dos_kpm(H::hamiltonian,E;E_max=10,nk=50,N=100,kernel="Jackson",mode="matrix")
  a1=H.g.inter_vector.x
  a2=H.g.inter_vector.y
  b1=2*pi/(a1[2]*a2[1]-a1[1]*a2[2])*[-a2[2],a2[1]]
  b2=2*pi/(a2[2]*a1[1]-a2[1]*a1[2])*[-a1[2],a1[1]]
  rho=0
  for i=1:nk
    for j=1:nk
      k=i/nk*b1+j/nk*b2
      Hk=get_Hk(H,k)
      rho+=dos_kpm(Hk,E,E_max=E_max,N=N,kernel=kernel,mode=mode)
    end
  end
  return rho
end

function get_ldos_kpm(H::hamiltonian,E,ii;E_max=10,nk=50,N=100,kernel="Jackson",mode="matrix")
  a1=H.g.inter_vector.x
  a2=H.g.inter_vector.y
  b1=2*pi/(a1[2]*a2[1]-a1[1]*a2[2])*[-a2[2],a2[1]]
  b2=2*pi/(a2[2]*a1[1]-a2[1]*a1[2])*[-a1[2],a1[1]]
  rho=0
  for i=1:nk
    for j=1:nk
      k=i/nk*b1+j/nk*b2
      Hk=get_Hk(H,k)
      rho+=ldos_kpm(H,E,ii,E_max=E_max,N=N,kernel=kernel,mode=mode)
    end
  end
  return rho
end

function ldos_ED(H,E,ii,eta)
  F=eigen(H)
  eigvals,eigvecs=real(F.values),F.vectors
  N=size(eigvals)[1]
  rho=0
  for i=1:N
    rho+=1/pi*eta/((E-eigvals[i])^2+eta^2)*abs(eigvecs[ii,i])^2     # the i th eigenvector is eigvecs[:,i]
  end
  return rho
end

function dos_ED(H,E,eta)
  F=eigen(H)
  eigvals,eigvecs=real(F.values),F.vectors
  N=size(eigvals)[1]
  rho=0
  for i=1:N
    rho+=1/pi*eta/((E-eigvals[i])^2+eta^2)
  end
  return rho
end
  
function get_E_max(H)
  return abs(eigs(H, which=:LM,nev=1)[1][1])
end 

function get_dos_ED(H::hamiltonian,E;nk=50,eta=1e-2)
  a1=H.g.inter_vector.x
  a2=H.g.inter_vector.y
  b1=2*pi/(a1[2]*a2[1]-a1[1]*a2[2])*[-a2[2],a2[1]]
  b2=2*pi/(a2[2]*a1[1]-a2[1]*a1[2])*[-a1[2],a1[1]]
  rho=0
  for i=1:nk
    for j=1:nk
      k=i/nk*b1+j/nk*b2
      Hk=get_Hk(H,k)
      rho+=dos_ED(Hk,E,eta)
    end
  end
  return rho
end

function get_ldos_ED(H::hamiltonian,E,ii,nk=50,eta=1e-2)
  a1=H.g.inter_vector.x
  a2=H.g.inter_vector.y
  b1=2*pi/(a1[2]*a2[1]-a1[1]*a2[2])*[-a2[2],a2[1]]
  b2=2*pi/(a2[2]*a1[1]-a2[1]*a1[2])*[-a1[2],a1[1]]
  rho=0
  for i=1:nk
    for j=1:nk
      k=i/nk*b1+j/nk*b2
      Hk=get_Hk(H,k)
      rho+=ldos_ED(Hk,E,eta)
    end
  end
  return rho
end

  





