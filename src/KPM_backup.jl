include("hamiltonian.jl")
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

function get_vn(H::Array,v::Array,N::Int)  # Chebyshev moments for matrix H, vector v, of order 0 to N-1, vn[n]=v_(n-1)
  #print(H)
  vn=Array{Array,1}(undef,N)
  vn[1]=v
  vn[2]=H*v
  for n=3:N
    vn[n]=2*H*vn[n-1]-vn[n-2]
  end
  return vn
end

function get_mu_n(H::Array,v::Array,N::Int)  # mu_n[n]=mu_(n-1)
  vn=get_vn(H,v,N)
  #print(vn,"\n")
  vT=conj(transpose(v))
  mu_n=zeros(Complex,N)
  for n=1:N
    mu_n[n]=(vT*vn[n])     
  end
  return mu_n
end

function get_mu_n_2(H::Array,v::Array,N::Int)  # mu_n[n]=mu_(n-1)
  vn=get_vn(H,v,Int(N/2+1))
  #print(vn,"\n")
  vT=conj(transpose(v))
  mu_n=zeros(N)
  mu_n[1]=(vT*vn[1])
  mu_n[2]=(vT*vn[2])
  for i=1:Int(N/2)-1
    mu_n[2*i+1]=2*(conj(transpose(vn[i+1]))*vn[i+1])-mu_n[1]
    mu_n[2*i+2]=2*(conj(transpose(vn[i+2]))*vn[i+1])-mu_n[2]      
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

function ldos_kpm(H,E,ii;E_max=10,N=100,kernel="Jackson")   # rho_ii(E), N is the number of moments
  if N<3
    N=3
    print("N need to be at least 3")
  end
  #E_max=abs(eigs(H, which=:LM,nev=1)[1][1])
  H=H/E_max
  E=E/E_max
  dim=size(H)[1]
  v=zeros(dim)
  v[ii]=1

  Tn=Cheb_poly_iter(E,N)
  mu_n=get_mu_n(H,v,N)
  #mu_n=get_mu_n_2(H,v,N)
  
  if kernel=="Jackson"
    gn=Jackson_kernel(N)
  else
    gn=Lorentz_kernel(N)
  end

  rho=gn[1]*mu_n[1]*Tn[1]
  for n=2:N
    rho+=2*gn[n]*mu_n[n]*Tn[n]     
  end
  rho=rho*1/(pi*sqrt(1-E^2))
  return real(rho)
end

function dos_kpm(H,E;E_max=10,N=100,kernel="Jackson")
  dim=size(H)[1]
  rho=0
  for i=1:dim
    rho+=ldos_kpm(H,E,i,E_max=E_max,N=N,kernel=kernel)
  end
  return rho
end

function get_dos_kpm(H::hamiltonian,E;E_max=10,nk=50,N=100,kernel="Jackson")
  a1=H.g.inter_vector.x
  a2=H.g.inter_vector.y
  b1=2*pi/(a1[2]*a2[1]-a1[1]*a2[2])*[-a2[2],a2[1]]
  b2=2*pi/(a2[2]*a1[1]-a2[1]*a1[2])*[-a1[2],a1[1]]
  rho=0
  for i=1:nk
    for j=1:nk
      k=i/nk*b1+j/nk*b2
      Hk=get_Hk(H,k)
      rho+=dos_kpm(Hk,E,E_max=E_max,N=N,kernel=kernel)
    end
  end
  return rho
end

function get_ldos_kpm(H::hamiltonian,E,ii;E_max=10,nk=50,N=100,kernel="Jackson")
  a1=H.g.inter_vector.x
  a2=H.g.inter_vector.y
  b1=2*pi/(a1[2]*a2[1]-a1[1]*a2[2])*[-a2[2],a2[1]]
  b2=2*pi/(a2[2]*a1[1]-a2[1]*a1[2])*[-a1[2],a1[1]]
  rho=0
  for i=1:nk
    for j=1:nk
      k=i/nk*b1+j/nk*b2
      Hk=get_Hk(H,k)
      rho+=ldos_kpm(H,E,ii,E_max=E_max,N=N,kernel=kernel)
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

  





