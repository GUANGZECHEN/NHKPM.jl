include("../../src/NH_hamiltonian.jl")
include("../../src/NHKPM.jl")

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
elseif symmetry=="Hatano"
  gamma=0.2
  v1=1-gamma
  v2=1+gamma
  w1=1-gamma
  w2=1+gamma
  u=0
elseif symmetry=="PT"
  t=pi/6
  v1=sin(t)
  v2=sin(t)
  w1=cos(t)
  w2=cos(t)
  u=0.3
elseif symmetry=="YSY"
  v1=0.5
  v2=1.5
  w1=1
  w2=1
  u=0
end

function get_Hermitrized_DOS(h)    # rho(0) for Hermitrized h-z
  n_E=100
  E_max=1.5
  N=(n_E-1)^2
  Exs=zeros(N)
  Eys=zeros(N)   # E=Ex+iEy
  rhos=zeros(N)
  for i=1:n_E-1
    for j=1:n_E-1
      ii = Int(j+(i-1)*(n_E-1))
      Exs[ii]=-E_max+2*E_max*i/n_E
      Eys[ii]=-E_max+2*E_max*j/n_E
      z=Exs[ii]+im*Eys[ii]
      H=Hermitrize(h-z*I)
      rhos[ii]=dos_ED(H,0,1e-2)
    end
  end
  open("rho0_Hermitrized_H.txt","w") do io
    writedlm(io, [Exs Eys rhos])
  end
end

N=8

H=get_NH_SSH(v1,v2,w1,w2,u,N)
h=H.intra+H.tx+H.tmx  # for OBC use h=H.intra

n_E=100
E_max=2

N=(n_E-1)^2
Exs=zeros(N)
Eys=zeros(N)   # E=Ex+iEy
rho1s=zeros(N)
rho2s=zeros(N)

for i=1:n_E-1
  for j=1:n_E-1
    ii = Int(j+(i-1)*(n_E-1))
    Exs[ii]=-E_max+2*E_max*i/n_E
    Eys[ii]=-E_max+2*E_max*j/n_E
    z=Exs[ii]+im*Eys[ii]
    rho1s[ii]=dos_ED_NH(h,z,2e-2)
    rho2s[ii]=real(dos_kpm_NH(h,z,E_max=10,N=400,kernel="Jackson"))
  end
end

open("DOS_NH_SSH.txt","w") do io
  writedlm(io, [Exs Eys rho1s rho2s])
end
