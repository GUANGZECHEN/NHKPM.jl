include("../../../src/NH_hamiltonian.jl")
include("../../../src/NHKPM.jl")
include("../../../src/plot_figs.jl")

symmetry="YSY"

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

N=20

H=get_NH_SSH(v1,v2,w1,w2,u,N)
#h=H.intra+H.tx+H.tmx  # for OBC use h=H.intra
h=H.intra

n_Ex=201
Ex_max=4
n_Ey=26
Ey_max=0.5

#Exs, Eys, rhos=plot_DOS_ED(h,n_Ex,Ex_max,n_Ey,Ey_max,1e-1)  # ED result
Exs, Eys, rhos=plot_DOS_KPM(h,n_Ex,Ex_max,n_Ey,Ey_max,E_max=10,npol=800)

open("DOS_NHSSH.OUT","w") do io
  writedlm(io, [Exs Eys rhos])
end

