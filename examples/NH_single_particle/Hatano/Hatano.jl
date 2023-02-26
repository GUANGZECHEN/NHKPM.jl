include("../../../src/NH_hamiltonian.jl")
include("../../../src/NHKPM.jl")
include("../../../src/plot_figs.jl")

symmetry="Hatano"

gamma=0.2   
t1=1-gamma
t2=1+gamma

N=20 # number of sites

H=get_NH_SSH(t1,t2,t1,t2,0,N)
h=H.intra+H.tx+H.tmx  # for OBC use h=H.intra
#h=H.intra

n_Ex=201
Ex_max=2.5
n_Ey=11
Ey_max=0.125

#Exs, Eys, rhos=plot_DOS_ED(h,n_Ex,Ex_max,n_Ey,Ey_max,1e-1)  # ED result
Exs, Eys, rhos=plot_DOS_KPM(h,n_Ex,Ex_max,n_Ey,Ey_max,E_max=5,npol=200)
#Exs, Eys, rhos=plot_DOS_KPM(h,n_Ex,Ex_max,n_Ey,Ey_max,E_max=5,npol=1200) # for OBC

open("DOS_Hatano.OUT","w") do io
  writedlm(io, [Exs Eys rhos])
end

