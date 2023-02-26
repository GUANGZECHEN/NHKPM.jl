include("../../../src/NH_hamiltonian.jl")
include("../../../src/NHKPM.jl")
include("../../../src/plot_figs.jl")

symmetry="YSY"

t1=1/6   # t1=2.5: trivial; t1=1: topological, gapless; t1=1/6: topological, gapped
t2=1
gamma=4/3   

v1=t1-gamma/2
v2=t1+gamma/2
w1=t2
w2=t2

N=20 # number of sites

H=get_NH_SSH(v1,v2,w1,w2,0,N)
h_PBC=H.intra+H.tx+H.tmx  
h_OBC=H.intra

# tune the energy range and energy mesh
n_Ex=121       
Ex_max=1.2
n_Ey=81
Ey_max=0.8

#Exs, Eys, rhos=plot_DOS_ED(h,n_Ex,Ex_max,n_Ey,Ey_max,1e-1)  # ED result
Exs, Eys, rhos=plot_DOS_KPM(h_PBC,n_Ex,Ex_max,n_Ey,Ey_max,E_max=4,npol=300) # tune E_max according to spectrum to have smallest width possible
#Exs, Eys, rhos=plot_DOS_KPM(h_OBC,n_Ex,Ex_max,n_Ey,Ey_max,E_max=6,npol=50000) # for OBC

open("DOS_YSY.OUT","w") do io
  writedlm(io, [Exs Eys rhos])
end

