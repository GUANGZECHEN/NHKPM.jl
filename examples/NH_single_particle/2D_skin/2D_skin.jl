include("../../../src/NH_hamiltonian.jl")
include("../../../src/NH_hamiltonian_2D.jl")
include("../../../src/NHKPM.jl")
include("../../../src/plot_figs.jl")

t1=1
t2=1
u=0.0
gamma=0.1
L=8

H=get_2D_square_NH(t1,t2,u,n=L,gamma=gamma)
h_OBC=H.intra
h_PBC=H.intra+H.tx+H.tmx+H.ty+H.tmy

n_Ex=201
Ex_max=4.5
#n_Ey=21
#Ey_max=0.45

n_Ey=5
Ey_max=0.09

#Exs, Eys, rhos=plot_DOS_ED(h,n_Ex,Ex_max,n_Ey,Ey_max,1e-1)  # ED result
#Exs, Eys, rhos=plot_DOS_KPM(h_PBC,n_Ex,Ex_max,n_Ey,Ey_max,E_max=10,npol=300)
Exs, Eys, rhos=plot_DOS_KPM(h_OBC,n_Ex,Ex_max,n_Ey,Ey_max,E_max=10,npol=800) # for OBC

open("DOS_2D_skin.OUT","w") do io
  writedlm(io, [Exs Eys rhos])
end

