include("../../../src/NH_hamiltonian.jl")
include("../../../src/NHKPM.jl")
include("../../../src/plot_figs.jl")

symmetry="YSY"

t=1
u=1
gamma=0

N=20 # number of sites

H=get_four_partite_chain(t,u,n=Int(N/4),gamma=gamma)

h=H.intra+H.tx+H.tmx  
n_Ex=201
Ex_max=2.5
n_Ey=81
Ey_max=1

#Exs, Eys, rhos=plot_DOS_ED(h,n_Ex,Ex_max,n_Ey,Ey_max,1e-1)  # ED result
Exs, Eys, rhos=plot_DOS_KPM(h,n_Ex,Ex_max,n_Ey,Ey_max,E_max=6,npol=200)

open("DOS_Timo.OUT","w") do io
  writedlm(io, [Exs Eys rhos])
end

