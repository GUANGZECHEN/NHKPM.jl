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

# energy range
n_Ex=121       
Ex_max=1.2
n_Ey=81
Ey_max=0.8

plot_LDOS_KPM(h_OBC,n_Ex,Ex_max,n_Ey,Ey_max,E_max=4,npol=3000) # tune E_max to have smallest spreading of peak


