include("../../../src/NH_hamiltonian.jl")
include("../../../src/NHKPM.jl")
include("../../../src/plot_figs.jl")

symmetry="Hatano"

gamma=0.2   
t1=1-gamma
t2=1+gamma

N=20 # number of sites

H=get_NH_SSH(t1,t2,t1,t2,0,N)
h_PBC=H.intra+H.tx+H.tmx  # for OBC use h=H.intra
h_OBC=H.intra

plot_eigen(h_OBC,label="OBC",mode="eigen")
plot_eigen(h_PBC,label="PBC",mode="eigen")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(fontsize=20)
plt.show()

