include("../../../src/NH_hamiltonian.jl")
include("../../../src/NHKPM.jl")
include("../../../src/plot_figs.jl")

t=1
u=0.5
gamma=0.01

N=20 # number of sites

H=get_four_partite_chain(t,u,n=Int(N/4),gamma=gamma)
h_PBC=H.intra+H.tx+H.tmx  # for OBC use h=H.intra
h_OBC=H.intra

plot_eigen(h_OBC,label="OBC",mode="eigen")
plot_eigen(h_PBC,label="PBC",mode="eigen")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(fontsize=20)
plt.show()

