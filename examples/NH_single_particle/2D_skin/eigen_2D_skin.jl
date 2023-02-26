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

plot_eigen(h_OBC,label="OBC",mode="eigen")
plot_eigen(h_PBC,label="PBC",mode="eigen")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(fontsize=20)
plt.show()

