include("../../../src/NH_hamiltonian.jl")
include("../../../src/NHKPM.jl")
include("../../../src/plot_figs.jl")

symmetry="YSY"

t1=1/6
t2=1
gamma=4/3   

v1=t1-gamma/2
v2=t1+gamma/2
w1=t2
w2=t2

N=20 # number of sites

H=get_NH_SSH(v1,v2,w1,w2,0,N)
h_PBC=H.intra+H.tx+H.tmx  # for OBC use h=H.intra
h_OBC=H.intra

plot_eigen(h_OBC,label="OBC",mode="eigen")
plot_eigen(h_PBC,label="PBC",mode="eigen")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(fontsize=20)
plt.show()

