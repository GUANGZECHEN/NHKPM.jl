include("../../src/special_spin_hamiltonian.jl")
include("../../src/NHKPM.jl")
include("../../src/KrylovSchur.jl")

using LinearAlgebra
using Arpack

using PyCall
plt = pyimport("matplotlib.pyplot")

t=1
N=20
U=5

H=zeros(N,N)
for i=1:(N-1)
  H[i,i+1]=t
  H[i+1,i]=t
  H[i,i]=(i-1)*U
end

H[N,N]=(N-1)*U

F=eigen(Matrix(H))
Es=F.values
psis=F.vectors

exp_ns=zeros(N,N)
for i=1:N
  for j=1:N
    exp_ns[i,j]=abs(psis[i,j])^2
  end
end


println(Matrix(H))
println(Es)

for i=1:N
  plt.plot(range(1,N,N),exp_ns[:,i])
end
plt.show()

