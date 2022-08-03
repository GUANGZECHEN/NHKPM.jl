# test performance of Krylov algorithm in Hermitian case, compared with DMRG

include("../../src/KrylovSchur.jl")

function test_krylov_tensor(N,J,B;boundary="OBC")
  H=get_generalized_SSH_spin_chain(J,B,N=N,mode="sparse_tensor",boundary=boundary)

  dim=2^N
  x=sprand(dim,0.5)+im*sprand(dim,0.5)
  E,eigvecs,err=get_eigs_krylov(H,x,which="SR",k=1,m=3,maxerror=1e-4,mode="sparse_tensor")
  println("krylov tensor energy: ",E[1])
  println("error vector norm: ",err)
  #println("krylov tensor vec: ",eigvecs[1])
end

function test_krylov_MPO(N,J,B;dim=1,m=5,tol=1e-4,maxerror=1e-4,maxdim=50,boundary="OBC")     # tol is the tolerance of the norm of the error vector due to inexactness of MPS + and *, for larger systems tol cannot be too small (otherwise convergence is slow) # tol and maxerror are relative
  H=get_generalized_SSH_spin_chain(J,B,N=N,mode="MPO",boundary=boundary)
  sites=firstsiteinds(H)

  psi = randomMPS(sites,dim)   # dim better to be 1 due to inexactness of + and * for MPS

  E,eigvecs,err=get_eigs_krylovMPO(H,psi,which="SR",k=1,m=m,mode="MPO",tol=tol,maxerror=maxerror,max_iter=10,maxdim=maxdim)

  psi=eigvecs[1]
  E2=inner(psi,H*psi)
  Id=get_Id_MPO(sites)
  err=norm(Id*(H*psi)-E2*psi)/norm(psi)
  println("krylov MPO energy: ",E[1])
  println("true eigval: ",E2)
  println("error in state: ", err)
end

function test_DMRG(N,J,B)   
  H2=get_generalized_SSH_spin_chain(J,B,N=N,mode="MPO")
  sites=firstsiteinds(H2)
  sweeps = Sweeps(10) # number of sweeps is 10
  maxdim!(sweeps,10,10,10,10,10) # gradually increase states kept
  cutoff!(sweeps,1E-10) # desired truncation error

  psi0 = randomMPS(sites,10)

  energy,psi = dmrg(H2,psi0,sweeps)
  println("DMRG result: ",energy)
end

function test_performance()
  N=8
  J=[1,1,1,1,1,1]
  B=[0,0,0,0,0,0]

  @time test_DMRG(N,J,B)
  #@time test_krylov_tensor(N,J,B,boundary="OBC")
  @time test_krylov_MPO(N,J,B,dim=1,m=5,boundary="OBC")   # PBC not friendly to MPO
end

test_performance()
