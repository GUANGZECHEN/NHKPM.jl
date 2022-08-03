include("hamiltonian.jl")
include("spin_hamiltonian.jl")
include("General_algorithms.jl")
using LinearAlgebra
using Arpack
using SparseArrays

using BenchmarkTools


struct krylov_decomp
  U                   # krylov space
  B                   # Rayleigh quotient, no matter the type of A (matrix, sp-matrix or MPO), B is always a dense Hessenberg matrix
  u                   # residue vector of keylov space
  b                   # residue, b is a row-vector, when b=0, A=UB is exact
  A                   # original matrix, A=UB+ub
end

function Arnoldi_finite(A,x;m=10,mode="matrix",maxdim=50)              # A: matrix or MPO to diagonalize; x: random initial vector or MPS; m: dim of krylov space
  U=Array{Any}(undef,m+1)
  B=zeros(Complex,(m+1,m))                     # A*U[i]=sum_j B_ji*U[j]; if A is matrix, then A*U=U*B
  U[1]=1/norm(x)*x

  #if mode=="sparse_tensor"
  #  droptol!(U[1],1e-6)
  #end

  for j=1:m
    r=general_product(A,U[j],maxdim=maxdim,mode=mode)
    
    for i=1:j
      B[i,j]=general_inner_product(U[i],r,mode=mode)
      y=-B[i,j]*U[i]                          # multiply with number works with MPS
      r=general_vec_sum(r,y,mode=mode,maxdim=maxdim)        
    end

    B[j+1,j]=norm(r)                           # norm works for both vector and MPS
    U[j+1]=1/B[j+1,j]*r
    #if norm(r)<1e-6                # stop when converged to an invariant space
    #  m=j
    #  break
    #end
  end

  b=B[m+1,1:m]
  b=transpose(b)
  u=U[m+1]

  B=B[1:m,1:m]
  U=U[1:m]

  K=krylov_decomp(U,B,u,b,A)
  return K
end

function verify_krylov_decomp(K::krylov_decomp;mode="matrix",maxdim=50)  # check if a Krylov decomp is correct, i.e. AU=UB+ub; return the average of norm of the error vectors
  U=K.U
  B=K.B
  A=K.A
  u=K.u
  b=K.b
  m=size(B)[2] 
  err=0                                         
  for j=1:m
    r=general_product(A,U[j],maxdim=maxdim,mode=mode)
    scale=norm(r)
    for i=1:m
      y=-B[i,j]*U[i]
      r=general_vec_sum(r,y,mode=mode,maxdim=maxdim)
    end
    y0=-b[j]*u
    r=general_vec_sum(r,y0,mode=mode,maxdim=maxdim)
    err+=norm(r)/scale
  end
  err=err/m
  return err
end

function reorder_schur_vals(vals;which="SR",E=0)                 # find indexes of wanted schur vals
  if which=="SR"
    vals=real(vals)
    perm=sortperm(vals)
  elseif which=="LR"
    vals=real(vals)
    perm=sortperm(vals,rev=true)
  elseif which=="SM"
    vals=abs.(vals)
    perm=sortperm(vals)
  elseif which=="LM"
    vals=abs.(vals)
    perm=sortperm(vals,rev=true)
  elseif which=="specific_E"
    vals_diff=[norm(vals[i]-E) for i=1:size(vals)[1]]
    perm=sortperm(vals_diff)
  end
  return perm
end

function truncate_krylov_decomp(K::krylov_decomp;k=5, which="SR", mode="matrix",maxdim=50, E=0)       # do a unitary trans to make B schur, then do unitary trans to move wanted eigvals top-left, totoal unitary trans is Z, then truncate to k dimension
  B=K.B
  U=K.U
  b=K.b
  u=K.u
  A=K.A

  dim=size(B)[2]
  S=schur(B)
  vals=S.values
  #println(vals)
  perm=reorder_schur_vals(vals,which=which,E=E)
  x=zeros(Bool,dim)
  x[perm[1:k]].=true
  O=ordschur(S,x)                # x is the vector that tells to retain values of indexes i with x[i]=true    
  M,Z=O.T,O.Z                    # the transformed Rayleigh quotient M, and the transformation matrix; B=Z*M*Z'
  b=b*Z                          # transformed residue, now AU=UM+ub; with wanted eigvals already on top-left, we can do the truncation
  
  U_new=Array{Any}(undef,k)      # the transformed and truncated krylov space
  for j=1:k
    U_new[j]=Z[1,j]*U[1]
    for i=2:dim
      y=Z[i,j]*U[i]
      U_new[j]=general_vec_sum(y,U_new[j],mode=mode,maxdim=maxdim)
    end
  end

  M=M[1:k,1:k]
  b=transpose(b[1:k])
  return krylov_decomp(U_new,M,u,b,A)          
end

function expand_krylov_decomp(K::krylov_decomp;m=10,mode="matrix",maxdim=50)        # expand a krylov decomp of dim k to dimension m
  k=size(K.B)[2]
  if !(k<m)
    println("krylov decomp has a larger dim than desired already, expanding Krylov decomp to dim+1")
    m=k+1
  end 
  
  A=K.A
  U=Array{Any}(undef,m+1)
  B=zeros(Complex,(m+1,m))
  U[1:k]=K.U
  U[k+1]=K.u
  B[1:k,1:k]=K.B
  B[k+1,1:k]=K.b

  for j=k+1:m
    r=general_product(A,U[j],maxdim=maxdim,mode=mode)
    for i=1:j
      B[i,j]=general_inner_product(U[i],r,mode=mode)
      y=-B[i,j]*U[i]
      r=general_vec_sum(r,y,mode=mode,maxdim=maxdim)        # multiply with number works with MPS
    end
    B[j+1,j]=norm(r)                           # norm works for both vector and MPS
    U[j+1]=1/B[j+1,j]*r
  end

  b=B[m+1,1:m]
  b=transpose(b)
  u=U[m+1]

  B=B[1:m,1:m]
  U=U[1:m]

  return krylov_decomp(U,B,u,b,A)
end

function get_eigs_krylov(A,x;which="SR",k=5,m=10,mode="matrix",maxerror=1e-3,maxdim=50,maxiter=50)             # get eigs of matrix A, x is random initial vector, k is number of eigs, m is dim of krylov space
  K=Arnoldi_finite(A,x,m=m,mode=mode,maxdim=maxdim)
  K=truncate_krylov_decomp(K,k=k,which=which,mode=mode,maxdim=maxdim)
  println("error in Krylov: ",norm(K.b)/norm(K.B)*sqrt(k))
  num=1
  while norm(K.b)/norm(K.B)*sqrt(k)>maxerror && num<maxiter                                             # if doesn't converge after 50 steps, stop (since it might converge bad) because of MPS truncation error
    K=expand_krylov_decomp(K,m=m,mode=mode,maxdim=maxdim)
    K=truncate_krylov_decomp(K, k=k, which=which,mode=mode,maxdim=maxdim)
    println("error in Krylov: ",norm(K.b)/norm(K.B)*sqrt(k))
    num+=1
  end

  println("steps taken to converge: ",num)
  # add later algorithm to deal with the problem in 3rd and 4th SR: use a smaller k, then re-expand (this helps to get rid of converged but wrong vectors)

  err=verify_krylov_decomp(K,mode=mode,maxdim=maxdim)
  B=K.B
  U=K.U
  F=eigen(B)
  eigvals=F.values
  Z=F.vectors         # Sigma=diagm(eigvals), B=Z*Sigma*inv(Z), AU=UB  -> A*(UZ)=(UZ)*Sigma, columns of UZ are eigvecs
  
  eigvecs=Array{Any}(undef,k)      # the eigvecs
  for j=1:k
    eigvecs[j]=Z[1,j]*U[1]
    for i=2:k
      y=Z[i,j]*U[i]
      eigvecs[j]=general_vec_sum(eigvecs[j],y,mode=mode,maxdim=maxdim)
    end

    if mode=="sparse_tensor"
      droptol!(eigvecs[j],1e-3)
    end
  end

  perm=sortperm(real(eigvals))  
  eigvals=eigvals[perm]
  eigvecs=eigvecs[perm]

  err=0
  for i=1:k
    err+=norm(general_product(A,eigvecs[i],maxdim=maxdim,mode=mode)-eigvals[i]*eigvecs[i])
  end
#  if abs(real(eigvals[1]-eigvals[2]))<1e-5                      # same real part, then choose the positive imaginary part
#    if imag(eigvals[2]-eigvals[1])>1e-5
#      E1,E2=eigvals[1],eigvals[2]
#      psi1,psi2=eigvecs[1],eigvecs[2]
#      eigvals[1],eigvals[2]=E2,E1
#      eigvecs[1],eigvecs[2]=psi2,psi1
#    end
#  end
  #if norm(K.b)/norm(K.B)*sqrt(k)>maxerror
 #   err=1
 # end
  
  return eigvals, eigvecs, err
end

function get_eigs_krylovMPO(H,psi0;which="SR",k=1,m=5,mode="MPO",tol=1e-4,maxerror=1e-4,max_iter=10,maxdim=50,factor=5)  # works best for k=1
  n_iter=0
  err=1
  ini_maxerror=1e-1
  ini_err=ini_maxerror/factor
  psi=psi0

  #E,eigvecs,err=1,psi0,1

  println("starting new iter with: maxerror=",ini_maxerror," tolerance=",ini_err)
  E,eigvecs,err=get_eigs_krylov(H,psi,which=which,k=k,m=m,mode=mode,maxerror=ini_maxerror,maxdim=maxdim)
  println("error in vector: ", err)
  while (err>tol || ini_maxerror>maxerror) && n_iter<max_iter
    ii=0
    while err>ini_err && ii<20   
      #psi=general_vec_sum(eigvecs[1],eigvecs[2],maxdim=maxdim)/2
      psi=eigvecs[1]
      #if mode=="MPO"
      #  psi=truncate(psi,maxdim=Int(maxdim/2),cutoff=1e-10)
      #end
      E,eigvecs,err=get_eigs_krylov(H,psi,which="SR",k=k,m=m,mode="MPO",maxerror=ini_maxerror,maxdim=maxdim)
      println("error in vector: ", err)
      ii+=1
    end
    #psi=general_vec_sum(eigvecs[1],eigvecs[2],maxdim=maxdim)/2
    psi=eigvecs[1]
    #if mode=="MPO"
    #  psi=truncate(psi,maxdim=maxdim/5,cutoff=1e-10)
    #end
    ini_maxerror=min(err,ini_err)  # set the new maxerror (Krylov error) as err (MPO inexact error)
    ini_err=max(ini_maxerror/factor,tol)
    println("starting new iter with: maxerror=",ini_maxerror," tolerance=",ini_err)
    E,eigvecs,err=get_eigs_krylov(H,psi,which=which,k=k,m=m,mode=mode,maxerror=ini_maxerror,maxdim=maxdim)
    println("error in vector: ", err)
    n_iter+=1
  end
  return E,eigvecs,err
end

function get_eigen_krylov_decomp(A,x;which="SR",k=5,m=10,mode="matrix",maxerror=1e-3,maxdim=50,maxiter=50)             # same as get_eigs_krylov, but return the whole decomp instead of only eigvals and eigvecs
  K=Arnoldi_finite(A,x,m=m,mode=mode,maxdim=maxdim)
  K=truncate_krylov_decomp(K,k=k,which=which,mode=mode,maxdim=maxdim)
  println("error in Krylov: ",norm(K.b)/norm(K.B)*sqrt(k))
  num=1
  while norm(K.b)/norm(K.B)*sqrt(k)>maxerror && num<maxiter                                             # if doesn't converge after 50 steps, stop (since it might converge bad) because of MPS truncation error
    K=expand_krylov_decomp(K,m=m,mode=mode,maxdim=maxdim)
    K=truncate_krylov_decomp(K, k=k, which=which,mode=mode,maxdim=maxdim)
    println("error in Krylov: ",norm(K.b)/norm(K.B)*sqrt(k))
    num+=1
  end

  println("steps taken to converge: ",num)
  # add later algorithm to deal with the problem in 3rd and 4th SR: use a smaller k, then re-expand (this helps to get rid of converged but wrong vectors)

  #err=verify_krylov_decomp(K,mode=mode,maxdim=maxdim)
  B=K.B
  U=K.U
  #b=K.b
  #u=K.u
  F=eigen(B)
  eigvals=F.values
  Sigma=diagm(eigvals)
  Z=F.vectors         # Sigma=diagm(eigvals), B=Z*Sigma*inv(Z), AU=UB+ub  -> A*(UZ)=(UZ)*Sigma+ubZ, columns of UZ are eigvecs
  #b=b*Z  

  eigvecs=Array{Any}(undef,k)      # the eigvecs
  for j=1:k
    eigvecs[j]=Z[1,j]*U[1]
    for i=2:k
      y=Z[i,j]*U[i]
      eigvecs[j]=general_vec_sum(eigvecs[j],y,mode=mode,maxdim=maxdim)
    end
  end

  #Sigma=Sigma[1:k,1:k]
  #b=transpose(b[1:k])
   

  err=0              # err will be the new u, which is the max MPO error vector
  err_vec=Array{Any}(undef,k)
  index_maxerr=1
  for i=1:k
    err_vec[i]=general_product(A,eigvecs[i],maxdim=maxdim,mode=mode)-eigvals[i]*eigvecs[i]  # AU_i=U_jB_ji+err_i, project all things to max(err_i)
    if norm(err_vec[i])>err
      err=norm(err_vec[i])
      index_maxerr=i
    end    
  end

  u=err_vec[index_maxerr]/norm(err_vec[index_maxerr])
  b=Array{Any}(undef,k)
  for i=1:k
    b[i]=dot(err_vec[i],u)
  end
  #println(b)
  K_new=krylov_decomp(eigvecs,Sigma,u,b,A) 
  
  return K_new, err
end

function careful_converge_krylov_MPO(A,x;which="SR",k=5,m=10,mode="matrix",tol=1e-2,maxerror=1e-5,maxdim=50,maxiter=50,dim_increase=30)    # this method has limitation: no matter how large maxdim is, it will have an error in vector (>1e-5 for L=8, which means the vector is bad)
  maxdim_0=10
  K,err=get_eigen_krylov_decomp(A,x,which=which,k=k,m=m,mode=mode,maxerror=maxerror,maxdim=maxdim_0,maxiter=maxiter)
  println("maxdim: ",maxdim_0," err: ",err)
  n_iter=0 
  while err>tol && maxdim_0<maxdim
    n_iter+=1
    if n_iter==10
      maxdim_0+=dim_increase
      n_iter=0
    end
    err_Krylov=1
    num=0
    while err_Krylov>maxerror && num<maxiter
      K=expand_krylov_decomp(K,m=m,mode=mode,maxdim=maxdim_0)
      K=truncate_krylov_decomp(K, k=k, which=which,mode=mode,maxdim=maxdim_0)
      err_Krylov=norm(K.b)/norm(K.B)*sqrt(k)
      num+=1
      println("error in Krylov: ",err_Krylov)
    end

    B=K.B
    U=K.U
    F=eigen(B)
    eigvals=F.values
    Sigma=diagm(eigvals)
    Z=F.vectors         # Sigma=diagm(eigvals), B=Z*Sigma*inv(Z), AU=UB+ub  -> A*(UZ)=(UZ)*Sigma+ubZ, columns of UZ are eigvecs  

    eigvecs=Array{Any}(undef,k)      # the eigvecs
    for j=1:k
      eigvecs[j]=Z[1,j]*U[1]
      for i=2:k
        y=Z[i,j]*U[i]
        eigvecs[j]=general_vec_sum(eigvecs[j],y,mode=mode,maxdim=maxdim_0)
      end
    end

    err=0              # err will be the new u, which is the max MPO error vector
    err_vec=Array{Any}(undef,k)
    index_maxerr=1
    for i=1:k
      err_vec[i]=general_product(A,eigvecs[i],maxdim=maxdim,mode=mode)-eigvals[i]*eigvecs[i]  # AU_i=U_jB_ji+err_i, project all things to max(err_i)
      if norm(err_vec[i])>err
        err=norm(err_vec[i])
        index_maxerr=i
      end    
    end
    println("err vec: ",[norm(err_vec[i]) for i=1:k])

    u=err_vec[index_maxerr]/norm(err_vec[index_maxerr])
    b=Array{Any}(undef,k)
    for i=1:k
      b[i]=dot(err_vec[i],u)
    end
   # println("b vector: ",b)

    K=krylov_decomp(eigvecs,Sigma,u,b,A)  
    println("maxdim: ",maxdim_0," err: ",err)  
  end
  return K
end

function refine_Krylov_eigvec(A,x,E;k=1,m=5,mode="MPO",maxerror=1e-5,maxdim=50,maxiter=100)    # Ax=Ex
  K=Arnoldi_finite(A,x,m=m,mode=mode,maxdim=maxdim)
  K=truncate_krylov_decomp(K,k=k,which="specific_E",mode=mode,maxdim=maxdim,E=E)
  U=K.U
  B=K.B
  err_vec=general_product(A,U[1],maxdim=maxdim,mode=mode)-B[1,1]*U[1]
  println("error in refined Krylov: ",norm(err_vec))
  num=1
  while norm(err_vec)>maxerror && num<maxiter                                             
    K=expand_krylov_decomp(K,m=m,mode=mode,maxdim=maxdim)
    K=truncate_krylov_decomp(K,k=k,which="specific_E",mode=mode,maxdim=maxdim,E=E)
    U=K.U
    B=K.B
    err_vec=general_product(A,U[1],maxdim=maxdim,mode=mode)-B[1,1]*U[1]
    println("error in refined Krylov: ",norm(err_vec))
    num+=1
  end
  return B[1,1], U[1]
end
