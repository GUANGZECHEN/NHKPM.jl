using ITensors
using LinearAlgebra

# basic operations

function get_zero_MPS(sites;maxdim=50)
  N=size(sites)[1]
  psi=randomMPS(sites,maxdim)
  for i=1:N
    psi[i].=0
  end
  return psi
end

function get_H_dag_MPO(H)
  H2=deepcopy(H)
  sites = firstsiteinds(H2)
  N=size(sites)[1]

  i1=inds(H2[1])[2]
  i2=inds(H2[1])[3]
  ITensors.swapind!(H2[1],i1,i2)
  H2[1]=conj(H2[1])
  for i=2:N-1
    i1=inds(H2[i])[3]
    i2=inds(H2[i])[4]
    ITensors.swapind!(H2[i],i1,i2)
    H2[i]=conj(H2[i])
  end
  i1=inds(H2[N])[2]
  i2=inds(H2[N])[3]
  ITensors.swapind!(H2[N],i1,i2)
  H2[N]=conj(H2[N])
  return H2
end

function get_Id_MPO(sites)      # sites should be siteinds object
  h=OpSum()
  h+="Id",1
  H=MPO(h,sites)
  return H
end

function permute_inds_MPO(H)  # without changing H
  H2=deepcopy(H)
  sites = firstsiteinds(H2)
  N=size(sites)[1]

  i1=inds(H2[1])[1]
  i2=inds(H2[1])[2]
  i3=inds(H2[1])[3]
  H2[1]=ITensors.permute(H2[1],i3,i2,i1)
  for i=2:N-1
    i1=inds(H2[i])[1]
    i2=inds(H2[i])[2]
    i3=inds(H2[i])[3]
    i4=inds(H2[i])[4]
    H2[i]=ITensors.permute(H2[i],i4,i3,i2,i1)
  end
  i1=inds(H2[N])[1]
  i2=inds(H2[N])[2]
  i3=inds(H2[N])[3]
  H2[N]=ITensors.permute(H2[N],i3,i2,i1)
  return H2
end

function permute_inds_MPS(psi)  # without changing psi
  psi2=deepcopy(psi)
  N=size(psi2)[1]

  i1=inds(psi2[1])[1]
  i2=inds(psi2[1])[2]
  psi2[1]=ITensors.permute(psi2[1],i2,i1)
  for i=2:N-1
    i1=inds(psi2[i])[1]
    i2=inds(psi2[i])[2]
    i3=inds(psi2[i])[3]
    psi2[i]=ITensors.permute(psi2[i],i3,i1,i2)
  end
  i1=inds(psi2[N])[1]
  i2=inds(psi2[N])[2]
  psi2[N]=ITensors.permute(psi2[N],i2,i1)
  return psi2
end


#algebraic operations


function general_product(A,x;mode="matrix",maxdim=100)
  if mode=="matrix" || mode=="sparse_tensor"
    if norm(x)<1e-10
      return x
    end
    y=A*x
    return y
  else
    sites=firstsiteinds(A)
    y=contract(A,x,cutoff=1e-12,maxdim=maxdim)
    Identity=get_Id_MPO(sites)
    y=contract(Identity,y,maxdim=maxdim)
    return y
  end
end

function general_inner_product(x,y;mode="matrix")   # inner product for MPS is same as for matrix
  z=dot(x,y)
  return z
end

function general_vec_sum(x,y;mode="matrix",maxdim=100)   # the first vector cannot be zero vector in the sum function
  if mode=="matrix" || mode=="sparse_tensor"
    z=x
    if typeof(y)!=SparseVector{Any, Int64}
      z=x+y
    end
    return z
  else
    z=x
    if norm(y)>0
      z=+(z,y,cutoff=1e-12,maxdim=maxdim)
    end
    return z
  end
end

function general_vec_sum_2(xs;mode="matrix",maxdim=100)  # summing over more than one vector
  N=size(xs)[1]
  if mode=="matrix" || mode=="sparse_tensor"
    z=xs[1]
    for i=2:N   
      if typeof(xs[i])!=SparseVector{Any, Int64}
        z=z+xs[i]
      end
    end
    return z
  else
    z=xs[1]
    for i=2:N
      if norm(xs[i])>1e-10
        z=+(z,xs[i],cutoff=1e-12,maxdim=maxdim)
      end
    end
    return z
  end
end


