include("Geometry.jl")
using DelimitedFiles
using Random

struct hamiltonian
    intra
    tx
    ty
    txy
    txmy
    g::geometry
end

struct hamiltonian_1D
    intra
    tx
end


function neighbor(site_a,site_b,lattice)
  d=norm(site_a-site_b)
  if lattice=="triangular" || lattice=="honeycomb"
    if abs(d-1)<0.01
      return 1  # nearest
    elseif abs(d-sqrt(3))<0.01
      return 2
    elseif abs(d-2)<0.01
      return 3
    else
      return 0  # larger than 3rd neighbor
    end
  elseif lattice=="square"
    if abs(d-1)<0.01
      return 1  # nearest
    elseif abs(d-sqrt(2))<0.01
      return 2
    elseif abs(d-2)<0.01
      return 3
    else
      return 0
    end
  else
    println("wrong lattice") 
  end
end

function get_H_component(R,r,t,lattice)
  N=size(R,1)
  H=zeros(Complex,(N,N))
  for i=1:N
    for j=1:N
      R1=R[i]
      R2=R[j]+r
      index=neighbor(R1,R2,lattice)
      if index>0
        H[i,j]=t[index]                  # H[i,j] is the hopping from j to i
      end
    end
  end
  return H
end

function get_H(g::geometry,t1,t2=0,t3=0)
  t=[t1,t2,t3]
  R=g.sites
  inter_vec=g.inter_vector
  lattice=g.lattice
  h_intra=get_H_component(R,[0,0],t,lattice)
  hx=get_H_component(R,inter_vec.x,t,lattice)
  hy=get_H_component(R,inter_vec.y,t,lattice)
  hxy=get_H_component(R,inter_vec.xy,t,lattice)
  hxmy=get_H_component(R,inter_vec.xmy,t,lattice)

  H=hamiltonian(h_intra,hx,hy,hxy,hxmy,g)
end

function get_Hk(H::hamiltonian,k::Array{Float64,1})
  inter_vec=H.g.inter_vector
  H0=exp(im*dot(k,inter_vec.x))*H.tx+exp(im*dot(k,inter_vec.y))*H.ty+exp(im*dot(k,inter_vec.xy))*H.txy+exp(im*dot(k,inter_vec.xmy))*H.txmy
  Hk=H.intra+H0+adjoint(H0)  
  return Hk 
end

function plot_H(H,R)
  plt.figure(figsize=(6,6),dpi=80)
  N=size(R,1)
  for i=1:N
    for j=1:N
      if real(H[i,j])>0.1
        plt.plot((R[i][1],R[j][1]), (R[i][2],R[j][2]),color="red",linewidth=norm(H[i,j])*1)
      elseif real(H[i,j])<-0.1
        plt.plot((R[i][1],R[j][1]), (R[i][2],R[j][2]),color="blue",linewidth=norm(H[i,j])*1)
      end
    end
  end
  plt.axis([-5,5,-5,5])
  plt.xlabel("")
  plt.ylabel("")
  plt.show()
end

function test_hamiltonian()
  n=2
  m=2
  lattice="triangular"

  g=get_geometry(lattice,n,m)
  t=1
  H=get_H(g,t)
  Hk=get_Hk(H,[pi/2,0.0])
  plot_H(Hk,g.sites)
end

#test_hamiltonian()
