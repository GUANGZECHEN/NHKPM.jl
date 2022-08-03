include("Geometry.jl")
include("hamiltonian.jl")
using DelimitedFiles
using Random

struct NH_hamiltonian_1D
    intra
    tx
    tmx
    g::geometry_1D
end

function get_four_partite_chain(t,u;n=2,gamma=0) # n is number of supercells
  N=4*n
  g = get_geometry_1D(N)
  h_intra = zeros(Complex,(N,N))
  for i=1:n-1
    h_intra[4*i-3,4*i-3]=+im*u*0
    h_intra[4*i-2,4*i-2]=-im*u
    h_intra[4*i-1,4*i-1]=-im*u
    h_intra[4*i,4*i]=+im*u*0

    h_intra[4*i-3,4*i-2]=t+gamma
    h_intra[4*i-2,4*i-3]=t-gamma
    h_intra[4*i-2,4*i-1]=t+gamma
    h_intra[4*i-1,4*i-2]=t-gamma
    h_intra[4*i-1,4*i]=t+gamma
    h_intra[4*i,4*i-1]=t-gamma
    h_intra[4*i,4*i+1]=t+gamma
    h_intra[4*i+1,4*i]=t-gamma
  end
  h_intra[4*n-3,4*n-3]=+im*u*0
  h_intra[4*n-2,4*n-2]=-im*u
  h_intra[4*n-1,4*n-1]=-im*u
  h_intra[4*n,4*n]=+im*u*0
  
  h_intra[4*n-3,4*n-2]=t+gamma
  h_intra[4*n-2,4*n-3]=t-gamma
  h_intra[4*n-2,4*n-1]=t+gamma
  h_intra[4*n-1,4*n-2]=t-gamma
  h_intra[4*n-1,4*n]=t+gamma
  h_intra[4*n,4*n-1]=t-gamma

  h_tx=zeros(Complex,(N,N))
  h_tx[N,1]=t

  h_tmx=zeros(Complex,(N,N))
  h_tmx[1,N]=t

  H=NH_hamiltonian_1D(h_intra,h_tx,h_tmx,g)
  return H
end

function get_NH_SSH(v1,v2,w1,w2,u,N)   # generic 1D bipartite NH model, from PHYSICAL REVIEW B 97, 045106 (2018), N is chain length
  g = get_geometry_1D(N)
  h_intra = zeros(Complex,(N,N))
  n = Int(N/2)
  for i=1:n-1
    h_intra[2*i-1,2*i-1]=+im*u
    h_intra[2*i,2*i]=-im*u
    h_intra[2*i-1,2*i]=v2               # hopping from 2*i to 2*i-1, i.e. B to A
    h_intra[2*i,2*i-1]=v1
    h_intra[2*i,2*i+1]=w2
    h_intra[2*i+1,2*i]=w1
  end
  h_intra[2*n-1,2*n-1]=+im*u
  h_intra[2*n,2*n]=-im*u
  h_intra[2*n-1,2*n]=v2               
  h_intra[2*n,2*n-1]=v1  
  
  h_tx = zeros(Complex,(N,N))
  h_tx[2*n,1]=w2     
  
  h_tmx = zeros(Complex,(N,N))
  h_tmx[1,2*n]=w1  

  H=NH_hamiltonian_1D(h_intra,h_tx,h_tmx,g)
  return H
end

function add_t2_1D(H,t2)     # add NNN hopping to 1D Hamiltonian
  h_intra=H.intra
  N=size(h_intra)[1]
  for i=1:N-2
    h_intra[i,i+2]=t2
    h_intra[i+2,i]=t2
  end
  h_tx=H.tx
  h_tmx=H.tmx
  h_tx[N-1,1]=t2
  h_tx[N,2]=t2
  h_tmx[1,N-1]=t2
  h_tmx[2,N]=t2

  H2=NH_hamiltonian_1D(h_intra,h_tx,h_tmx,H.g)
  return H2
end

function get_Hk_NH_1D(H::NH_hamiltonian_1D,k)
  inter_vec=H.g.x
  Hk=H.intra+exp(im*k*inter_vec)*H.tx+exp(-im*k*inter_vec)*H.tmx
  return Hk 
end


function Hermitrize(h)
  dim=size(h)[1]
  Dim=2*dim
  H=zeros(Complex,(Dim,Dim))
  for i=1:dim
    for j=1:dim
      H[i,j+dim]=h[i,j]
      H[j+dim,i]=conj(h[i,j])
    end
  end
  return H
end


function test_NH_hamiltonian()
  H=get_NH_SSH(1,2,1,2,0,4)
  print(H.intra,"\n")
  H2=Hermitrize(H.intra)
  print(H2)
end

#test_NH_hamiltonian()
