include("Geometry.jl")
include("hamiltonian.jl")

struct NH_hamiltonian_2D
    intra
    tx
    ty
    tmx
    tmy
    g::geometry
end

function get_2D_square_NH(t1,t2,u;n=2, gamma=0)
  N=n
  g=get_geometry("square",N,N)
  R=g.sites
  n_sites=size(R)[1]
  h_intra=zeros(Complex,n_sites,n_sites)
  h_tx=zeros(Complex,n_sites,n_sites)
  h_ty=zeros(Complex,n_sites,n_sites)
  h_tmx=zeros(Complex,n_sites,n_sites)
  h_tmy=zeros(Complex,n_sites,n_sites)
  for i=1:n_sites
    x=R[i][1]+R[i][2]
    if mod(x,2)==0
      h_intra[i,i]=im*u
    else
      h_intra[i,i]=-im*u
    end
  end
  
  for i=2:n_sites
    for j=1:i-1
      r=R[i]-R[j]
      if 0.99<r[1]<1.01 && abs(r[2])<0.01
        h_intra[i,j]=t1+gamma
        h_intra[j,i]=t1-gamma
      elseif 0.99<-r[1]<1.01 && abs(r[2])<0.01
        h_intra[i,j]=t1-gamma
        h_intra[j,i]=t1+gamma
      elseif 0.99<r[2]<1.01 && abs(r[1])<0.01
        h_intra[i,j]=t2+gamma
        h_intra[j,i]=t2-gamma
      elseif 0.99<-r[2]<1.01 && abs(r[1])<0.01
        h_intra[i,j]=t2-gamma
        h_intra[j,i]=t2+gamma              
      end
    end
  end

  inter_vector=g.inter_vector
  #println(inter_vector)
  A1=inter_vector.x
  A2=inter_vector.y
  
  for i=1:n_sites
    for j=1:n_sites
      r1=R[i]+A1-R[j]
      r2=R[i]+A2-R[j]
      if 0.99<r1[1]<1.01 && abs(r1[2])<0.01
        h_tx[i,j]=t1+gamma
        h_tmx[j,i]=t1-gamma
      elseif 0.99<r2[2]<1.01 && abs(r2[1])<0.01
        h_ty[i,j]=t2+gamma
        h_tmy[j,i]=t2-gamma
      end
    end
  end

  H=NH_hamiltonian_2D(h_intra,h_tx,h_ty,h_tmx,h_tmy,g)

  return H
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
