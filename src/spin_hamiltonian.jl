include("Geometry.jl")
include("General_algorithms.jl")
using DelimitedFiles
using Random

using ITensors
using SparseArrays

function get_vec_S(dof)  # get spin operators for spin-S, dof=2S+1
  if dof==2
    sp=[0 2;0 0]/2        # S+, divide by 2 because of spin-1/2
    sm=[0 0;2 0]/2        # S-
    sz=[1 0;0 -1]/2
    id=[1 0;0 1]
  elseif dof==3
    sp=[0 2 0;0 0 2;0 0 0]/sqrt(2)       # S+
    sm=[0 0 0;2 0 0;0 2 0]/sqrt(2)     # S-
    sz=[1 0 0;0 0 0;0 0 -1]
    id=[1 0 0;0 1 0;0 0 1]
  end    
  return sp,sm,sz,id
end

function get_op_S(i,dof_sites)  # given dof on each site, create the S_i operators in tensor form
  N=size(dof_sites)[1]
  
  sp,sm,sz,id=get_vec_S(dof_sites[1])
  if i==1
    Sp=sp
    Sm=sm
    Sz=sz
  else
    Sp=id
    Sm=id
    Sz=id   
  end  
  
  for k=2:N
    sp,sm,sz,id=get_vec_S(dof_sites[k])
    if k==i
      Sp=kron(Sp,sp)
      Sm=kron(Sm,sm)
      Sz=kron(Sz,sz)
    else
      Sp=kron(Sp,id)
      Sm=kron(Sm,id)
      Sz=kron(Sz,id)
    end
  end   
  return Sp,Sm,Sz      
end 

function get_truncated_Haldane_chain(J,hz,N)
  dim=3^(N-2)*2^2
  H=spzeros(dim,dim)
  dof_sites=ones(N)*3
  dof_sites[1]=2
  dof_sites[N]=2  
  
  for i=1:N-1
    Sp,Sm,Sz=get_op_S(i,dof_sites)
    Sp_1,Sm_1,Sz_1=get_op_S(i+1,dof_sites)
    X=(Sp*Sm_1/2+Sm*Sp_1/2+Sz*Sz_1)
    H+=J*X+J/3*X*X
    if i==2 || i==3 || i==6 || i==7
      H+=-im*hz*Sz
    end    
  end
  return H
end

function Si(i,N;S=1/2,mode="sparse_tensor")    # tensor form of Si, for spin chain of length N
  if S==1/2
    sp=[0 2;0 0]/2        # S+, divide by 2 because of spin-1/2
    sm=[0 0;2 0]/2        # S-
    sz=[1 0;0 -1]/2
    id=[1 0;0 1]
  elseif S==1
    sp=[0 2 0;0 0 2;0 0 0]/sqrt(2)       # S+
    sm=[0 0 0;2 0 0;0 2 0]/sqrt(2)     # S-
    sz=[1 0 0;0 0 0;0 0 -1]
    id=[1 0 0;0 1 0;0 0 1]
  end    

  if mode=="sparse_tensor"
    sp=sparse(sp)
    sm=sparse(sm)
    sz=sparse(sz)
    id=sparse(id)
  end

  if i==1
    Sp=sp
    Sm=sm
    Sz=sz
  else
    Sp=id
    Sm=id
    Sz=id   
  end

  for k=2:N
    if k==i
      Sp=kron(Sp,sp)
      Sm=kron(Sm,sm)
      Sz=kron(Sz,sz)
    else
      Sp=kron(Sp,id)
      Sm=kron(Sm,id)
      Sz=kron(Sz,id)
    end
  end   
  return Sp,Sm,Sz  
end

function get_Haldane_chain(J,hz,N)
  dim=3^N
  H=spzeros(dim,dim)
  for i=1:N-1
    Sp,Sm,Sz=Si(i,N,S=1)
    Sp_1,Sm_1,Sz_1=Si(i+1,N,S=1)
    X=(Sp*Sm_1/2+Sm*Sp_1/2+Sz*Sz_1)
    H+=J*X+J/3*X*X
    if i==2 || i==3 || i==6 || i==7
      H+=-im*hz*Sz
    end  
  end
  return H
end

function Si_Sj(i,j,N;spin=1/2,mode="sparse_tensor")    # tensor form of Si*Sj, for spin chain of length N; Mp=S+i*S-j, Mm=S-i*S+j, Mz=Szi*Szj
  sp=[0 2;0 0]/2        # S+, divide by 2 because of spin-1/2
  sm=[0 0;2 0]/2        # S-
  sz=[1 0;0 -1]/2
  id=[1 0;0 1]

  if mode=="sparse_tensor"
    sp=sparse(sp)
    sm=sparse(sm)
    sz=sparse(sz)
    id=sparse(id)
  end

  if i==1
    Mp=sp
    Mm=sm
    Mz=sz
  else
    Mp=id
    Mm=id
    Mz=id   
  end

  for k=2:N
    if k==i
      Mp=kron(Mp,sp)
      Mm=kron(Mm,sm)
      Mz=kron(Mz,sz)
    elseif k==j
      Mp=kron(Mp,sm)
      Mm=kron(Mm,sp)
      Mz=kron(Mz,sz)
    else
      Mp=kron(Mp,id)
      Mm=kron(Mm,id)
      Mz=kron(Mz,id)
    end
  end
  return Mp,Mm,Mz
end

function get_generalized_SSH_spin_chain(J,B;N=10,mode="MPO",boundary="OBC",sites=0,boundary_factor=1)      # generalized dimerized spin chain, J and B are vectors; N is chain length; mode switch between MPO and Tensor representation.
  J2=zeros(Complex,6)
  B2=zeros(Complex,6)
  J2[1:size(J)[1]].=J[:]
  B2[1:size(B)[1]].=B[:]

  J=J2
  B=B2
  if mode=="MPO"
    if sites==0
      sites = siteinds("S=1/2",N)
    end
    h = OpSum()
    for i=1:Int(N/2)-1
      h += 0.5*J[2],"S+",2*i-1,"S-",2*i                                                  # be careful I devide 1/2 here, which means the parameter I put is J_perp
      h += 0.5*J[3],"S-",2*i-1,"S+",2*i
      h += J[1],"Sz",2*i-1,"Sz",2*i

      h += B[1],"Sz",2*i-1
      h += B[2],"S+",2*i-1
      h += B[3],"S-",2*i-1

      h += 0.5*J[5],"S+",2*i,"S-",2*i+1
      h += 0.5*J[6],"S-",2*i,"S+",2*i+1
      h += J[4],"Sz",2*i,"Sz",2*i+1

      h += B[4],"Sz",2*i
      h += B[5],"S+",2*i
      h += B[6],"S-",2*i
    end
    h += 0.5*J[2],"S+",N-1,"S-",N
    h += 0.5*J[3],"S-",N-1,"S+",N
    h += J[1],"Sz",N-1,"Sz",N

    h += B[1],"Sz",N-1
    h += B[2],"S+",N-1
    h += B[3],"S-",N-1
    h += B[4],"Sz",N
    h += B[5],"S+",N
    h += B[6],"S-",N

    if boundary=="PBC"
      h += boundary_factor*0.5*J[6],"S+",1,"S-",N
      h += boundary_factor*0.5*J[5],"S-",1,"S+",N
      h += boundary_factor*J[4],"Sz",1,"Sz",N
    end

    H=MPO(h,sites)

    return H
  else
    dim=2^N
    H=spzeros(dim,dim)
    for i=1:Int(N/2)-1
      Mp,Mm,Mz=Si_Sj(2*i-1,2*i,N,mode=mode)
      Sp,Sm,Sz=Si(2*i-1,N,mode=mode)
      H+=0.5*J[2]*Mp+0.5*J[3]*Mm+J[1]*Mz
      H+=B[1]*Sz+B[2]*Sp+B[3]*Sm

      Mp,Mm,Mz=Si_Sj(2*i,2*i+1,N)
      Sp,Sm,Sz=Si(2*i,N)
      H+=0.5*J[5]*Mp+0.5*J[6]*Mm+J[4]*Mz
      H+=B[4]*Sz+B[5]*Sp+B[6]*Sm
    end
    Mp,Mm,Mz=Si_Sj(N-1,N,N)
    Sp,Sm,Sz=Si(N-1,N)
    H+=0.5*J[2]*Mp+0.5*J[3]*Mm+J[1]*Mz
    H+=B[1]*Sz+B[2]*Sp+B[3]*Sm

    Sp,Sm,Sz=Si(N,N)
    H+=B[4]*Sz+B[5]*Sp+B[6]*Sm

    if boundary=="PBC"
      Mp,Mm,Mz=Si_Sj(1,N,N)
      H+=boundary_factor*(0.5*J[5]*Mm+0.5*J[6]*Mp+J[4]*Mz)                     # here mp and mm are reversed
    end

    return H
  end  
end

function add_J2_1D(H,J2;mode="MPO",boundary="OBC")
  if mode=="MPO"
    sites=firstsiteinds(H)
    N=size(sites)[1]
    h=OpSum()
    for i=1:N-2
      h += J2[2]/2,"S+",i,"S-",i+2
      h += J2[3]/2,"S-",i,"S+",i+2
      h += J2[1],"Sz",i,"Sz",i+2
    end

    if boundary=="PBC"
      h += 0.5*J2[2],"S+",2,"S-",N
      h += 0.5*J2[3],"S-",2,"S+",N
      h += J2[1],"Sz",2,"Sz",N

      h += 0.5*J2[2],"S+",1,"S-",N-1
      h += 0.5*J2[3],"S-",1,"S+",N-1
      h += J2[1],"Sz",1,"Sz",N-1
    end

    H2=MPO(h,sites)
    H=H+H2
    return H
  else
    dim=size(H)[1]
    N=Int(log2(dim))
    for i=1:N-2
      Sp_0,Sm_0,Sz_0=Si(i,N)
      Sp_2,Sm_2,Sz_2=Si(i+2,N)
      H += J2[2]/2*Sp_0*Sm_2
      H += J2[3]/2*Sm_0*Sp_2
      H += J2[1]*Sz_0*Sz_2
    end

    if boundary=="PBC"
      Sp_0,Sm_0,Sz_0=Si(N-1,N)
      Sp_1,Sm_1,Sz_1=Si(N,N)
      Sp_2,Sm_2,Sz_2=Si(1,N)
      Sp_3,Sm_3,Sz_3=Si(2,N)

      H += J2[2]/2*Sp_0*Sm_2
      H += J2[3]/2*Sm_0*Sp_2
      H += J2[1]*Sz_0*Sz_2

      H += J2[2]/2*Sp_1*Sm_3
      H += J2[3]/2*Sm_1*Sp_3
      H += J2[1]*Sz_1*Sz_3
    end

    return H    
  end  
end

function add_t2_JW_1D(H,t2;mode="MPO")
  if mode=="MPO"
    sites=firstsiteinds(H)
    N=size(sites)[1]
    h=OpSum()
    for i=1:N-2
      h += -2*t2,"S+",i,"Sz",i+1,"S-",i+2
      h += -2*t2,"S-",i,"Sz",i+1,"S+",i+2
    end
    H2=MPO(h,sites)
    H=H+H2
    return H
  else
    dim=size(H)[1]
    N=Int(log2(dim))
    for i=1:N-2
      Sp_0,Sm_0,Sz_0=Si(i,N)
      Sp_1,Sm_1,Sz_1=Si(i+1,N)
      Sp_2,Sm_2,Sz_2=Si(i+2,N)
      H += -2*t2*Sp_0*Sz_1*Sm_2
      H += -2*t2*Sm_0*Sz_1*Sp_2
    end
    return H    
  end
  return 0
end

function get_four_partite_imag_onsite(H,mu;mode="MPO")
  if mode=="MPO"
    sites=firstsiteinds(H)
    N=size(sites)[1]
    h=OpSum()
    for i=1:Int(N/4)
      h += -im*mu,"Sz",4*i-2
      h += -im*mu,"Sz",4*i-1
    end
    H2=MPO(h,sites)
    H=H+H2
    return H
  else
    dim=size(H)[1]
    N=Int(log2(dim))
    for i=1:Int(N/4)
      Sp,Sm,Sz=Si(4*i-2,N)
      H+=-im*mu*Sz
      Sp,Sm,Sz=Si(4*i-1,N)
      H+=-im*mu*Sz
    end
    return H    
  end  
end

function get_four_partite_spin_chain(J,mu;B=0,N=2,mode="MPO",BC="OBC")  # N is number of supercells
  if mode=="MPO"
    sites = siteinds("S=1/2",N*4)
    h = OpSum()
    for i=1:N-1
      h += 0.5*J[2],"S+",4*i-3,"S-",4*i-2
      h += 0.5*J[3],"S-",4*i-3,"S+",4*i-2
      h += J[1],"Sz",4*i-3,"Sz",4*i-2

      h += 0.5*J[2],"S+",4*i-2,"S-",4*i-1
      h += 0.5*J[3],"S-",4*i-2,"S+",4*i-1
      h += J[1],"Sz",4*i-2,"Sz",4*i-1

      h += 0.5*J[2],"S+",4*i-1,"S-",4*i
      h += 0.5*J[3],"S-",4*i-1,"S+",4*i
      h += J[1],"Sz",4*i-1,"Sz",4*i

      h += 0.5*J[2],"S+",4*i,"S-",4*i+1
      h += 0.5*J[3],"S-",4*i,"S+",4*i+1
      h += J[1],"Sz",4*i,"Sz",4*i+1

      #h += im*mu,"Sz",4*i-3
      h += -im*mu,"Sz",4*i-2
      h += -im*mu,"Sz",4*i-1
      #h += im*mu,"Sz",4*i
    end
    h += 0.5*J[2],"S+",4*N-3,"S-",4*N-2
    h += 0.5*J[3],"S-",4*N-3,"S+",4*N-2
    h += J[1],"Sz",4*N-3,"Sz",4*N-2

    h += 0.5*J[2],"S+",4*N-2,"S-",4*N-1
    h += 0.5*J[3],"S-",4*N-2,"S+",4*N-1
    h += J[1],"Sz",4*N-2,"Sz",4*N-1

    h += 0.5*J[2],"S+",4*N-1,"S-",4*N
    h += 0.5*J[3],"S-",4*N-1,"S+",4*N
    h += J[1],"Sz",4*N-1,"Sz",4*N

    #h += im*mu,"Sz",4*N-3
    h += -im*mu,"Sz",4*N-2
    h += -im*mu,"Sz",4*N-1
    #h += im*mu,"Sz",4*N

    for i=1:4*N
      h+=B,"Sz",i
    end

    H=MPO(h,sites)
    return H
  else
    N=N*4
    dim=2^N
    H=spzeros(dim,dim)
    for i=1:Int(N/4)-1
      Mp,Mm,Mz=Si_Sj(4*i-3,4*i-2,N)
      H+=0.5*(J[2]*Mp+J[3]*Mm)+J[1]*Mz

      Mp,Mm,Mz=Si_Sj(4*i-2,4*i-1,N)
      H+=0.5*(J[2]*Mp+J[3]*Mm)+J[1]*Mz

      Mp,Mm,Mz=Si_Sj(4*i-1,4*i,N)
      H+=0.5*(J[2]*Mp+J[3]*Mm)+J[1]*Mz

      Mp,Mm,Mz=Si_Sj(4*i,4*i+1,N)
      H+=0.5*(J[2]*Mp+J[3]*Mm)+J[1]*Mz

      Sp,Sm,Sz=Si(4*i-3,N)
     # H+=im*mu*Sz
      Sp,Sm,Sz=Si(4*i-2,N)
      H+=-im*mu*Sz
      Sp,Sm,Sz=Si(4*i-1,N)
      H+=-im*mu*Sz
      Sp,Sm,Sz=Si(4*i,N)
     # H+=im*mu*Sz
    end
    Mp,Mm,Mz=Si_Sj(N-3,N-2,N)
    H+=0.5*(J[2]*Mp+J[3]*Mm)+J[1]*Mz

    Mp,Mm,Mz=Si_Sj(N-2,N-1,N)
    H+=0.5*(J[2]*Mp+J[3]*Mm)+J[1]*Mz

    Mp,Mm,Mz=Si_Sj(N-1,N,N)
    H+=0.5*(J[2]*Mp+J[3]*Mm)+J[1]*Mz

    Sp,Sm,Sz=Si(N-3,N)
    #H+=im*mu*Sz
    Sp,Sm,Sz=Si(N-2,N)
    H+=-im*mu*Sz
    Sp,Sm,Sz=Si(N-1,N)
    H+=-im*mu*Sz
    Sp,Sm,Sz=Si(N,N)
    #H+=im*mu*Sz

    for i=1:N
      Sp,Sm,Sz=Si(i,N)
      H+=B*Sz
    end
    
    if BC=="PBC"
      Sp_1,Sm_1,Sz_1=Si(N,N)
      Sp_2,Sm_2,Sz_2=Si(1,N)
      
      H += J[2]/2*Sp_1*Sm_2
      H += J[3]/2*Sm_1*Sp_2
      H += J[1]*Sz_1*Sz_2
    end

    return H
  end    
end

function get_SSH_spin_chain(J1,J2;N=10,mode="MPO")      # spin chain with dimerized J1,J2 NN exchange; N is chain length; mode switch between MPO and Tensor representation.
  if mode=="MPO"
    sites = siteinds("S=1/2",N)
    h = OpSum()
    for i=1:Int(N/2)-1
      h += 0.5*J1,"S+",2*i-1,"S-",2*i
      h += 0.5*J1,"S-",2*i-1,"S+",2*i
      h += J1,"Sz",2*i-1,"Sz",2*i

      h += 0.5*J2,"S+",2*i,"S-",2*i+1
      h += 0.5*J2,"S-",2*i,"S-",2*i+1
      h += J2,"Sz",2*i,"Sz",2*i+1
    end
    h += 0.5*J1,"S+",N-1,"S-",N
    h += 0.5*J1,"S-",N-1,"S+",N
    h += J1,"Sz",N-1,"Sz",N
    H=MPO(h,sites)
    return H
  else
    dim=2^N
    H=spzeros(dim,dim)
    for i=1:Int(N/2)-1
      Mp,Mm,Mz=Si_Sj(2*i-1,2*i,N)
      H+=0.5*J1*(Mp+Mm)+J1*Mz

      Mp,Mm,Mz=Si_Sj(2*i,2*i+1,N)
      H+=0.5*J2*(Mp+Mm)+J2*Mz
    end
    Mp,Mm,Mz=Si_Sj(N-1,N,N)
    H+=0.5*J1*(Mp+Mm)+J1*Mz
    return H
  end  
end

function get_Sp_MPO(sites,i)
  sp=OpSum()
  sp+=1,"S+",i
  Sp=MPO(sp,sites)
  return Sp
end

function get_Sz_MPO(sites,i)
  sp=OpSum()
  sp+=1,"Sz",i
  Sp=MPO(sp,sites)
  return Sp
end

function get_Sm_MPO(sites,i)
  sp=OpSum()
  sp+=1,"S-",i
  Sp=MPO(sp,sites)
  return Sp
end

function get_Sx_MPO(sites,i)
  sx=OpSum()
  sx+=1,"Sx",i
  Sx=MPO(sx,sites)
  return Sx
end

function get_Sy_MPO(sites,i)
  sy=OpSum()
  sy+=1,"Sy",i
  Sy=MPO(sy,sites)
  return Sy
end

function multiply_JW_string(psi,i;mode="sparse_tensor",sites=0,N=0)
  if i==1
    return psi
  else
    if mode=="MPO"
      for j=1:i-1
        F=-2*get_Sz_MPO(sites,j)
        psi=contract(F,psi)
      end
      return psi
    else
      for j=1:i-1
        Sp,Sm,Sz=Si(j,N,spin=1/2,mode=mode)
        psi=-2*Sz*psi
      end
      return psi
    end      
  end
end


function test_spin_H()
  N=2
  J=[0,im,-im,0,0,0]
  B=[0,0,0,0,0,0]
  H=get_generalized_SSH_spin_chain(J,B,N=N,mode="tensor")
  println(H)
  println(eigen(H).values)
  println(minimum(eigen(H).values))
  psi=eigen(H).vectors[:,1]
  sp,sm,sz=Si(1,N;spin=1/2)
  println(psi'*sz*psi)
  
  H2=get_generalized_SSH_spin_chain(J,B,N=N,mode="MPO")
  sites=firstsiteinds(H2)
  sweeps = Sweeps(10) # number of sweeps is 5
  maxdim!(sweeps,10,10,10,10,10) # gradually increase states kept
  cutoff!(sweeps,1E-10) # desired truncation error

  psi0 = randomMPS(sites,10)
  psi0=psi0*2

  energy,psi = dmrg(H2,psi0,sweeps)
  println(energy)
  println(psi)
  psi1=contract(H2,psi,cutoff=1e-3,maxdim=10)
  println(inner(psi,psi1))
  println(expect(psi,"Sz"))
end

function test_NH_spin_H()
  N=2
  J=[0,1,2,0,0,0]
  B=[0,0,0,0,0,0]
  H=get_generalized_SSH_spin_chain(J,B,N=N,mode="MPO")
  sites=firstsiteinds(H)

  H2=get_generalized_SSH_spin_chain(J,B,N=N,mode="tensor")
  A1=[1 0]
  A2=[0 1] 
  A=kron(A1,A2)
  A=A/norm(A)
  println(H[1]*H[2])
  println(H2)
  println(A,"\n")
  #println(H2*A')
  psi = MPS(A,sites;cutoff=1e-8,maxdim=2)
  println(psi[1]*psi[2])
  psi1=*(H,psi)
  #println(psi1[1]*psi1[2])

  Id=get_Id_MPO(sites)
  #println(Id[1]*Id[2])
  psi2=*(Id,psi)
  #println(psi2[1]*psi2[2])
  
  psi3=+(psi1,psi2,cutoff=1e-10,maxdim=10)
  #println(psi3[1]*psi3[2])
end
