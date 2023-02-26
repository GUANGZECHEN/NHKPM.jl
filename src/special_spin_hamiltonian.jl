include("spin_hamiltonian.jl")

function get_HdagH_generalized_spin_chain(J;N=8,mode="MPO",boundary="OBC",sites=0)       # not bi-partite
  J2=zeros(Complex,3)
  J2[1:size(J)[1]].=J[:]
  J=J2
  
  if mode=="MPO"
    if sites==0
      sites = siteinds("S=1/2",N)
    end
    h = OpSum()
    for i=1:N-1
      for j=1:N-1
        h += 0.5*J[2]*0.5*J[2],"S-",j,"S+",j+1,"S+",i,"S-",i+1                           # be careful I devide 1/2 here, which means the parameter I put is J_perp
        h += 0.5*J[3]*0.5*J[3],"S+",j,"S-",j+1,"S-",i,"S+",i+1
        h += J[1]*J[1],"Sz",j,"Sz",j+1,"Sz",i,"Sz",i+1
        
        h += 0.5*J[2]*0.5*J[3],"S-",j,"S+",j+1,"S-",i,"S+",i+1
        h += 0.5*J[2]*0.5*J[3],"S+",j,"S-",j+1,"S+",i,"S-",i+1
        
        h += 0.5*J[2]*J[1],"S-",j,"S+",j+1,"Sz",i,"Sz",i+1
        h += 0.5*J[2]*J[1],"Sz",j,"Sz",j+1,"S+",i,"S-",i+1
        
        h += 0.5*J[3]*J[1],"S+",j,"S-",j+1,"Sz",i,"Sz",i+1
        h += 0.5*J[3]*J[1],"Sz",j,"Sz",j+1,"S-",i,"S+",i+1
      end
    end
    
    # add PBC later
    
    H=MPO(h,sites)
    return H
  else
    dim=2^N
    H=spzeros(dim,dim)
    if boundary=="OBC"
      for i=1:N-1
        for j=1:N-1
          Sp_i,Sm_i,Sz_i=Si(i,N,mode=mode)
          Sp_i2,Sm_i2,Sz_i2=Si(i+1,N,mode=mode)
          Sp_j,Sm_j,Sz_j=Si(j,N,mode=mode)
          Sp_j2,Sm_j2,Sz_j2=Si(j+1,N,mode=mode)
                                
          H += 0.5*J[2]*0.5*J[2]*Sm_j*Sp_j2*Sp_i*Sm_i2                           # be careful I devide 1/2 here, which means the parameter I put is J_perp
          H += 0.5*J[3]*0.5*J[3]*Sp_j*Sm_j2*Sm_i*Sp_i2
          H += J[1]*J[1]*Sz_j*Sz_j2*Sz_i*Sz_i2
        
          H += 0.5*J[2]*0.5*J[3]*Sm_j*Sp_j2*Sm_i*Sp_i2
          H += 0.5*J[2]*0.5*J[3]*Sp_j*Sm_j2*Sp_i*Sm_i2
        
          H += 0.5*J[2]*J[1]*Sm_j*Sp_j2*Sz_i*Sz_i2
          H += 0.5*J[2]*J[1]*Sz_j*Sz_j2*Sp_i*Sm_i2
        
          H += 0.5*J[3]*J[1]*Sp_j*Sm_j2*Sz_i*Sz_i2
          H += 0.5*J[3]*J[1]*Sz_j*Sz_j2*Sm_i*Sp_i2
        end
      end
    else
      for i=1:N
        for j=1:N
          i2=i+1
          j2=j+1
          if i==N
            i2=1
          end
          
          if j==N
            j2=1
          end            
        
          Sp_i,Sm_i,Sz_i=Si(i,N,mode=mode)
          Sp_i2,Sm_i2,Sz_i2=Si(i2,N,mode=mode)
          Sp_j,Sm_j,Sz_j=Si(j,N,mode=mode)
          Sp_j2,Sm_j2,Sz_j2=Si(j2,N,mode=mode)
                                
          H += 0.5*J[2]*0.5*J[2]*Sm_j*Sp_j2*Sp_i*Sm_i2                           # be careful I devide 1/2 here, which means the parameter I put is J_perp
          H += 0.5*J[3]*0.5*J[3]*Sp_j*Sm_j2*Sm_i*Sp_i2
          H += J[1]*J[1]*Sz_j*Sz_j2*Sz_i*Sz_i2
        
          H += 0.5*J[2]*0.5*J[3]*Sm_j*Sp_j2*Sm_i*Sp_i2
          H += 0.5*J[2]*0.5*J[3]*Sp_j*Sm_j2*Sp_i*Sm_i2
        
          H += 0.5*J[2]*J[1]*Sm_j*Sp_j2*Sz_i*Sz_i2
          H += 0.5*J[2]*J[1]*Sz_j*Sz_j2*Sp_i*Sm_i2
        
          H += 0.5*J[3]*J[1]*Sp_j*Sm_j2*Sz_i*Sz_i2
          H += 0.5*J[3]*J[1]*Sz_j*Sz_j2*Sm_i*Sp_i2
        end
      end     
    end
    # add PBC later
    
    return H
  end
end

function get_generalized_SSH_spin_chain(J,B;N=10,mode="MPO",boundary="OBC",sites=0)      # generalized dimerized spin chain, J and B are vectors; N is chain length; mode switch between MPO and Tensor representation.
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
      h += 0.5*J[6],"S+",1,"S-",N
      h += 0.5*J[5],"S-",1,"S+",N
      h += J[4],"Sz",1,"Sz",N
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
      H+=0.5*J[5]*Mm+0.5*J[6]*Mp+J[4]*Mz                     # here mp and mm are reversed
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


function get_Liouvillian_spin_chain(Jx,Jy,gamma;Bx=0,N=4,mode="sparse_tensor",boundary="OBC",sites=0)      # Liouvillian model in arXiv:1812.10373
  if mode=="MPO"
    L=2*N
    if sites==0
      sites = siteinds("S=1/2",L)
    end
    h = OpSum()
    for i=1:Int(N/2)-1
      h += im*Jx,"Sx",4*i-3,"Sx",4*i-1                                                  
      h += im*Jy,"Sy",4*i-1,"Sy",4*i+1
      h += -im*Jx,"Sx",4*i-2,"Sx",4*i                                                  
      h += -im*Jy,"Sy",4*i,"Sy",4*i+2
    end

    h+=im*Jx,"Sx",2*N-3,"Sx",2*N-1 
    h+=-im*Jx,"Sx",2*N-2,"Sx",2*N 
    
    for i=1:N
      h+=gamma,"Sz",2*i-1,"Sz",2*i
    end

    H=MPO(h,sites)
    Identity=get_Id_MPO(sites)
    H=+(H,-N/4*gamma*Identity)
    H=4*H

    return H
  else
    L=2*N
    dim=2^L
    H=spzeros(dim,dim)
    for i=1:Int(N/2)-1
      Mx,My,Mz=Si_Sj_xyz(2*i-1,2*i,L,mode=mode)
      Mx2,My2,Mz2=Si_Sj_xyz(2*i,2*i+1,L,mode=mode)
      
      Mx3,My3,Mz3=Si_Sj_xyz(2*i-1+N,2*i+N,L,mode=mode)
      Mx4,My4,Mz4=Si_Sj_xyz(2*i+N,2*i+1+N,L,mode=mode)
      
      H+=-Jx*Mx-Jy*My2+Jx*Mx3+Jy*My4
    end
    
    Mx,My,Mz=Si_Sj_xyz(N-1,N,L,mode=mode)
    Mx3,My3,Mz3=Si_Sj_xyz(2*N-1,2*N,L,mode=mode)
    H+=-Jx*Mx+Jx*Mx3
    
    for i=1:N
      Mx,My,Mz=Si_Sj_xyz(i,2*N+1-i,L,mode=mode)
      H+=im*gamma*Mz
      
      Sp,Sm,Sz=Si(i,L)
      Sx=(Sp+Sm)/2
      
      Sp2,Sm2,Sz2=Si(L+1-i,L)
      Sx2=(Sp2+Sm2)/2
      
      H+=Bx*Sx-Bx*Sx2
    end
    H=-im*H-N/4*gamma*I
    H=4*H
    return H
  end  
end

function get_Liouvillian_GS(N)
  x=[ITensor() for i=1:2*N]
  M=[1 0;0 1]
    
  bonds=[Index(2,"link") for i=1:(2*N-1)]
  for i=1:(N-1)
    bonds[2*i]=Index(1,"single_link")
  end
    
  x[1]=ITensor(M,sites[1],bonds[1])
  x[2*N]=ITensor(M,bonds[2*N-1],sites[2*N])
  for i=2:2*N-1
    x[i]=ITensor(M,bonds[i-1],sites[i],bonds[i])
  end
  psi_1=MPS(x)
  return psi_1
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
