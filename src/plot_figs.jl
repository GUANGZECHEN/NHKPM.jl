using PyCall
plt = pyimport("matplotlib.pyplot")
np = pyimport("numpy")

function plot_DOS_ED(h,n_Ex,Ex_max,n_Ey,Ey_max,eta)
  N=n_Ex*n_Ey    
  Exs=zeros(N)
  Eys=zeros(N)   # E=Ex+iEy
  rhos=zeros(N)

  for i=1:n_Ex
    for j=1:n_Ey
      ii = Int((j-1)+(i-1)*(n_Ey)+1)
      Exs[ii]=-Ex_max+2*Ex_max*(i-1)/(n_Ex-1)
      Eys[ii]=-Ey_max+2*Ey_max*(j-1)/(n_Ey-1)
      z=Exs[ii]+im*Eys[ii]
      rhos[ii]=dos_ED_NH(h,z,eta) 
    end
  end
  
  Exs_p=zeros(n_Ex)
  Eys_p=zeros(n_Ey)
  rhos_p=zeros((n_Ey,n_Ex))

  for i=1:n_Ex
    Exs_p[i]=-Ex_max+2*Ex_max*(i-1)/(n_Ex-1)
  end

  for j=1:n_Ey
    Eys_p[j]=-Ey_max+2*Ey_max*(j-1)/(n_Ey-1)
  end

  for i=1:n_Ex
    for j=1:n_Ey
      ii = Int((j-1)+(i-1)*(n_Ey)+1)
      rhos_p[j,i]=rhos[ii]
    end
  end

  X,Y=np.meshgrid(Exs_p,Eys_p)
  #println(size(X),size(Y),size(rhos))
  plt.pcolormesh(X,Y,rhos_p)  # find a better plot algorithm
  plt.show()

  return Exs, Eys, rhos
end

function plot_DOS_KPM(h,n_Ex,Ex_max,n_Ey,Ey_max;E_max=10,npol=100,kernel="Jackson")
  N=n_Ex*n_Ey    
  Exs=zeros(N)
  Eys=zeros(N)   # E=Ex+iEy
  rhos=zeros(N)

  for i=1:n_Ex
    for j=1:n_Ey
      ii = Int((j-1)+(i-1)*(n_Ey)+1)
      Exs[ii]=-Ex_max+2*Ex_max*(i-1)/(n_Ex-1)
      Eys[ii]=-Ey_max+2*Ey_max*(j-1)/(n_Ey-1)
      z=Exs[ii]+im*Eys[ii]
      rhos[ii]=abs(dos_kpm_NH(h,z,E_max=E_max,N=npol,kernel=kernel))
    end
  end
  
  Exs_p=zeros(n_Ex)
  Eys_p=zeros(n_Ey)
  rhos_p=zeros((n_Ey,n_Ex))

  for i=1:n_Ex
    Exs_p[i]=-Ex_max+2*Ex_max*(i-1)/(n_Ex-1)
  end

  for j=1:n_Ey
    Eys_p[j]=-Ey_max+2*Ey_max*(j-1)/(n_Ey-1)
  end

  for i=1:n_Ex
    for j=1:n_Ey
      ii = Int((j-1)+(i-1)*(n_Ey)+1)
      rhos_p[j,i]=rhos[ii]
    end
  end

  X,Y=np.meshgrid(Exs_p,Eys_p)
  #println(size(X),size(Y),size(rhos))
  plt.pcolormesh(X,Y,rhos_p)  # find a better plot algorithm
  plt.show()

  return Exs, Eys, rhos
end

function plot_LDOS_KPM(h,n_Ex,Ex_max,n_Ey,Ey_max;E_max=10,npol=100,kernel="Jackson")
  N=size(h)[1]  # number of sites
  NN=n_Ex*n_Ey*N
  
  # all DOS  
  sites=zeros(NN)
  Exs=zeros(NN)
  Eys=zeros(NN)
  rhos=zeros(NN)
  
  # tDOS
  N_tDOS=n_Ex*n_Ey
  tDOSs=zeros(N_tDOS)
  Exs_tDOS=zeros(N_tDOS)
  Eys_tDOS=zeros(N_tDOS)
  
  # LDOS_real
  N_LDOS=n_Ex*N
  LDOS_real=zeros(N_LDOS)
  Es_LDOS_real=zeros(N_LDOS)
  sites_LDOS_real=zeros(N_LDOS)

  for i=0:n_Ex-1
    for j=0:n_Ey-1
      Ex=-Ex_max+2*Ex_max*i/(n_Ex-1)
      Ey=-Ey_max+2*Ey_max*j/(n_Ey-1)
      z=Ex+im*Ey
      local tDOS=0
      
      for site_index=1:N
        ii=site_index+N*(j+i*n_Ey)
        rho=ldos_kpm_NH(h,z,site_index;E_max=E_max,N=npol,kernel=kernel)
        
        Exs[ii]=Ex
        Eys[ii]=Ey
        sites[ii]=site_index
        rhos[ii]=real(rho)
        
        tDOS+=rho        
        
        Es_LDOS_real[site_index+i*N]=Ex
        sites_LDOS_real[site_index+i*N]=site_index
        LDOS_real[site_index+i*N]+=real(rho)
      end
      Exs_tDOS[j+i*n_Ey+1]=Ex
      Eys_tDOS[j+i*n_Ey+1]=Ey
      tDOSs[j+i*n_Ey+1]=abs(tDOS)      
    end
  end
  
  Exs_p=zeros(n_Ex)
  sites_p=zeros(N)
  rhos_p=zeros((n_Ex,N))

  for i=1:n_Ex
    Exs_p[i]=-Ex_max+2*Ex_max*(i-1)/(n_Ex-1)
  end

  for j=1:N
    sites_p[j]=j
  end

  for i=1:n_Ex
    for j=1:N
      ii = Int((j-1)+(i-1)*N+1)
      rhos_p[i,j]=LDOS_real[ii]
    end
  end

  X,Y=np.meshgrid(sites_p,Exs_p)
  #println(size(X),size(Y),size(rhos))
  plt.pcolormesh(X,Y,rhos_p)  # find a better plot algorithm
  plt.show()
  
  open("all_DOS.OUT","w") do io
    writedlm(io, [Exs Eys sites rhos])
  end

  open("tDOS.OUT","w") do io
    writedlm(io, [Exs_tDOS Eys_tDOS tDOSs])
  end

  open("LDOS_real_projected.OUT","w") do io
    writedlm(io, [Es_LDOS_real sites_LDOS_real LDOS_real])
  end
end



function plot_eigen(h;label=" ",mode="eigen")
  if mode=="eigen"
    F=eigen(h)
    eigvalues,eigenvecs=F.values,F.vectors
    Exs=real(eigvalues)
    Eys=imag(eigvalues)
    plt.scatter(Exs,Eys,label=label)
  elseif mode=="svd"
    svd=svdvals(h)
    Exs=real(svd)
    N=size(Exs)[1]
    Eys=ones(N)
    if label=="PBC"
      Eys=-Eys
    end
    plt.scatter(Exs,Eys,label=label)
  else
    println("compute either eigen or svd")
  end
end



