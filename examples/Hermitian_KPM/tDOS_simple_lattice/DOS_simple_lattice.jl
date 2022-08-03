using PyCall
plt = pyimport("matplotlib.pyplot")

include("../../../src/KPM.jl")

#g=get_geometry("triangular",2,2)
#g=get_geometry("honeycomb",2,2)
g=get_geometry("square",2,2)
H=get_H(g,1)

E_max=5
Es=zeros(101)
rhos_ED=zeros(101)
rhos_Jackson=zeros(101)
rhos_Lorentz=zeros(101)
for i=0:100
  Es[i+1]=real(-E_max+2*E_max*i/100)
  rhos_ED[i+1]=get_dos_ED(H,Es[i+1],nk=50,eta=5e-1)
end

rhos_Jackson=get_dos_kpm(H,Es,nk=50,N=200,kernel="Jackson",E_max=6)
rhos_Lorentz=get_dos_kpm(H,Es,nk=50,N=200,kernel="Lorentz",E_max=6)

plt.figure(figsize=(6,6),dpi=80)
plt.scatter(Es,rhos_Jackson/maximum(rhos_Jackson),color="red",s=10,label="Jackson")
plt.scatter(Es,rhos_Lorentz/maximum(rhos_Lorentz),color="blue",s=10,label="Lorentz")
plt.scatter(Es,rhos_ED/maximum(rhos_ED),color="black",s=10,label="ED")
plt.legend()
plt.show()
