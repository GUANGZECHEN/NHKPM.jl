using LinearAlgebra
using DelimitedFiles
using Statistics

using PyCall
plt = pyimport("matplotlib.pyplot")
np = pyimport("numpy")
# all the functions here work in 2d

struct unit_vector
    x
    y
    xy
    xmy
end

struct geometry
    lattice::String
    sites
    inter_vector::unit_vector
    dimension::Int
end

struct geometry_1D
  sites
  x
end

function get_geometry_1D(N)
  a=[1,0]
  R=Array{Any}(undef,N)
  for i=1:N
    R[i]=i*a
  end
  g=geometry_1D(R,N*a)
end

function get_inter_vector(n,m,a1,a2)
    inter_vector=unit_vector(n*a1,m*a2,m*a2+n*a1,n*a1-m*a2)
    return inter_vector
end

function get_geometry(lattice,n,m)
    a1=[1,0]
    a2=[0,1]
    R=geometry_simple(n,m,a1,a2)
    if lattice=="triangular"
        a1=[1,0]
        a2=[1/2,-sqrt(3)/2]
        R=geometry_simple(n,m,a1,a2)
    elseif lattice=="square"
        a1=[1,0]
        a2=[0,1]
        R=geometry_simple(n,m,a1,a2)
    elseif lattice=="honeycomb"
        a1=[sqrt(3),0]
        a2=[sqrt(3)/2,-3/2]
        a3=[sqrt(3)/2,-1/2]
        R=geometry_bipartite(n,m,a1,a2,a3)
    else
        println("invalid geometry")
    end

    inter_vector=get_inter_vector(n,m,a1,a2)
    dim=2
    
    return geometry(lattice,R,inter_vector,dim)
end

function geometry_simple(n,m,a1,a2)
    R=Array{Any}(undef,n*m)
    ii=1
    for i=0:n-1
        for j=0:m-1
            R[ii]=i*a1+j*a2
            ii+=1
        end
    end
    return R
end

function geometry_bipartite(n,m,a1,a2,a3)
    R=Array{Any}(undef,2*n*m)
    ii=1
    for i=0:n-1
        for j=0:m-1
            for k=0:1
                R[ii]=i*a1+j*a2+k*a3
                ii+=1
            end
        end
    end
    return R
end


function get_lattice_honeycomb_nanoribbon(n,m)
    a1=[sqrt(3),0]
    a2=[0,-3]
    b1=[sqrt(3)/2,-1/2]
    b2=[sqrt(3)/2,-3/2]
    b3=[0,-2]

    R=[]
    for i=0:n-1
        for j=0:m-1
            push!(R,i*a1+j*a2)
            push!(R,i*a1+j*a2+b1)
            push!(R,i*a1+j*a2+b2)
            push!(R,i*a1+j*a2+b3)
        end
    end
    inter_vector=[[0,0],n*a1,-n*a1]
    return R,inter_vector
end

function get_lattice_triangular_nanoribbon(n,m)
    a1=[1,0]
    a2=[0,-sqrt(3)]
    b1=[1/2,-sqrt(3)/2]

    R=[]
    for i=0:n-1
        for j=0:m-1
            push!(R,i*a1+j*a2)
            push!(R,i*a1+j*a2+b1)
        end
    end
    inter_vector=[[0,0],m*a2,-m*a2]
    return R,inter_vector
end

function geometry_triangular(n,a1,a2)
    R=[]
    for i=0:n-1
        for j=0:n-1-i
            push!(R,i*a1+j*a2)
        end
    end
    return R
end

function geometry_zigzag(size)
    a1=[1,0]
    a2=[-1/2,sqrt(3)/2]
    a3=[-1/2,-sqrt(3)/2]
    R=[]
    for i=0:size
        for j=0:size
            if abs(i)+abs(j)<size+2
                push!(R,i*a1+j*a2)
            end
        end

        for l=1:size
            if abs(i)+abs(l)<size+2
                push!(R,i*a1+l*a3)
            end
        end
    end

    for j=1:size
        for l=1:size
            if abs(j)+abs(l)<size+2
                push!(R,l*a3+j*a2)
            end
        end
    end
    return R
end


function nearest(site_a,site_b)
    return 0.1<norm(site_a-site_b)<1.1
end

function next_nearest(site_a,site_b)
    return 1.1<norm(site_a-site_b)<1.9
end

function next_next_nearest(site_a,site_b)
    return 1.9<norm(site_a-site_b)<2.1
end

function plot_R(R)
    N=size(R,1)
    x=[]
    y=[]
    for i=1:N
        push!(x,R[i][1])
        push!(y,R[i][2])
    end
    plt.scatter(x,y,color="red",s=10)
    plt.axis([-40,40,-40,40])
    plt.xlabel("")
    plt.ylabel("")
    plt.show()
end

function add_vacancy(R,vacant_site)
    deleteat!(R,vacant_site)
    return R
end

function plot_DOS(Ex,Ey,rho)
  plt.figure(figsize=(6,6),dpi=80)
  sc1=plt.scatter(Exs,Eys,c=rho)   
  plt.show()
end

#n=2
#m=2
#g=get_geometry("triangular",n,m)
#R,inter_vector=g.sites,g.inter_vector
#println(R,inter_vector)
#plot_R(R)



