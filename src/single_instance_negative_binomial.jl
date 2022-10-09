using Distributions, LightGraphs, LinearAlgebra, SparseArrays, StatsBase
using SpecialFunctions
using QuadGK, PyCall, QuadGK, Roots, JLD, DelimitedFiles

nx = pyimport("networkx")

function conf_model(N, α, c)
    s=0
    pp =α/(c+α)
    E=zeros(Int64,N);

    for i in 1:N
        E[i] = rand(NegativeBinomial(α,pp))
        s += E[i]
    end
    
    if iseven(s) == false
        k = rand(DiscreteUniform(1,N))
        E[k] += 1
    end
    
    G = nx.configuration_model(E)
    a = nx.adjacency_matrix(G, nodelist = 0:N-1)
    
    a.A  ###Return Adjacency Matrix
end

function weights(N, α, c)
    A = conf_model(N, α, c)
    J = rand(Normal(0, 1/sqrt(c)), N,N)
    M = zeros(N,N);
    for i in 1:N
        for j in 1:i-1
            M[i,j] = A[i,j]*J[i,j]
            M[j,i] = M[i,j]
        end
    end

    indices = findall(x-> x .!= 0.0 , M);
    nei = [[l[2] for l in indices[findall(x-> x[1] == k, indices)]] for k in 1:N];
    
    M, nei, indices
end



function local_error(lambda::Float64,  neighs::Array{Array{Int64,1},1},
                     Omega_old:: Matrix{ComplexF64}, epsilon::Float64, N::Int64, 
                     Omega_sum::Matrix{ComplexF64}, Weights::Matrix{Float64})
    
    for k in 1:N
        value = neighs[k]
        for j in value
            a = 0.0+0.0*im
            for l in value
                if l != j
                    a += Weights[l, k]^2*Omega_old[l, k]
                end
            end
            @inbounds Omega_sum[k,j] = 1/((lambda - im*epsilon) - a)
        end
    end
    
    error = sum(abs.(Omega_sum .- Omega_old))  

    @inbounds fill!(Omega_old, zero(Complex{Float64}))

    return error, Omega_sum, Omega_old
end

function fixed_point_bp(lambda::Float64, epsilon::Float64, nei::Array{Array{Int64,1},1}, indices, Weights::Matrix{Float64} ; tolerance = 0.1)

    N = size(Weights)[1]
    error2 = 10.0*tolerance*N*N
    Omegas = zeros(Complex{Float64}, N, N);

    Omegas[indices] = rand(Complex{Float64}, length(indices));
    Omegap = zeros(Complex{Float64}, N, N);

    while error2 > tolerance
         error2, Omegas, Omegap = local_error(lambda, nei,  Omegas, epsilon, N, Omegap, Weights)
    end
    
    return Omegas
    
end

function dens(lambda::Float64, epsilon::Float64, neighs::Array{Array{Int64,1},1},
        Omegas::Matrix{ComplexF64}, Weights::Matrix{Float64})
    
    n = length(neighs)
  
    sum_var = 0.
    for j in 1:n
        sum_j = 0.
        for k in neighs[j]
            @inbounds sum_j += Weights[k,j]^2*Omegas[k,j]
        end
        
        omega_j = 1/((lambda - epsilon*im) - sum_j)
        sum_var +=  imag(omega_j)
    end
    
    return 1/(pi*n)*sum_var
end

    

function local_dos(lambda::Float64, epsilon::Float64, neighs::Array{Array{Int64,1},1},
        Omegas::Matrix{ComplexF64}, Weights::Matrix{Float64})

    n = length(neighs)
    ldos = zeros(n)
    #sum_ipr = 0.
    for j in 1:n
        sum_j = 0.
        for k in neighs[j]
             @inbounds sum_j += Weights[k,j]^2*Omegas[k,j]
        end
        
        omega_j = 1/((lambda - epsilon*im) - sum_j)
        ldos[j] =  imag(omega_j)
    end

    1/pi*ldos

end


