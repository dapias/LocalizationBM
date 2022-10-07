using Distributions, LightGraphs, LinearAlgebra, SparseArrays, StatsBase

###Generates the instances
function generate_sparse_barrat_matrix(n::Int64,c::Int64)
    dene = Exponential()
    L = random_regular_graph(n,c)
    adj = adjacency_matrix(L);
    energies = rand(dene, n);
    #M = 1.0*adj;
    list_neighbours = LightGraphs.SimpleGraphs.fadj(L)

    return energies, list_neighbours, adj, L
end

function barrat_matrix(energies::Array, list::Array, adj::SparseMatrixCSC, n::Int64, c::Int64, T::Float64)
    beta = 1/T
    M = 1.0*adj;

    a = zeros(n)
    
    for i in 1:n
        for j in 1:c
            k = list[i][j]
            M[i,k] = 1.0/c*(1.0/(1.0+exp(-beta*(energies[i] - energies[k]))))
            a[k] += M[i,k]
        end
    end
    
    for i in 1:n
        M[i,i] = -a[i]
    end
    
    peq = exp.(beta*energies)
    peq /= sum(peq)
    p_inv = peq.^(-1/2)
    p_dir = peq.^(1/2)
    Ms = Diagonal(p_inv)*M *Diagonal(p_dir);

   
    #    return  Symmetric(Ms)
    return  Symmetric(Matrix(Ms))
end


function omega_generator(adj::SparseMatrixCSC{Int64,Int64}, n::Int64,
                         c::Int64)
    col = repeat(1:n, inner = c)
    sparse(adj.rowval, col, rand(Complex{Float64}, c*n), n, n)
end


function symmetric_f(e1::Float64, beta::Float64, e2::Float64)
    exp(beta*(e1 + e2)/2)/(2*cosh(beta*(e1-e2)/2.)) 
end

function local_error(lambda::Float64,  neighs::Array{Array{Int64,1},1},
        Omega_old::SparseMatrixCSC{Complex{Float64},Int64}, epsilon::Float64, 
        beta::Float64, energies::Array{Float64,1}, N::Int64, 
        Omega_sum::SparseMatrixCSC{Complex{Float64},Int64})

    for k in 1:N
        value = neighs[k]
        for j in value
            a = 0.0+0.0*im
            for l in value
                if l != j
                  fs = symmetric_f(energies[k], beta, energies[l])
                  a += im*fs*Omega_old[l, k]/(im*fs +Omega_old[l,k])
                end
            end
            @inbounds Omega_sum[k,j] = a + im*(lambda - im*epsilon)/(exp(-beta*energies[k])/c)
        end
    end
    
    error = sum(abs.(Omega_sum .- Omega_old))  

#    @inbounds fill!(Omega_old, zero(Complex{Float64}))

    return error, Omega_sum, Omega_old
end

function fixed_point_bp(lambda::Float64, c::Int64, epsilon::Float64, energies::
                        Array{Float64,1}, nei::Array{Array{Int64,1},1}, n::Int64, 
        adj::SparseMatrixCSC{Int64,Int64}, T::Float64; tolerance = 0.1)

    beta  = 1/T
    error2 = 10.0*tolerance*n*n
    Omegas = omega_generator(adj, n, c)
    #Omegap = spzeros(Complex{Float64},n,n)
    Omegap = similar(Omegas)
    
    
    while error2 > tolerance
         error2, Omegas, Omegap = local_error(lambda, nei,  Omegas, epsilon, beta,  energies, n, Omegap)
    end
    
    return Omegas
    
end


function dens(lambda::Float64, c::Int64, epsilon::Float64, energies::
              Array{Float64,1}, nei::Array{Array{Int64,1},1}, n::Int64, Omegas::SparseMatrixCSC{Complex{Float64},Int64}, T::Float64 )
    beta = 1/T
    sum_var = 0.
    for j in 1:n
        sum_j = 0.
        for k in nei[j]
            fs = symmetric_f(energies[k], beta, energies[j])
            @inbounds sum_j += im*Omegas[k,j]*fs/( im*fs + Omegas[k,j])
        end
        
        omega_j = im*(lambda - epsilon*im)/(exp(-beta*energies[j])/c) + sum_j
        sum_var +=  real( (1/(omega_j*exp(-beta*energies[j])/c)))
    end
    
    return 1/(pi*n)*sum_var
end


function local_dos(lambda::Float64, c::Int64, epsilon::Float64, energies::
             Array{Float64,1}, nei::Array{Array{Int64,1},1}, n::Int64, Omegas::SparseMatrixCSC{Complex{Float64},Int64} , T::Float64)
    beta = 1/T
    ldos = zeros(n)
    #sum_ipr = 0.
    for j in 1:n
        sum_j = 0.
        for k in nei[j]
            fs = symmetric_f(energies[k], beta, energies[j])
            @inbounds sum_j += im*Omegas[k,j]*fs/( im*fs + Omegas[k,j])
        end
        
        omega_j = im*(lambda - epsilon*im)/(exp(-beta*energies[j])/c) + sum_j
        ldos[j] = real( (1/(omega_j*exp(-beta*energies[j])/c)))
            #sum_var +=  real( (1/(omega_j*exp(-beta*energies[j])/c)))
            #sum_ipr += abs2(1/(omega_j*exp(-beta*energies[j])/c))
    end

    #rho =  1/(pi*n)*sum_var
    #ipr = epsilon/(pi*rho*n)*sum_ipr
    1/pi*ldos

end




function ipr(lambda::Float64, c::Int64, epsilon::Float64, energies::
             Array{Float64,1}, nei::Array{Array{Int64,1},1}, n::Int64, Omegas::SparseMatrixCSC{Complex{Float64},Int64} , T::Float64)
    beta = 1/T
    sum_var = 0.
    sum_ipr = 0.
    for j in 1:n
        sum_j = 0.
        for k in nei[j]
            fs = symmetric_f(energies[k], beta, energies[j])
            @inbounds sum_j += im*Omegas[k,j]*fs/( im*fs + Omegas[k,j])
        end
        
        omega_j = im*(lambda - epsilon*im)/(exp(-beta*energies[j])/c) + sum_j
        sum_var +=  real( (1/(omega_j*exp(-beta*energies[j])/c)))
        sum_ipr += abs2(1/(omega_j*exp(-beta*energies[j])/c))
    end

    rho =  1/(pi*n)*sum_var
    ipr = epsilon/(pi*rho*n)*sum_ipr
end


function all_omegas(lambda::Float64, c::Int64, epsilon::Float64, energies::
             Array{Float64,1}, nei::Array{Array{Int64,1},1}, n::Int64, Omegas::SparseMatrixCSC{Complex{Float64},Int64} , T::Float64)
    beta = 1/T
    ldos = zeros(Complex{Float64}, n)
    #sum_ipr = 0.
    for j in 1:n
        sum_j = 0.
        for k in nei[j]
            fs = symmetric_f(energies[k], beta, energies[j])
            @inbounds sum_j += im*Omegas[k,j]*fs/( im*fs + Omegas[k,j])
        end
        
        omega_j = im*(lambda - epsilon*im)/(exp(-beta*energies[j])/c) + sum_j
        ldos[j] =  1/(omega_j*exp(-beta*energies[j])/c)  ##real part gives the LDOS whereas abs2 leads to the IPR
    end
    ldos
end



function ipr_tikhonov(lambda::Float64, c::Int64, epsilon::Float64, energies::
             Array{Float64,1}, nei::Array{Array{Int64,1},1}, n::Int64, Omegas::SparseMatrixCSC{Complex{Float64},Int64} , T::Float64)
    beta = 1/T
    sum_var = 0.
    sum_ipr = 0.
    for j in 1:n
        sum_j = 0.
        for k in nei[j]
            fs = symmetric_f(energies[k], beta, energies[j])
            @inbounds sum_j += im*Omegas[k,j]*fs/( im*fs + Omegas[k,j])
        end
        
        omega_j = im*(lambda - epsilon*im)/(exp(-beta*energies[j])/c) + sum_j
        sum_var +=  real( (1/(omega_j*exp(-beta*energies[j])/c)))
        sum_ipr += real(1/(omega_j*exp(-beta*energies[j])/c))^2
    end

    rho_square =  (sum_var/n)^2
    ipr = 3/n*(sum_ipr/n)/rho_square
end


# function ipr_biroli(lambda::Float64, c::Int64, epsilon::Float64, energies::
#              Array{Float64,1}, nei::Array{Array{Int64,1},1}, n::Int64, Omegas::SparseMatrixCSC{Complex{Float64},Int64} , T::Float64)
#     beta = 1/T
#     sum_var = 0.
#     sum_ipr = 0.
#     for j in 1:n
#         sum_j = 0.
#         for k in nei[j]
#             fs = symmetric_f(energies[k], beta, energies[j])
#             @inbounds sum_j += im*Omegas[k,j]*fs/( im*fs + Omegas[k,j])
#         end
        
#         omega_j = im*(lambda - epsilon*im)/(exp(-beta*energies[j])/c) + sum_j
#         sum_var +=  real( (1/(omega_j*exp(-beta*energies[j])/c)))
#         sum_ipr += abs2(1/(omega_j*exp(-beta*energies[j])/c))
#     end

#     rho_square =  (sum_var/(pi*n))^2
#     ipr = 3/n*(sum_ipr/n)/rho_square
# end

function biroli(lambda::Float64, c::Int64, epsilon::Float64, energies::
             Array{Float64,1}, nei::Array{Array{Int64,1},1},
        n::Int64, Omegas::SparseMatrixCSC{Complex{Float64},Int64} , T::Float64, 
    pref)
    ##epsilon = pref/(N*ρ(λ))  By default pref = 3*pi
    beta = 1/T
    sum_var = 0.
    sum_ipr = 0.
    for j in 1:n
        sum_j = 0.
        for k in nei[j]
            fs = symmetric_f(energies[k], beta, energies[j])
            @inbounds sum_j += im*Omegas[k,j]*fs/( im*fs + Omegas[k,j])
        end
        
        omega_j = im*(lambda - epsilon*im)/(exp(-beta*energies[j])/c) + sum_j
        sum_var +=  real( (1/(omega_j*exp(-beta*energies[j])/c)))
        sum_ipr += abs2(1/(omega_j*exp(-beta*energies[j])/c))
    end

    rho_square =  (sum_var/(pi*n))^2
    ipr = pref/(pi*n)*(sum_ipr/n)/rho_square
end
