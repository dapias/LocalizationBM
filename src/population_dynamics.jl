struct Population
    zetas::Array{Complex{Float64},1}
    energies::Array{Float64,1}
end

function symmetric_f(e1::Float64, beta::Float64, e2::Array{Float64})
    exp.(beta*(e1 .+ e2)/2)./(2*cosh.(beta*(e1.-e2)./2.)) 
end

function init_population(Np::Int64)
    dist_energy = Exponential()
    Population(rand(Complex{Float64},Np), rand(dist_energy, Np))
end

function generate_population(lambda::Float64,c::Int64,T::Float64,Np::Int64, nsteps::Int64, epsilon::Float64)
    beta = 1.0/T
    poparray = init_population(Np)
    dist_energy = Exponential()
    for i in 1:nsteps
        random_elements = rand(1:Np,c-1)
        zetas = poparray.zetas[random_elements]
        energies = poparray.energies[random_elements]
        e1 = rand(dist_energy)
        fsym = symmetric_f(e1, beta, energies)
        zeta_cm1 = sum(im*fsym.*zetas./(im*fsym .+ zetas))
        zeta_cm1 += im*(lambda-epsilon*im)*exp(beta*e1)*c 
        relement = rand(1:Np)
        poparray.zetas[relement] = zeta_cm1
        poparray.energies[relement] = e1
    end
    poparray
end

function population_update!(lambda::Float64,c::Int64,T::Float64,Np::Int64, nsteps::Int64, epsilon::Float64, poparray)
    beta = 1.0/T
    dist_energy = Exponential()
    for i in 1:nsteps
        random_elements = rand(1:Np,c-1)
        zetas = poparray.zetas[random_elements]
        energies = poparray.energies[random_elements]
        e1 = rand(dist_energy)
        fsym = symmetric_f(e1, beta, energies)
        zeta_cm1 = sum(im*fsym.*zetas./(im*fsym .+ zetas))
        zeta_cm1 += im*(lambda-epsilon*im)*exp(beta*e1)*c
        relement = rand(1:Np)
        poparray.zetas[relement] = zeta_cm1
        poparray.energies[relement] = e1
    end
    poparray
end

function DOS(lambda::Float64, c::Int64, T::Float64, Np::Int64, epsilon ::Float64; nsteps = Np*10^3, ensemble = Np*10^2)
    poparray = init_population(Np)
    population_update!(lambda, c, T, Np, nsteps, epsilon, poparray);
    
    beta = 1.0/T
    res = zeros(ensemble)
    dist_energy = Exponential()
    random_elements = rand(1:Np,c)
    zetas_sample = poparray.zetas[random_elements]
    energies = poparray.energies[random_elements]  
    e1 = rand(dist_energy)
    fsym = symmetric_f(e1, beta, energies)
    sum_term = sum(im*fsym.*zetas_sample./(im*fsym .+ zetas_sample))
    zeta_c =  im* (lambda-epsilon*im)*exp(beta*e1)*c  + sum_term
    res[1] = real(c*exp(beta*e1)/zeta_c)
    population_update!(lambda,c,T,Np, 1, epsilon, poparray)

    for j in 2:ensemble
        random_elements = rand(1:Np,c)
        zetas_sample = poparray.zetas[random_elements]
        energies = poparray.energies[random_elements]  
        e1 = rand(dist_energy)
        fsym = symmetric_f(e1, beta, energies)
        sum_term = sum(im*fsym.*zetas_sample./(im*fsym .+ zetas_sample))
        zeta_c =   im*(lambda-epsilon*im)/(exp(-beta*e1)/c)  + sum_term
        res[j] = real(c*exp(beta*e1)/zeta_c)
        ##update the population with new values
        population_update!(lambda,c,T,Np, 1, epsilon, poparray)
        ##########################
    end

    mean(res)*1/pi
    
end

