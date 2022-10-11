struct Population
    zetas::Array{Complex{Float64},1}
    couplings::Array{Float64,1}
end


function init_population(Np::Int64, c::Int64)
    #dist_energy = Exponential()
    dist_energy = Normal(0, 1/sqrt(c))
    Population(rand(Complex{Float64},Np), rand(dist_energy, Np))
end

function generate_population(lambda::Float64,c::Int64, alpha::Float64,Np::Int64, nsteps::Int64, epsilon::Float64)
    poparray = init_population(Np, c)
    dist_energy = Normal(0, 1/sqrt(c))
    p = alpha/(c+ alpha);
    dist_degree = NegativeBinomial(a+1, p)
    for i in 1:nsteps
        k = rand(dist_degree) + 1  ### k ~ k/c*p_k
        random_elements = rand(1:Np, k-1)
        if length(random_elements) > 0
            zetas = poparray.zetas[random_elements]
            couplings = poparray.couplings[random_elements]
            partial_sum = sum(couplings.^2 .* zetas)
            zeta_cm1 = 1/((lambda-epsilon*im) - partial_sum)
            relement = rand(1:Np)
            e1 = rand(dist_energy)
            poparray.zetas[relement] = zeta_cm1
            poparray.couplings[relement] = e1
        end
    end
    poparray
end

function population_update!(lambda::Float64,c::Int64, alpha::Float64, Np::Int64, nsteps::Int64, epsilon::Float64, poparray)
    dist_energy = Normal(0, 1/sqrt(c))
    p = alpha/(c+ alpha);
    dist_degree = NegativeBinomial(alpha+1, p)
    for i in 1:nsteps
        k = rand(dist_degree) + 1  ### k ~ k/c*p_k
        random_elements = rand(1:Np, k -1)
        if length(random_elements) > 0
            zetas = poparray.zetas[random_elements]
            couplings = poparray.couplings[random_elements]
            partial_sum = sum(couplings.^2 .* zetas)
            zeta_cm1 = 1/((lambda-epsilon*im) - partial_sum)
            relement = rand(1:Np)
            e1 = rand(dist_energy)
            poparray.zetas[relement] = zeta_cm1
            poparray.couplings[relement] = e1
        end
    end
    poparray
end

function DOS(lambda::Float64, c::Int64, alpha::Float64, Np::Int64, epsilon ::Float64; nsteps = Np*10^3, ensemble = Np*10^2)
    poparray = init_population(Np, c)
    population_update!(lambda, c, alpha, Np, nsteps, epsilon, poparray);

    dist_energy = Normal(0, 1/sqrt(c))
    p = alpha/(c+ alpha);
    dist_degree = NegativeBinomial(alpha, p)
  
    res = zeros(ensemble)

    k = rand(dist_degree)
    random_elements = rand(1:Np,k)
    partial_sum = 0.0

    if length(random_elements) > 0
        zetas = poparray.zetas[random_elements]
        couplings = poparray.couplings[random_elements]  
        partial_sum = sum(couplings.^2 .* zetas)
    end

    zeta_c = 1/((lambda-epsilon*im) - partial_sum)
    res[1] = imag(zeta_c)

    #population_update!(lambda,c,T,Np, 1, epsilon, poparray)

    for j in 2:ensemble
        k = rand(dist_degree)
        random_elements = rand(1:Np,k)
        partial_sum = 0.0
        
        if length(random_elements) > 0
            zetas = poparray.zetas[random_elements]
            couplings = poparray.couplings[random_elements]  
            partial_sum = sum(couplings.^2 .* zetas)
        end

        zeta_c = 1/((lambda-epsilon*im) - partial_sum)
        res[j] = imag(zeta_c)
        ##update the population with new values
        #population_update!(lambda,c,T,Np, 1, epsilon, poparray)
        ##########################
    end
    
    mean(res)*1/pi
    
end


function resolvent(lambda::Float64, c::Int64, alpha::Float64, Np::Int64, epsilon ::Float64; nsteps = Np*10^3, ensemble = Np*10^2)
    
    poparray = init_population(Np, c)
    population_update!(lambda, c, alpha, Np, nsteps, epsilon, poparray);

    dist_energy = Normal(0, 1/sqrt(c))
    p = alpha/(c+ alpha);
    dist_degree = NegativeBinomial(alpha, p)
  
    res = zeros(ensemble)

    k = rand(dist_degree)
    random_elements = rand(1:Np,k)
    partial_sum = 0.0

    if length(random_elements) > 0
        zetas = poparray.zetas[random_elements]
        couplings = poparray.couplings[random_elements]  
        partial_sum = sum(couplings.^2 .* zetas)
    end

    zeta_c = 1/((lambda-epsilon*im) - partial_sum)
    res[1] = imag(zeta_c)

    population_update!(lambda,c, alpha,Np, 1, epsilon, poparray)

    for j in 2:ensemble
        k = rand(dist_degree)
        random_elements = rand(1:Np,k)
        partial_sum = 0.0
        
        if length(random_elements) > 0
            zetas = poparray.zetas[random_elements]
            couplings = poparray.couplings[random_elements]  
            partial_sum = sum(couplings.^2 .* zetas)
        end

        zeta_c = 1/((lambda-epsilon*im) - partial_sum)
        res[j] = imag(zeta_c)
        ##update the population with new values
        population_update!(lambda,c,alpha,Np, 1, epsilon, poparray)
        ##########################
    end
    
    #mean(res)*1/pi
    res
end

