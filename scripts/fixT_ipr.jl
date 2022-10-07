using SpecialFunctions, Distributions, LinearAlgebra
using QuadGK, PyCall, QuadGK, Roots, JLD, DelimitedFiles

nx = pyimport("networkx")

function conf_model(N, α, c)
    s=0
  #  c=(N)*0.05
  #  c= N/2
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
    
    a.A
    

end

function weights(N, α, c)
   # c=(N)*0.05
    #c= N/2
    A = conf_model(N, α, c)
    J = rand(Normal(0, 1/sqrt(c)), N,N)
    M = zeros(N,N);
    for i in 1:N
        for j in 1:i-1
            M[i,j] = A[i,j]*J[i,j]
            M[j,i] = M[i,j]
        end
    end
    M
end
    


epsilon = 1e-4;
z = 0.0 .- epsilon*im;

new_f(g,z,a) = g^2-z^2/(-a*exp(-a*g)*expint(1+a,-a*g))
dens(g, z, a) = 1/pi*imag(exp(-a*g)*(-a*g)^a*gamma(1-a, -g*a)/z);

N = 2^10
a = 0.2;
c = N*0.05
M = weights(N, a, c);
esystem = eigen(M);

as = [0.2, 0.3, 0.4, 0.6,0.7, 0.8]
as = [1.2, 1.3, 1.4, 1.6, 1.7, 1.8]
#as = 0.5
ns = [2^10,2^11, 2^12, 2^13, 2^14]
#ns = [2^10,2^11, 2^12]
nsims = Int64.(2^17 ./ns)



epsilon = 1e-4;
z = 0.0 .- epsilon*im;
rhos = zeros(length(as))

for k in 1:length(as)
    a = as[k]
    rho_e = 0.0
    while rho_e <= 0
        x0 = rand(Uniform(-1,1)) + im*rand(Uniform(-1,1))
        dd = Roots.muller(x-> new_f(x,z, a), x0)
        rho_e = dens(dd, z, a)
        rhos[k] = rho_e
    end
end

iprs = zeros(length(ns))
typ_iprs = zeros(length(ns))
stds = zeros(length(ns))
for r in 1:length(as)
    a = as[r]
    rho_e = rhos[r]
    for k in 1:length(ns)
        ips = zeros(nsims[k])
        N = ns[k]
        c = N*0.05
         epsilon = 3/(N*pi*rho_e)
        for nrun in 1:nsims[k]
            println("a=$(a), n=$(N), nrun=$(nrun)")
            M = weights(N, a, c);
            esystem = eigen(M);
            mid_els = findall(x-> -epsilon < x< epsilon, esystem.values);
            evecs = esystem.vectors[:, mid_els];
            mean_i2 = mean([sum(evecs[:,i].^4) for i in 1:length(mid_els)])
            ips[nrun] = mean_i2
        end
        #iprs[r, k] = mean(ips[findall(x->!isnan(x), ips)])
        #typ_iprs[r,k] = exp(mean(log.(ips[findall(x->!isnan(x), ips)])))
        #stds[r,k] = std(ips[findall(x->!isnan(x), ips)])
        #save("../data/fix_a=$(a)_doble.jld", "ipr", iprs, "ns", ns, "stds", stds)
        iprs[k] = mean(ips[findall(x->!isnan(x), ips)])
        typ_iprs[k] = exp(mean(log.(ips[findall(x->!isnan(x), ips)])))
        stds[k] = std(ips[findall(x->!isnan(x), ips)])
        save("../data/fix_a=$(a).jld", "ipr", iprs, "ns", ns, "stds", stds, "typ_ipr", typ_iprs)
    end
end
        



