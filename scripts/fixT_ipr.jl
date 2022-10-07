include("../src/single_instance.jl")
include("../src/population_dynamics.jl")

c =3 ;
Np = 10^4;
lambda = -1/3;
epsilon = 1e-4;  ###To determine ρ_ϵ(λ)

Ts = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
ns = [2^10,2^11, 2^12, 2^13, 2^14]
nsims = Int64.(2^17 ./ns)

rhos = zeros(length(Ts))

for k in 1:length(Ts)
    T = Ts[k]
    rhos[k] = DOS(lambda, c, T, Np, epsilon)
end

iprs = zeros(length(ns))
typ_iprs = zeros(length(ns))
stds = zeros(length(ns))
for r in 1:length(Ts)
    T = Ts[r]
    rho_e = rhos[r]
    for k in 1:length(ns)
        ips = zeros(nsims[k])
        N = ns[k]
        epsilon = 3/(N*pi*rho_e)
        for nrun in 1:nsims[k]
            println("T=$(T), n=$(N), nrun=$(nrun)")
            energies, list_neighbours, adjacency, L = generate_sparse_barrat_matrix(N,c);
            M = barrat_matrix(energies, list_neighbours, adjacency, N, c, T);
            esystem = eigen(M);
            mid_els = findall(x-> -epsilon+ lambda < x< lambda + epsilon, esystem.values);
            evecs = esystem.vectors[:, mid_els];
            mean_i2 = mean([sum(evecs[:,i].^4) for i in 1:length(mid_els)])
            ips[nrun] = mean_i2
        end
        iprs[k] = mean(ips[findall(x->!isnan(x), ips)])
        typ_iprs[k] = exp(mean(log.(ips[findall(x->!isnan(x), ips)])))
        stds[k] = std(ips[findall(x->!isnan(x), ips)])
        save("../data/fix_T=$(T).jld", "ipr", iprs, "ns", ns, "stds", stds, "typ_ipr", typ_iprs)
    end
end
        




