using JLD, Arpack, SparseArrays, Distributions, Graphs, LinearAlgebra

include("../../SparseBarratMezard/src/single_instance.jl")


c = 3 ;
lambda = -1/3

Ts = collect(0.3:0.1:1.2)
ns = [2^11, 2^12, 2^13, 2^14, 2^15]
nsims = Int64.(2^19 ./ns)
nvectors = 10  ###Number of vectors per instance taken around lambda = 0.0

entropies = zeros(length(ns))
iprs = zeros(length(ns))
std1s = zeros(length(ns))
std2s = zeros(length(ns))
for r in 1:length(Ts)
    T = Ts[r]
    for k in 1:length(ns)
        ips = zeros(nsims[k])
        ss = zeros(nsims[k])
        N = ns[k]
        for nrun in 1:nsims[k]
            println("T =$(T), n=$(N), nrun=$(nrun)")
            enes, list, adja, graph = generate_sparse_barrat_elements(N, c);
            Ms = generate_sparse_barrat_matrix(enes, list, adja,N, c, T);
            a = eigs(Ms, sigma = lambda, nev = nvectors)
            evals = real.(a[1]);
            evecs = real.(a[2]);
            mean_s = -mean([sum(evecs[:,i].^2.0.*log.(evecs[:,i].^2)) for i in 1:length(evals)])
            mean_i2 = mean([sum(evecs[:,i].^4) for i in 1:length(evals)])
            ips[nrun] = mean_i2
            ss[nrun] = mean_s
        end
        entropies[k] = mean(ss)
        iprs[k] = exp(mean(log.(ips)))
        std1s[k] = std(ss)
        std2s[k]= exp(std(log.(ips)))
        save("../data/fix_T=$(T).jld", "ipr", iprs, "ns", ns, "stds_ipr", std2s, "entropies", entropies)
    end
end
        




