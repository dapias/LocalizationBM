using JLD, Arpack

include("../src/single_instance_negative_binomial.jl")


c =500 ;
lambda = 0.0

alphas = collect(0.5:0.2:2.3)
ns = [2^10, 2^11, 2^12, 2^13]
#ns=[2^13]
nsims = Int64.(2^15 ./ns)
nvectors = 10  ###Number of vectors per instance taken around lambda = 0.0

entropies = zeros(length(ns))
iprs = zeros(length(ns))
std1s = zeros(length(ns))
std2s = zeros(length(ns))
for r in 1:length(alphas)
    alpha = alphas[r]
    for k in 1:length(ns)
        ips = zeros(nsims[k])
        ss = zeros(nsims[k])
        N = ns[k]
        for nrun in 1:nsims[k]
            println("alpha=$(alpha), n=$(N), nrun=$(nrun)")
            M, neis, is, adj, nnodes = weights(N, alpha, c)
            evals, evecs = eigs(M, sigma = 1e-300, nev = nvectors)
            mean_s = -mean([sum(evecs[:,i].^2.0.*log.(evecs[:,i].^2)) for i in 1:length(evals)])
            mean_i2 = mean([sum(evecs[:,i].^4) for i in 1:length(evals)])
            ips[nrun] = mean_i2
            ss[nrun] = mean_s
        end
        entropies[k] = mean(ss[findall(x->!isnan(x), ss)])
        iprs[k] = exp(mean(log.(ips)))
        std1s[k] = std(ss[findall(x->!isnan(x), ss)])
        std2s[k]= exp(std(log.(ips)))
        save("../data/fix_alpha=$(alpha).jld", "ipr", iprs, "ns", ns, "stds_ipr", std2s, "entropies", entropies)
    end
end
        




