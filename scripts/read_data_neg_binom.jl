using JLD, Plots, LaTeXStrings


alphas = collect(0.5:0.2:2.3)
alphas = collect(0.5:0.2:0.9)


a = load("../data/fix_alpha=$(alphas[1]).jld")

ns = a["ns"]
ents = a["entropies"]
d1 = a["entropies"]./(log.(ns))
ipr = a["ipr"]
d2 = -log.(ipr)./log.(ns)



for k in 2:length(alphas)
    a = load("../data/fix_alpha=$(alphas[k]).jld")
    ents = a["entropies"]
    global d1 = hcat(d1, a["entropies"]./(log.(ns)))
    ipr = a["ipr"]
    global d2 = hcat(d2, -log.(ipr)./log.(ns))
end

plot(alphas, d1[1,:], color = :white, label=L"N", xlabel=L"\alpha",ylabel=L"D_1")
plot!(alphas, d1', label=ns')

savefig("../figures/nb_d1_alpha.pdf")
plot(alphas, d2[1,:], color = :white, label=L"N", xlabel=L"\alpha",ylabel=L"D_2")
plot!(alphas, d2', label=ns')
savefig("../figures/nb_d2_alpha.pdf")
