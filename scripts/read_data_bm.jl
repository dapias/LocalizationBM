using JLD, Plots, LaTeXStrings

Ts = collect(0.3:0.1:1.2)
#Ts = collect(0.3:0.1:1.1)

a = load("../data/fix_T=$(Ts[1]).jld")

ns = a["ns"]
ents = a["entropies"]
d1 = a["entropies"]./(log.(ns))
ipr = a["ipr"]
d2 = -log.(ipr)./log.(ns)

try
for k in 2:length(Ts)
    a =  load("../data/fix_T=$(Ts[k]).jld")
    ents = a["entropies"]
    global d1 = hcat(d1, a["entropies"]./(log.(ns)))
    ipr = a["ipr"]
    global d2 = hcat(d2, -log.(ipr)./log.(ns))
end
catch
end

plot(Ts, d1[1,:], color = :white, label=L"N", xlabel=L"T",ylabel=L"D_1")
plot!(Ts, d1', label=ns')

savefig("../figures/bm_d1_T.pdf")
plot(Ts, d2[1,:], color = :white, label=L"N", xlabel=L"T",ylabel=L"D_2")
plot!(Ts, d2', label=ns')
savefig("../figures/bm_d2_T.pdf")
