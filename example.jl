using MSCSimPlot, Plots, NewickTree, StatsBase
theme(:default)

S = nw"(((A:20,B:20):10,C:30):5,(D:25,E:25):10);"
N = Dict(id(n)=>rand(10:30) for n in postwalk(S))
S.data.distance=20  # root branch

d = MSCSimPlot.msc_layout(S, N)
P, X = MSCSimPlot.get_ancestry(S, d)

p = plot(legend=false, size=(1200,900), grid=false, framestyle=:none)
for x in P
    plot!(x, marker=4, color=:lightgray, markerstrokecolor=:lightgray)
end
paths = []
for (i,l) in enumerate(getleaves(S))
    xs, ys = d[id(l)]
    x0s = sample(xs, 3, replace=false)
    for x0 in x0s
        pth = MSCSimPlot.traceback(X, (x0,0))
        plot!(pth, color=i, lw=2, alpha=0.5)
        push!(paths, (pth, i))
    end
end
plot(p)

for (i,l) in enumerate(getleaves(S))
    xs, ys = d[id(l)]
    x0s = sample(xs, 3, replace=false)
    for x0 in x0s
        pushMSCSimPlot.traceback(X, (x0,0)), color=i, lw=2, alpha=0.5)
    end
end
