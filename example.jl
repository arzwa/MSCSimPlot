using MSCSimPlot, Plots, NewickTree, StatsBase
theme(:default)

S = nw"(((A:20,B:20):10,C:30):5,(D:25,E:25):10);"
N = Dict(id(n)=>rand(10:30) for n in postwalk(S))
S.data.distance=20  # root branch

d = MSCSimPlot.msc_layout(S, N)
P, X = MSCSimPlot.get_ancestry(S, d)


p = plot(legend=false, size=(1200,900), gridstyle=:dot, xshowaxis=false,
         xticks=false, ylabel="generations in the past", ylim=(-0.5,Inf))
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


# make an inset of the species tree and gene tree
plot!(p,
    inset = (1, bbox(0.01, 0.03, 0.35, 0.25, :top, :left)),
    subplot = 1,
)
p2 = p[2]
MSCSimPlot.treeshape!(p2, S, d, color="#ECECEC", linealpha=0., legend=false, framestyle=:none)
pths = copy.(first.(paths))
for t=1:length(pths[1])
    coal = []
    for i=2:length(pths)
        for j=1:i-1
            if pths[i][t] == pths[j][t]
                push!(coal, (i, j, t))
            end
        end
    end
    todel = Int[]
    for (i,j,t) in coal
        plot!(p2, [pths[i][1], pths[i][t]], color=:black)
        plot!(p2, [pths[j][1], pths[j][t]], color=:black)
        pths[i][1] = pths[i][t]
        push!(todel, j)
    end
    deleteat!(pths, todel)
end

plot(p, margin=5mm)


