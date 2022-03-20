module MSCSimPlot

# Simulate a WF population along a species tree for an MSC model
using NewickTree, Plots, StatsBase, Distributions

# Get a layout 
function msc_layout(tree, Ns, pad=5)
    xleaf = 0
    coord = Dict()
    function walk(n)
        if isleaf(n)
            N = Ns[id(n)]
            t = Int(distance(n))
            coord[id(n)] = (xleaf:(xleaf+N), 0:t-1) 
            xleaf += N+pad
        else
            x1, y1 = walk(n[1])
            x2, y2 = walk(n[2])
            N = Ns[id(n)]
            t = Int(distance(n))
            xmid = (x2[end] + x1[1]) ÷ 2
            x = xmid - (N÷2)
            y = y1[end] + 1
            coord[id(n)] = (x:(x+N), y:(y+t-1))
        end
        return coord[id(n)]
    end
    walk(tree)
    return coord
end

function get_ancestry(tree, coords)
    allpaths  = Vector{Tuple{Int,Int}}[]
    fullpaths = Dict()
    function walk(n, t, currpaths)
        xs, ys = coords[id(n)]
        if isroot(n) && t == 1
            y = ys[end]
            currpaths = [[(x,y)] for x=xs]
            walk(n, t+1, currpaths)
        elseif t == distance(n) + 1
            if !isleaf(n)
                walk(n[1], 1, copy(currpaths))
                walk(n[2], 1, copy(currpaths))
            else
                push!(allpaths, currpaths...)
            end
        else
            noff = rand(Multinomial(length(xs), length(currpaths)))
            xoff = xs[1]
            newpths = []
            for (pth, m) in zip(currpaths, noff)
                xp, yp = pth[end]
                if m == 0 
                    push!(allpaths, pth)
                else
                    # prolong the path
                    push!(allpaths, pth)
                    for j=1:m
                        xy = (xoff, yp-1)
                        push!(newpths, [pth[end] ; xy])
                        fullpaths[xy] = [copy(pth) ; (xoff, yp-1)]
                        xoff += 1
                    end
                end
            end
            walk(n, t+1, newpths)
        end
    end
    walk(tree, 1, Vector{Tuple{Int,Int}}[])
    return allpaths, fullpaths
end

function traceback(paths, x)
    pth = [x]
    while haskey(paths, x)
        push!(pth, reverse(paths[x])...)
        x = paths[x][1]
    end
    return pth
end

end
