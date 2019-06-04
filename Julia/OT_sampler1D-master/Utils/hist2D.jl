# Draws histogram for a given one dimesional data
function hist2D(data,binsizex)
    xmax = maximum(data[:,1])
    xmin = minimum(data[:,1])
    xmax = xmax+0.15*(xmax-xmin)
    xmin = xmin-0.15*(xmax-xmin)

    # X grid point locations
    xpoints = linspace(xmin, xmax, binsizex+1)
    # Bins
    bins = zeros(binsizex,1)
    # Bin allocation
    for i = 1:size(data)[1]
        tmp1 = data[i,1]*ones(binsizex,1) - xpoints[1:(end-1),1]
        [tmp1[i]=(tmp1[i] >= 0 ? 1 : 0) for i in 1:binsizex]
        xbinindx = Int(sum(tmp1))
        # print(tmp1)
        # if xbinindx==0
        #     error("Zero in column 1 not possible")
        # end
        bins[xbinindx,1] = bins[xbinindx,1] + 1
    end
    return bins , xpoints
end
