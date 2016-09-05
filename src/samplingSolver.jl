
#a massively stripped-down version of the samplingSolver file in the original repository. The complete file is at samplingSolver-old

function grid3(n::Int64)
    gr2 = grid2(n);
    a = kron(speye(n), gr2);
    b = kron(gr2, speye(n));
    gr3 = sparse(a + b);
    gr3.nzval = Float64[1 for i in 1:nnz(gr3)]
    gr3 = tril(gr3) + tril(gr3)'
    return gr3
end

function samplingLDL{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti})

    n = a.n
    eps = 0.1
    rho::Tv = ceil(0.2 * log(n) ^ 2 / eps ^ 2)
    ord = reverse!(dfsOrder(akpw(a), start = n));
    a = a[ord,ord]

    auxVal = zeros(Tv, n)                       # used to sum weights from multiedges
    auxMult = zeros(Tv, n)                      # used to count the number of multiedges

    wNeigh = zeros(Tv, n)
    multNeigh = zeros(Tv, n)
    indNeigh = zeros(Ti, n)
    asdfg = 0
    u = spzeros(n,n) 
    d = zeros(Tv, n)                           

    # neigh[i] = the list of neighbors for vertex i with their corresponding weights
    # note neigh[i] only stores neighbors j such that j > i
    # neigh[i][1] is weight, [2] is number of multi-edges, [3] is neighboring vertex

    neigh = llsInit(a)

    # gather the info in a and put it into neigh and w
    for i in 1:length(a.colptr) - 1

        for j in a.colptr[i]:a.colptr[i + 1] - 1
            if a.rowval[j] > i
                llsAdd(neigh, i, (a.nzval[j], rho, a.rowval[j]))
            end
        end
    end

    # Now, for every i, we will compute the i'th column in U
    for i in 1:(n-1)

        # We will get rid of duplicate edges
        # wSum - sum of weights of edges
        # multSum - sum of number of edges (including multiedges)
        # numPurged - the size in use of wNeigh, multNeigh and indNeigh
        # wNeigh - list of weights correspongind to each neighbors
        # multNeigh - list of number of multiedges to each neighbor
        # indNeigh - the indices of the neighboring vertices
        wSum, multSum, numPurged = llsPurge(neigh, i, auxVal, auxMult, wNeigh, multNeigh, indNeigh, rho = rho)

        # need to divide weights by the diagonal entry
        for j in 1:numPurged
	    u[indNeigh[j],i] = -wNeigh[j] / wSum
        end
	u[i,i] = 1	
        d[i] = wSum

        multSum = ceil(Int64, multSum)
        wSamp = FastSampler(wNeigh[1:numPurged])
        multSamp = FastSampler(multNeigh[1:numPurged])
	if multSum > asdfg
		println(multSum)
		asdfg = multSum
	end        
        jSamples = sampleMany(wSamp, multSum)
        kSamples = sampleMany(multSamp, multSum)
        
        # now propagate the clique to the neighbors of i
        for l in 1:multSum
            
            j = jSamples[l]
            k = kSamples[l]

            if j != k
                posj = indNeigh[j]
                posk = indNeigh[k]

                # swap so posj is smaller
                if posk < posj  
                    j, k = k, j
                    posj, posk = posk, posj
                end

                wj = wNeigh[j]                
                wk = wNeigh[k]

                sampScaling = wj * multNeigh[k] + wk * multNeigh[j]
                
                llsAdd(neigh, posj, (wj * wk / sampScaling, 1.0, posk))
            end
        end  
    end

    u[n,n] = 1
    d[n] = 0
    return u, d, ord
end

