"""Measure the pair correlation function (density-density) exploiting translational symmetry """
function pair_correlation(L::Int, N::Int, d::Vector{ComplexF64}, Cycle_leaders::Vector{Int64})

    # setup the needed vectors
    g2 = zeros(Float64, L)

    # number of basis states  
    D = num_vectors(N, L)
    # measure the pair correlation function (density-density) exploiting
    # translational symmetry
    for i = 1:L
        for v in BasisKetGenerator(L, N)#b=1:D
            #config = v      #basis.vectors[b]
            weight = d[lookup_cylceID_translationq0R1P1Cycles(v, Cycle_leaders, L)]       #d[InvCycles_Id[b]]
            n0 = CheckSite(v, 1) #CheckSite(config,1)
            ni = CheckSite(v, i) #CheckSite(config,i)
            g2[i] += n0 * ni * abs(weight)^2
        end
    end

    return g2
end
