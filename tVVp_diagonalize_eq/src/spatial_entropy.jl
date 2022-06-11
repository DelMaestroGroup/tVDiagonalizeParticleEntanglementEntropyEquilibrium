"""
Calculate both the spatial and the operational entanglement entropies of a
region A, using the SVD. The latter is the "entanglement of particles"
introduced by Wiseman and Vaccaro in 2003.
"""
function spatial_entropy(L::Int, N::Int, A, d::Vector{ComplexF64},Cycle_leaders::Vector{Int64})
    B = setdiff(1:L, A)
    A = convert(Array{Int,1},A)
    # Matrices to SVD
    Amatrices = []
    for i=0:N
        DimA = num_vectors(i, length(A))
        DimB = num_vectors(N-i, length(B))

        push!(Amatrices, zeros(ComplexF64, DimA, DimB))
    end

#    norms = zeros(Float64, basis.N+1)

    for v in BasisKetGenerator(L,N)
        braA = SubKet(v, A)
        braB = SubKet(v, B)

        row = serial_num(length(A), count_ones(braA), braA)
        col = serial_num(length(B), count_ones(braB), braB)

        #element is d[InvCycles_Id[i]] 
        Amatrices[1 + count_ones(braA)][row, col] = d[lookup_cylceID_translationq0R1P1Cycles(v,Cycle_leaders,L)]

#        norms[1 + count_ones(braA)] += abs(d[i])^2

    end

#    norm_err = abs(sum(norms) - 1.0)

#    if norm_err > 1e-12
#        warn("norm error ", norm_err)
#    end

    Ss_raw = [svdvals(Amatrix) for Amatrix in Amatrices]


    # Spatial.
    S_sp = vcat(Ss_raw...)
    λ = abs.(S_sp.^2)
    # sort lambda for numerical stability
    sort!(λ)
    err_sp = abs(sum(λ) - 1.0) 
    if err_sp > 1e-12
        warn("RDM eigenvalue error ", err_sp)
    end
    Sα = zeros(Float64,12)  
    for k=1:length(λ)
        if λ[k] >0
            Sα[1] -= λ[k] *log(λ[k])
	end
    end
    for α = 2:10
        Sα[α] = -log(abs(sum(λ.^(α))))/(α-1)
    end
    Sα[11] = 2.0*log(abs(sum(S_sp)))
    # Renyi for alpha->inf
    Sα[12] = -1.0*log(maximum(λ))

    # Operational.
    #Ss_op = [S / sqrt(n) for (S, n) in zip(Ss_raw, norms)]
    #errs_op = [abs(sum(S.^2) - 1.0) for S in Ss_op]

    #if any(errs_op .> 1e-12)
      #  warn("RDM  eigenvalue error ", maximum(errs_op))
     #   warn("RDM eigenvalue error ", S2_sp)

    #end

    #S2s_op = [-log(sum(S.^4)) for S in Ss_op]
    #S2_op = dot(norms, S2s_op)

    Sα
end

spatial_entropy(L::Int, N::Int, Asize::Int, d::Vector{ComplexF64},Cycle_leaders::Vector{Int64}) = spatial_entropy(L,N, 1:Asize, d, Cycle_leaders)

 