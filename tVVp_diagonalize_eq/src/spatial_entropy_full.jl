"""
Calculate both the spatial and the operational entanglement entropies of a
region A, using the SVD. The latter is the "entanglement of particles"
introduced by Wiseman and Vaccaro in 2003.
"""
function spatial_entropy_full(basis::AbstractFermionsbasis, A, d::Vector{ComplexF64})
    B = setdiff(1:basis.K, A)
    A = convert(Array{Int,1},A)
    # Matrices to SVD
    Amatrices = []
    for i=0:basis.N
        DimA = num_vectors(basis, i, length(A))
        DimB = num_vectors(basis, basis.N-i, length(B))

        push!(Amatrices, zeros(ComplexF64, DimA, DimB))
    end

#    norms = zeros(Float64, basis.N+1)

    for i=1: basis.D
        braA = SubKet(basis.vectors[i], A)
        braB = SubKet(basis.vectors[i], B)

        row = serial_num(basis, length(A), count_ones(braA), braA)
        col = serial_num(basis, length(B), count_ones(braB), braB)

        Amatrices[1 + count_ones(braA)][row, col] = d[i]

#        norms[1 + count_ones(braA)] += abs(d[i])^2

    end

#    norm_err = abs(sum(norms) - 1.0)

#    if norm_err > 1e-12
#        warn("norm error ", norm_err)
#    end

    Ss_raw = [svdvals(Amatrix) for Amatrix in Amatrices]

    # Spatial.
    S_sp = vcat(Ss_raw...)
    err_sp = abs(sum(S_sp.^2) - 1.0)
#        warn("RDM eigenvalue error ", S_sp)
    if err_sp > 1e-12
        warn("RDM eigenvalue error ", err_sp)
    end
    S1_sp = 0.0
    for k=1:length(S_sp)
        if abs(S_sp[k])>0
            S1_sp -= abs(S_sp[k])^2*log(abs(S_sp[k])^2)
	end
    end
    S2_sp = -log(abs(sum(S_sp.^4)))

    # Operational.
    #Ss_op = [S / sqrt(n) for (S, n) in zip(Ss_raw, norms)]
    #errs_op = [abs(sum(S.^2) - 1.0) for S in Ss_op]

    #if any(errs_op .> 1e-12)
      #  warn("RDM  eigenvalue error ", maximum(errs_op))
     #   warn("RDM eigenvalue error ", S2_sp)

    #end

    #S2s_op = [-log(sum(S.^4)) for S in Ss_op]
    #S2_op = dot(norms, S2s_op)

    S1_sp, S2_sp
end

spatial_entropy_full(basis::AbstractFermionsbasis, Asize::Int, d::Vector{ComplexF64}) = spatial_entropy_full(basis, 1:Asize, d)
