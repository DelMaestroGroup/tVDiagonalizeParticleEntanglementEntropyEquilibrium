"""Measure the pair correlation function (density-density) exploiting translational symmetry """
function pair_correlation(basis::AbstractFermionsbasis, d::Vector{ComplexF64}, InvCycles_Id::Vector{Int64})
    
        # setup the needed vectors
        g2 = zeros(Float64, basis.K)
    
        # measure the pair correlation function (density-density) exploiting
        # translational symmetry
        for i=1:basis.K
            for b=1:basis.D
                config = basis.vectors[b]
                weight = d[InvCycles_Id[b]]
                n0 = CheckSite(config,1)
                ni = CheckSite(config,i)
                g2[i] += n0*ni*abs(weight)^2
            end
        end
    
        return g2
    end