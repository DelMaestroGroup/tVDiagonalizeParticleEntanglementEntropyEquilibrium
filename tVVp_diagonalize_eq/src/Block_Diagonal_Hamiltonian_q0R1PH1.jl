"""
Create the translational, reflection and particle-hole symmetry block H_(q=0,R=1,P=1) of the hamiltonian of fermionic 1D chains with PBC/APBC.
"""
function Block_Diagonal_Hamiltonian_q0R1PH1(basis::AbstractFermionsbasis, Cycles:: Array{Int64,2}, CycleSize:: Vector{Int64}, NumOfCycles::Int64, InvCycles_Id:: Vector{Int64}, InvCycles_order:: Vector{Int64}, t::Float64, V::Float64, Vp::Float64)
    if basis.K!=2*basis.N
        warn("particle-hole symmetry works only at half-filling,", "  quit")
        quit()
    end
    N=basis.N

    #Creating the block H_(q=0,R=1,P=1) of the hamiltonian.
    H011=zeros(Float64, NumOfCycles, NumOfCycles)
    end_site = basis.K
    for CycleId =1: NumOfCycles
        # Diagonal part
        Vsum = 0.0
        Vpsum = 0.0
        bra=basis.vectors[Cycles[CycleId,1]]
        for j=1:end_site
            j_next = j % basis.K + 1
            j_next_next = j_next % basis.K + 1
            Vsum += CheckSite(bra,j) * CheckSite(bra,j_next)
            Vpsum += CheckSite(bra,j) * CheckSite(bra, j_next_next)
        end
        H011[CycleId, CycleId]= Vsum*V+Vpsum*Vp
        # Off-diagonal part

        CycleId_CycleSize=CycleSize[CycleId]
        for j=1:end_site
             j_next = j % basis.K + 1
             # Tunnel right, tunnel left.
             for (site1, site2) in [(j, j_next), (j_next, j)]
                 if CheckSite(bra,site1) == 1
                     ket = copy(bra)
                      if CheckSite(bra,site2) == 0
                         ket =EmptySite(ket,site1)
                         ket =OccupySite(ket,site2)
                         kId=serial_num(basis, ket)
                         CycleId1 = InvCycles_Id[kId]
                         kIdcy = InvCycles_order[kId]
                         CycleId1_CycleSize=CycleSize[CycleId1]
                         k1=CycleId1
                         factor=-t*sqrt(CycleId_CycleSize/CycleId1_CycleSize)
                         H011[k1, CycleId]+= factor
                     end
                 end
             end
         end
    end
        return H011, NumOfCycles
end