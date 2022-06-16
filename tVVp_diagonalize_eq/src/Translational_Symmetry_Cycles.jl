"""
    translational_symmetry_cycles_biartition(L::Int, subSize::Int)
Create a list of leaders of translational symmetry cycle in the occupation basis 
for each for fermionic 1D chains with PBC/APBC For odd/even N.
"""
function translational_symmetry_cycles_biartition(L::Int, subSize::Int)

    for x = [:NumOfCycles, :Num_cycles_max, :MemberID, :i_next]
        @eval $x = Int64
    end

    basis = Fermionsbasis(L, subSize)

    if subSize == 1
        return Vector{Int64}([basis.vectors[1]]), Vector{Int64}([L]), 1
    end

    IdStatus = Bool
    Num_cycles_max = round(Int64, basis.D / basis.K * 1.5)
    Cycle_leaders = zeros(Int64, Num_cycles_max)
    Cycle_sizes = zeros(Int64, Num_cycles_max)
    Status = trues(basis.D)
    NumOfCycles = 0
    for i = 1:basis.D
        if Status[i]
            NumOfCycles += 1
            Cycle_leaders[NumOfCycles] = basis.vectors[i]
            MemberID = 0
            IdStatus = true
            i_next = i
            while IdStatus
                Status[i_next] = false
                MemberID += 1
                i_next = serial_num(basis, CircshiftKet(basis.vectors[i_next], basis.K))
                IdStatus = Status[i_next]
            end
            Cycle_sizes[NumOfCycles] = MemberID
        end
    end
    return Cycle_leaders[1:NumOfCycles], Cycle_sizes[1:NumOfCycles], NumOfCycles
end

"""
    translational_symmetry_cycles_biartition_noleaders(L::Int, subSize::Int)
Same as translational_symmetry_cycles_biartition but does not store and return the 
leaders.
"""
function translational_symmetry_cycles_biartition_noleaders(L::Int, subSize::Int)

    for x = [:NumOfCycles, :Num_cycles_max, :MemberID, :i_next]
        @eval $x = Int64
    end

    basis = Fermionsbasis(L, subSize)

    if subSize == 1
        return Vector{Int64}([L]), 1
    end

    IdStatus = Bool
    Num_cycles_max = round(Int64, basis.D / basis.K * 1.5)
    Cycle_sizes = zeros(Int64, Num_cycles_max)
    Status = trues(basis.D)
    NumOfCycles = 0
    for i = 1:basis.D
        if Status[i]
            NumOfCycles += 1
            MemberID = 0
            IdStatus = true
            i_next = i
            while IdStatus
                Status[i_next] = false
                MemberID += 1
                i_next = serial_num(basis, CircshiftKet(basis.vectors[i_next], basis.K))
                IdStatus = Status[i_next]
            end
            Cycle_sizes[NumOfCycles] = MemberID
        end
    end
    return Cycle_sizes[1:NumOfCycles], NumOfCycles
end
