"""
    symmetry_cycles(K::Int,N::Int)
Construc the fermion basis with N fermions on K sites, constuct the symmetry cycles
and only return the cycle leaders and the vector of cycle sizes.
"""
function symmetry_cycles_q0R1PH1(K::Int, N::Int)

    if K != 2 * N
        warn("particle-hole symmetry works only at half-filling,", "  quit")
        quit()
    end
    # Makes a Symbol for each, of type Int64
    for x = [:NumOfCycles, :Num_cycles_max, :MemberID, :i_next]
        @eval $x = Int64
    end
    # construct the basis
    basis = Fermionsbasis(K, N)
    # length of the basis
    D = basis.D
    # This is a rough upper bound on the number of cycles. It could me smaller.
    Num_cycles_max = round(Int64, basis.D / basis.K * 1.5)
    # Only store cylce leader and sizes of cycles 
    Cycles_leaders = zeros(Int64, Num_cycles_max)
    Cycle_sizes = zeros(Int64, Num_cycles_max)

    # to check whether state was already seen before (=false) or is new (=true)
    IdStatus = Bool
    Status = trues(D)
    # Counter variables (total number of cylces NumOfCycles, length of each cycle MemberID)
    NumOfCycles = 0
    MemberID = 0
    for i = 1:D
        # next basis state
        i_next = i
        # j=4 begin of potential new cycle
        j = 4
        # start cycle search
        while j > 0
            # if element was never seen before
            if Status[i_next]
                # if new cycle
                if j == 4
                    NumOfCycles += 1
                    MemberID = 0
                    # add cycle leader (will be always the smallest int in the cylce, as this is how the basis is sorted)
                    Cycles_leaders[NumOfCycles] = basis.vectors[i_next]
                end
                # we end up here whenever we have a new subcylce to
                # map with translations T
                IdStatus = true
                while IdStatus
                    # mark i_next as seen
                    Status[i_next] = false
                    # count i_next towards cycle members
                    MemberID += 1
                    # shift to get the next element from this cycle
                    i_next = serial_num(basis, CircshiftKet(basis.vectors[i_next], K))
                    # stop the loop as we are back at the beginning
                    # of the cycle, at an element we already saw
                    IdStatus = Status[i_next]
                end
            else
                # if element has already been seen go directly
                # to the next basis element, as the complete cycle
                # is already known
                j = 0
            end
            # we get here only if element was new previously, and we rotated
            # i_next back to the cycle leader
            j -= 1
            if j == 3
                # cylce leader is i_next, we shifted through the cycle with T 
                # try now to reverse the cycle leader
                i_next = serial_num(basis, ReverseKet(basis.vectors[i], K))
                # already seen?
                if ~Status[i_next]
                    # the element is in a known cycle
                    # try particle hole P|i>
                    i_next = serial_num(basis, FlipKet(basis.vectors[i], K))
                    # already seen?
                    if ~Status[i_next]
                        # yes, P and R map to a known cycle, go down
                        # and store the size of this now finished cylce
                        # then use next basis state
                        j = 0
                    else
                        # no, P yields a new cycle but R does not, go on with
                        # the reversed and flipped i_next below
                        j = 1
                    end
                end
            end
            if j == 2
                # we are here if R|l> was in a new cycle
                # after mapping out all tranlations TR|i>
                # try also to flip the basis state P|i>
                i_next = serial_num(basis, FlipKet(basis.vectors[i], K))
                if ~Status[i_next]
                    # P does not map to a new cycle, so we are done
                    # with this cycle (only T and R)
                    j = 0
                end
                # otherwise we need to map this cycle above, so
                # we would continue with j=2           
            elseif j == 1
                # We end up here when eighter only P maps
                # to a new cycle or both P and R do and R 
                # has already been mapped out by T
                i_next = serial_num(basis, ReverseKet(FlipKet(basis.vectors[i], K), K))
                # continue with j=1 and shift with T above
            end
            if j == 0
                # end of cycle
                Cycle_sizes[NumOfCycles] = MemberID
            end
        end
    end
    return Cycles_leaders[1:NumOfCycles], Cycle_sizes[1:NumOfCycles], NumOfCycles
    #|============================================================================
    #| Cycles_leaders:stores the serial numbers of leader in each cycle.         |
    #| Cycle_sizes:stores the number of spatial Kets in each cycle. 
    #| NumOfCycles:stores the total number of cycles len(Cycle_sizes[Cycle_sizes.>0])|
    #| Note that Cycles_leaders is sorted and leaders have the smallest int in   |
    #| the cycle by construction of the basis                                    |
    #|============================================================================
end

"""
Finds the cycle ID of a bra of the basis from which the cycle leaders were constructed.
"""
function lookup_cylceID_translationq0R1P1Cycles(bra::Int64, Cycle_leaders::Vector{Int64}, L::Int)
    ket_leader = bra
    for iop = 1:4
        # respective symmetry operation
        if iop == 1
            candidate = bra
        elseif iop == 2
            candidate = FlipKet(bra, L)
        elseif iop == 3
            candidate = ReverseKet(bra, L)
        else
            candidate = FlipKet(ReverseKet(bra, L), L)
        end
        # ket leader is the minimum in the cycle
        if candidate < ket_leader
            ket_leader = candidate
        end
        # all shifts
        for ishift = 1:(L-1)
            candidate = CircshiftKet(candidate, L)
            if candidate < ket_leader
                ket_leader = candidate
            end
        end
    end
    # return postion of the leader of bra's cycle in cycle leaders vector
    return searchsortedfirst(Cycle_leaders, ket_leader)
end
