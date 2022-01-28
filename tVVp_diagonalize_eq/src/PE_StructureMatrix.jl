"""
Compute the structure matrix.
"""
function PE_StructureMatrix(L::Int,N::Int,Cycle_leaders::Vector{Int64}, Asize::Int)
    for x = [:Flips,:SumFlips]
       @eval $x = Int64
    end 
    Bsize = N - Asize
    if Asize>Bsize
        Asize,Bsize=Bsize,Asize
    end 
    if Asize==0
       return zeros(Int64,1,1,1)
    end

    Cycle_leadersA, Cycle_sizesA, NumOfCyclesA = translational_symmetry_cycles_biartition(L,Asize)
    Cycle_leadersB, Cycle_sizesB, NumOfCyclesB = translational_symmetry_cycles_biartition(L,Bsize)

    AmatrixStructure=zeros(Int64,NumOfCyclesA, NumOfCyclesB,L)

    braA = Int
    braB = Int
    bra  = Int

    Aparity= Asize%2
    Bparity= Bsize%2
   # constructs the AmatrixStructure
    OcupationOverlap =2
    for i=1: NumOfCyclesB
        braB=Cycle_leadersB[i]
 
        for j=1: NumOfCyclesA 
            minSize=Cycle_sizesA[j]
            maxSize=Cycle_sizesB[i]
            if Cycle_sizesA[j]>Cycle_sizesB[i]
                minSize,maxSize=maxSize, minSize
            end
            phase0=1
            braA = Cycle_leadersA[j]
            for k=1: minSize           
                # obtain basisA.vectors[CyclesA[j,k]]  
                braA = CircshiftKet(braA,L) 
                 bra=braA+braB
                    # absorbing phase changes due to a particle crossing the system boundary.
                    if Bparity==1 && j>1
                        phase0*= 1-2*CheckSite(braA,1)

                    end
                 if braA & braB==0
                    Flips=0
                    SumFlips=0
                    for Index=1:L
                       Flips = Flips +(1-2* Flips)* CheckSite(braA, Index)
                       SumFlips += Flips* CheckSite(braB,Index)
                    end
                    phase=(-1)^SumFlips*phase0
                    # need InvCycles_Id[serial_num(basis, bra)]
                    AmatrixStructure[j,i,k]= lookup_cylceID_translationq0R1P1Cycles(bra,Cycle_leaders,L) * phase
                end
            end
        end
    end
    return AmatrixStructure

end
