"""
Calculate the particle entanglement entropy for an eigenstate of the one-site-translation operator with eigenvalue q=0.
"""
function PE_StructureMatrix(basis::AbstractFermionsbasis, Asize::Int, InvCycles_Id::Vector{Int64})
    for x = [:Flips,:SumFlips]
       @eval $x = Int64
    end
    L = basis.K
    N = basis.N
    Bsize = N - Asize
    if Asize>Bsize
        Asize,Bsize=Bsize,Asize
    end 
    if Asize==0
       return zeros(Int64,1,1,1)
    end
    basisA = Fermionsbasis(L, Asize)
    basisB = Fermionsbasis(L, Bsize)

    CyclesA, CycleSizeA, NumOfCyclesA =Translational_Symmetry_Cycles(basisA)
    CyclesB, CycleSizeB, NumOfCyclesB =Translational_Symmetry_Cycles(basisB)
    AmatrixStructure=zeros(Int64,NumOfCyclesA, NumOfCyclesB,L)


    braA = Int
    braB = Int
    bra  = Int


    Aparity= Asize%2
    Bparity= Bsize%2
   # constructs the AmatrixStructure
    OcupationOverlap =2
    for i=1: NumOfCyclesB
        braB=basisB.vectors[CyclesB[i,1]]
 
        for j=1: NumOfCyclesA 
            minSize=CycleSizeA[j]
            maxSize=CycleSizeB[i]
            if CycleSizeA[j]>CycleSizeB[i]
                minSize,maxSize=maxSize, minSize
            end
            phase0=1
            for k=1: minSize             
                braA= basisA.vectors[CyclesA[j,k]]
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
                    AmatrixStructure[j,i,k]= InvCycles_Id[serial_num(basis, bra)]* phase
                end
            end
        end
    end
    return AmatrixStructure

end
