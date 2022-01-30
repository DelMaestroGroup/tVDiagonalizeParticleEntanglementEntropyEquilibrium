using Printf

"""
Calculate structure matrix and entropy in one function to ensure that memory from structure matrix is freed.
"""
function particle_entropy_Ts_and_structureMatrix(L::Int, N::Int, Asize::Int, d::Vector{ComplexF64}, Cycle_leaders::Vector{Int64}, measure_obdm::Bool)
    AmatrixStructure = PE_StructureMatrix(L, N, Cycle_leaders, Asize) 
    return particle_entropy_Ts(L,N, Asize, d, measure_obdm, AmatrixStructure)
end


"""
Calculate the particle entanglement entropy for an eigenstate of the one-site-translation operator with eigenvalue q=0.
"""
function particle_entropy_Ts(L::Int, N::Int, Asize::Int, d::Vector{ComplexF64}, measure_obdm::Bool, AmatrixStructure:: Array{Int64,3})

    SRen = zeros(Float64,11)
    DimA=Int64
    DimB=Int64
    DimAdA=Int64
    for x = [:fN,:Wa,:Wb]
       @eval $x = Int64
    end
    norm = Float64
    facto = Float64
    AdA_elem= ComplexF64 

    Bsize = N - Asize
    if Asize>Bsize
        Asize,Bsize=Bsize,Asize
    end
    if Asize==0
       return  SRen
    end
    # Dimensions of partition Hilbert spaces
    facto=1.0
    for i=1:Asize
    	facto *=((1.0*L-i+1)/(1.0*i))
    end
    DimA = Int(round(facto))
    facto=1.0
    for i=1:Bsize
    	facto *=((1.0*L-i+1)/(1.0*i))
    end
    DimB = Int(round(facto))

    DimAdA= DimA 
    # Get cycle sizes and number of cylces in both subspaces but skip leader caclulation/stroage
    CycleSizeA, NumOfCyclesA = translational_symmetry_cycles_biartition_noleaders(L, Asize) 
    CycleSizeB, NumOfCyclesB = translational_symmetry_cycles_biartition_noleaders(L, Bsize) 
 
    λ=zeros(Float64, NumOfCyclesA*L)

    # Weight factors
    Wa=factorial(Asize)
    Wb=factorial(Bsize)

    # Normalization coefficient
    norm=sqrt(Wa*Wb/factorial(N))
    Aparity= Asize%2
    Bparity= Bsize%2
    element= ComplexF64

    #find the spectrum of the reduced density matrix 
    for q=0:L-1
        Amatrixq=zeros(ComplexF64, NumOfCyclesA, NumOfCyclesB)
        for i=1: NumOfCyclesB
            for j=1: NumOfCyclesA 
                minSize=CycleSizeA[j]
                maxSize=CycleSizeB[i]
                minparity= Aparity
                maxparity= Bparity
                minCycleparity= div(Asize* CycleSizeA[j],L)%2 
                maxCycleparity= div(Bsize* CycleSizeB[i],L)%2 
                if CycleSizeA[j]>CycleSizeB[i]
                    minSize,maxSize=maxSize, minSize
                    minparity, maxparity= maxparity ,minparity
                    minCycleparity, maxCycleparity= maxCycleparity, minCycleparity
                end
                shiftmin=div(div(L,minSize)-minparity,2)* maxparity* minCycleparity
                shiftmax=div(div(L,maxSize)-maxparity,2)* minparity* maxCycleparity
                element=0.0+0.0im
               if  (q-shiftmin)%div(L,minSize)==0 && (q-shiftmax)%div(L,maxSize)==0
                    halfe=0.5*maxparity * minCycleparity
                    factor1=norm*sqrt(maxSize)/sqrt(minSize)
                    phasefactor=(0.0+1.0im)*2*pi/minSize
                    qValue=(q-shiftmin)/div(L,minSize)
                    for k=1:minSize
                        InvCycles_Id = AmatrixStructure[j,i,k]
                        if InvCycles_Id != 0 
                            phase= sign(InvCycles_Id)
                            InvCycles_Id =abs(InvCycles_Id)
                            element += factor1*phase*d[InvCycles_Id]*exp((1.0+0.0im)*(k-1)*(qValue +halfe)* phasefactor)
                        end 
                    end
                end
                Amatrixq[j,i]= element
          end
       end
       S = svdvals!(Amatrixq)
       S.=S.^2
       λ[q*NumOfCyclesA+1:q*NumOfCyclesA+NumOfCyclesA]=S[1:NumOfCyclesA]
    end

   # construct the spatial OBDM
   if measure_obdm && Asize == 1
       obdm = zeros(Float64, DimA)
       phase=(0.0+1.0im)*2*pi/L
       shift=(0.5* Bparity)
       j_index=Int(L/2)
       for i=1:L
           AdA_elem=0.0+0.0im
           for k=1:L
               AdA_elem += exp(phase*(i-j_index)*(k-1+shift))*λ[k]/L
           end
           obdm[i]=real(AdA_elem)
       end
   end

    err = abs(sum(λ.^1) - 1.0)
    if err > 1e-12
        print("RDM eigenvalue error ", err)
    end
#  println("  λ ", λ )

    LogNn=log(factorial(N)/factorial(Bsize)/factorial(Asize))
    SRen[1]=0 
    for k=1: NumOfCyclesA*L
        if λ[k]>0
            SRen[1] += λ[k]*log(λ[k]) 
	    end
    end
    SRen[1]=-SRen[1]-LogNn
    # all Renyis from 2 to 10
    for α = 2:10
        SRen[α]=-log(sum(λ.^α))/(α-1)-LogNn
    end 
    # Renyi of 1/2: 
    SRen[11] =2*log(sum(sqrt.(λ)))-LogNn

    if measure_obdm && Asize==1
        return SRen,obdm
    else
        return SRen
    end
end


function save_obdm(obdm::AbstractArray{Float64},M::Int64,N::Int64,V::Float64,Vp::Float64,storage_folder::String="./")
    obdm_name = @sprintf "obdm_%02d_%02d_%+5.3f_%+5.3f.dat" M N V Vp
    path = joinpath(storage_folder,obdm_name)
                obdm_f = open(path, "w")
                write(obdm_f, @sprintf "#%11s" "|i-j|") 
                write(obdm_f, "\n")
                flush(obdm_f)

                for i = 1:M
                    write(obdm_f, @sprintf "%16d" (i-Int(M/2))) 
                    write(obdm_f, @sprintf "%16.6E" obdm[i]) 
                    write(obdm_f, "\n")
                    flush(obdm_f)
                end
                close(obdm_f)
end
