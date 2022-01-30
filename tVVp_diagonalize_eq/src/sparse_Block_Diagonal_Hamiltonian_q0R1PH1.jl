using Printf

struct FileHeader_Ham
    M::Int64
    N::Int64 
    len_rows:: Int64 
    len_cols::Int64 
    len_elements::Int64
    t::Float64 
end

"""
    sparse_Block_Diagonal_Hamiltonian_q0R1PH1(L::Int, N::Int ,Cycle_leaders:: Vector{Int64}, Cycle_sizes:: Vector{Int64}, NumOfCycles::Int64, t::Float64, V::Float64, Vp::Float64, load_offdiagonal::Bool=false, save_offdiagonal::Bool=false, storage_path::String="./")
Create the sparse translational, reflection and particle-hole symmetry block H_(q=0,R=1,P=1) of the hamiltonian of fermionic 1D chains with PBC/APBC.
"""
function sparse_Block_Diagonal_Hamiltonian_q0R1PH1(L::Int, N::Int ,Cycle_leaders:: Vector{Int64}, Cycle_sizes:: Vector{Int64}, NumOfCycles::Int64, t::Float64, V::Float64, Vp::Float64, load_offdiagonal::Bool=false, storage_path::String="./")
    if L!=2*N
        warn("particle-hole symmetry works only at half-filling,", "  quit")
        quit()
    end 

    if load_offdiagonal  
        storage_path = joinpath(storage_path, @sprintf "ham_offdiag_N%d_M%d.dat" N L)
    end

    #Creating the block H_(q=0,R=1,P=1) of the hamiltonian.
    end_site = L

    if load_offdiagonal
        rows, cols, elements = load_offdiagonal_terms(storage_path,N,L,t)
    else
        rows = Int64[]
        cols = Int64[]
        elements = Float64[]
        # Off-diagonal part
        for CycleId = 1: NumOfCycles  
            bra = Cycle_leaders[CycleId] 
            CycleId_CycleSize=Cycle_sizes[CycleId]
            for j=1:end_site
                j_next = j % L + 1
                # Tunnel right, tunnel left.
                for (site1, site2) in [(j, j_next), (j_next, j)]
                    if CheckSite(bra,site1) == 1
                        ket = copy(bra)
                        if CheckSite(bra,site2) == 0
                            ket =EmptySite(ket,site1)
                            ket =OccupySite(ket,site2) 
                            CycleId1 = lookup_cylceID_translationq0R1P1Cycles(ket,Cycle_leaders,L)  
                            CycleId1_CycleSize=Cycle_sizes[CycleId1]
                            k1=CycleId1
                            factor=-t*sqrt(CycleId_CycleSize/CycleId1_CycleSize)*.5
                            push!(rows, k1)
                            push!(cols, CycleId)
                            push!(elements, factor)
                            push!(rows, CycleId)
                            push!(cols, k1)
                            push!(elements, conj(factor)) 
                        end
                    end
                end
            end
        end
 
    end

    # Diagonal terms
    for CycleId = 1: NumOfCycles 
        Vsum = 0.0
        Vpsum = 0.0
        bra = Cycle_leaders[CycleId]

        for j=1:end_site
            j_next = j % L + 1
            j_next_next = j_next % L + 1
            Vsum += CheckSite(bra,j) * CheckSite(bra,j_next)
            Vpsum += CheckSite(bra,j) * CheckSite(bra, j_next_next)

        end
        push!(rows, CycleId)
        push!(cols, CycleId)
        push!(elements, Vsum*V+ Vpsum*Vp)
    end 

        return sparse(rows, cols, elements, NumOfCycles, NumOfCycles), NumOfCycles 
end

"""
save_offdiag_sparse_Block_Diagonal_Hamiltonian_q0R1PH1(L::Int, N::Int ,Cycle_leaders:: Vector{Int64}, Cycle_sizes:: Vector{Int64}, NumOfCycles::Int64, t::Float64, storage_path::String="./")
Create the sparse translational, reflection and particle-hole symmetry block H_(q=0,R=1,P=1) of the hamiltonian of fermionic 1D chains with PBC/APBC
and only save the offdiagonal elements to disk.
"""
function save_offdiag_sparse_Block_Diagonal_Hamiltonian_q0R1PH1(L::Int, N::Int ,Cycle_leaders:: Vector{Int64}, Cycle_sizes:: Vector{Int64}, NumOfCycles::Int64, t::Float64, storage_path::String="./")
    if L!=2*N
        warn("particle-hole symmetry works only at half-filling,", "  quit")
        quit()
    end 
    storage_path = joinpath(storage_path, @sprintf "ham_offdiag_N%d_M%d.dat" N L)
    #Creating the block H_(q=0,R=1,P=1) of the hamiltonian.
    end_site = L
 
        rows = Int64[]
        cols = Int64[]
        elements = Float64[]
        # Off-diagonal part
        for CycleId = 1: NumOfCycles  
            bra = Cycle_leaders[CycleId] 
            CycleId_CycleSize=Cycle_sizes[CycleId]
            for j=1:end_site
                j_next = j % L + 1
                # Tunnel right, tunnel left.
                for (site1, site2) in [(j, j_next), (j_next, j)]
                    if CheckSite(bra,site1) == 1
                        ket = copy(bra)
                        if CheckSite(bra,site2) == 0
                            ket =EmptySite(ket,site1)
                            ket =OccupySite(ket,site2) 
                            CycleId1 = lookup_cylceID_translationq0R1P1Cycles(ket,Cycle_leaders,L)  
                            CycleId1_CycleSize=Cycle_sizes[CycleId1]
                            k1=CycleId1
                            factor=-t*sqrt(CycleId_CycleSize/CycleId1_CycleSize)*.5
                            push!(rows, k1)
                            push!(cols, CycleId)
                            push!(elements, factor)
                            push!(rows, CycleId)
                            push!(cols, k1)
                            push!(elements, conj(factor)) 
                        end
                    end
                end
            end
        end
 
            save_offdiagonal_terms(storage_path,L,N,t, rows, cols, elements)
    
        return nothing

end

function save_offdiagonal_terms(storage_path::String, M::Int64,N::Int64,t::Float64, rows::Vector{Int64}, cols::Vector{Int64}, elements::Vector{Float64})
    len_rows = length(rows)
    len_cols = length(cols)
    len_elements = length(elements)

    file_header= FileHeader_Ham(M, N, len_rows, len_cols, len_elements, t)
            H_out=open(storage_path, "w")
            write(H_out, file_header.M, file_header.N, file_header.len_rows, file_header.len_cols, file_header.len_elements,file_header.t,rows,cols,elements)
            flush(H_out)
            close(H_out) 
end

function load_offdiagonal_terms(storage_path::String, N::Int64,M::Int64,t::Float64) 
    file_header =zeros(Int64,5) 
    H_in=open(storage_path, "r")
        read!(H_in,file_header) 
        M_f=file_header[1]
        N_f=file_header[2] 
        len_rows = file_header[3]
        len_cols = file_header[4]
        len_elements = file_header[5]
        file_header2 = zeros(Float64,1)
        read!(H_in,file_header2) 
        t_f = file_header2[1]
        if (M_f!=M) || (N_f!=N) 
            println("the file of states is not compatible with the input parameters" )
            println("M=",M," N=",N)
            println("M_f=",M_f," N_f=",N_f)
            exit(1)
        end 
        rows=zeros(Int64, len_rows)  
        read!(H_in, rows)
        cols=zeros(Int64, len_cols)  
        read!(H_in, cols)
        elements=zeros(Float64, len_elements)  
        read!(H_in, elements)
    close(H_in) 
     
    if (abs(t_f- t)> 1.0E-12)
        println("parameter t is different than in save file" )
        println("t=",t)
        println("t_f=",t_f)
        println("adjust for this by multipying elements with t/t_f ..." ) 
        elements .= t/t_f.*elements
    end

    return rows, cols, elements
end
