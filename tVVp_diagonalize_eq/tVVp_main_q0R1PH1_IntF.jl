# Entanglement entropy of the t-V-V` Model at half-filling after a quantum quench 

push!(LOAD_PATH, joinpath(dirname(@__FILE__), "src"))

using tVVpDiagonalize
using ArgParse
using IntFermionicbasis
using Arpack
using Printf
using LinearAlgebra

struct FileHeader
    M::Int64
    N::Int64 
    basis_num:: Int64 
    V::Float64 
    Vp::Float64 
end

# ------------------------------------------------------------------------------
function parse_commandline()
    s = ArgParseSettings()
    s.autofix_names = true
    @add_arg_table s begin
        "M"
            help = "number of sites"
            arg_type = Int
            required = true
        "N"
            help = "number of particles"
            arg_type = Int
            required = true
        "--out"
            metavar = "FILE"
            help = "path to output file"
        "--obdm"
            help = "output the spatial dependence of the OBDM"
            action = :store_true
        "--g2"
            help = "output the pair correlation function ⟨ρ_iρ_0⟩"
            action = :store_true
        "--spatial"
            help = "output the spatial entanglement entropy for ℓ = M/2"
            action = :store_true
        "--save-states"
            help = "save the time-evolved state to disk"
            action = :store_true
        "--load-states"
            help = "load the time-evolved state from disk"
            action = :store_true 
        "--states-file"
            metavar = "FILE"
            help = "path to I/O file for states"

    end
    add_arg_group(s, "boundary conditions")
    @add_arg_table s begin
        "--pbc"
            help = "periodic boundary conditions (default)"
            arg_type = BdryCond
            action = :store_const
            dest_name = "boundary"
            constant = PBC
            default = PBC
        "--obc"
            help = "open boundary conditions"
            arg_type = BdryCond
            action = :store_const
            dest_name = "boundary"
            constant = OBC
            default = PBC
    end
    add_arg_group(s, "tV parameters")
    @add_arg_table s begin 
        "--V"
            metavar = "V"
            help = "final V"
            arg_type = Float64
            default = 1.0 
        "--Vp"
            metavar = "Vp"
            help = "final Vp"
            arg_type = Float64
            default = 0.0
        "--t"
            metavar = "t"
            help = "t value"
            arg_type = Float64
            default = 1.0
    end
    add_arg_group(s, "entanglement entropy")
    @add_arg_table s begin
        "--ee"
            metavar = "ℓ"
            help = "compute all EEs with partition size ℓ"
            arg_type = Int
            required = true
    end

        return parse_args(s, as_symbols=true)
end

function main()
    # _______________Parameter_Setup________________
    c=parse_commandline()
    # Number of sites
    M = c[:M]
    # Number of particles
    N = c[:N]
    if M!=2N
        println("Not at half-filling: the number of sites =", M," and the number of particles =",N )
        exit(1)
    end
    if c[:save_states] && c[:load_states]
        println("only one option can be true.")
        exit(1)
    end

    # Boundary conditions
    boundary = c[:boundary]
    # Size of region A
    Asize = c[:ee]

    # Interaction paramers V and V'
    V = c[:V]
    Vp = c[:Vp]

    # _______________Output_Setup________________
    # Output file
    if c[:out] === nothing
        output = @sprintf "partEE_%02d_%02d_%+5.3f_%+5.3f_%1d.dat" M N V Vp Asize
    else
        output = c[:out]
    end

    # output file if we are measuring the spatial entanglement entropy
    if c[:spatial]
        spat_output = @sprintf "spatEE_%02d_%02d_%+5.3f_%+5.3f_%1d.dat" M N V Vp Asize
    end

    # output file if we are measuring the pair correlation function 
    if c[:g2]
        g2_output = @sprintf "g2_%02d_%02d_%+5.3f_%+5.3f_%1d.dat" M N V Vp Asize
    end

    # state output file
    if c[:save_states] || c[:load_states]
        if c[:states_file] === nothing
            Ψt_output = @sprintf "psioft_%02d_%02d_%+5.3f_%+5.3f.dat" M N V Vp
        else
            Ψt_output=c[:states_file]
        end
    end

    # open and prepare files for output
    if ~c[:save_states]
        f_part = open(output, "w")
        write(f_part, "# M=$(M), N=$(N), V=$(V), Vp=$(Vp), $(boundary)\n")
        write(f_part,@sprintf "#%24s%24s\n" "S₁(n=$(Asize))" "S₂(n=$(Asize))")
        if c[:spatial]
            ℓsize = div(M, 2)
            f_spat = open(spat_output, "w")
            write(f_spat, "# M=$(M), N=$(N), V=$(V), Vp=$(Vp), $(boundary)\n")
            write(f_spat,@sprintf "#%24s%24s\n" "S₁(ℓ=$(ℓsize))" "S₂(ℓ=$(ℓsize))")
        end
        if c[:g2]
            f_g2 = open(g2_output, "w")
            write(f_g2, "# M=$(M), N=$(N), V=$(V), Vp=$(Vp), $(boundary)\n") 
            for x=0:M-1
                write(f_g2,@sprintf "%24d" x )
            end
            write(f_g2,"\n")
        end

    end 

    #_______________________________________________________________________________
    basis = Fermionsbasis(M, N)
    ll=length(basis)
    μ=zeros(Float64, M)
    exp_q=zeros(ComplexF64, basis.K)
    #  Initial wave function in terms of the spatial basis
    Ψn=ComplexF64
    # the one body density matrix
    obdm=zeros(Float64,M,1)
    # Exploit symmetries of the hamiltonian to perform a bloack diagonalization
    Cycles, CycleSize, NumOfCycles, InvCycles_Id, InvCycles_order =Symmetry_Cycles_q0R1PH1(basis)
        #println("size(Cycles) = ",Base.summarysize(Cycles)/1024^3," gb")
    #_______________________________________________________________________________


    if ~c[:load_states]  
        #---------------------------------------find the states---------------------------------------
        # Create the Hamiltonian 
        #H = sparse_hamiltonian(basis,c[:t],V0,Vp0 ,boundary=boundary) 
        H, HRank = sparse_Block_Diagonal_Hamiltonian_q0R1PH1(basis, Cycles, CycleSize, NumOfCycles, InvCycles_Id, InvCycles_order, c[:t],V,Vp) 
        #println("size(H) = ",Base.summarysize(H)/1024^3," gb")
        print(" sparse_hamiltonian finish\n ")
        # H0 = full_hamiltonian(basis, c[:t], V0,boundary=boundary)
        # EigenValues, EigenVectors = eig(H0)
        # wft0 = EigenVectors[:,1]

        #Perform the Lanczos diagonalization to obtain the lowest eigenvector
        # http://docs.julialang.org/en/release-0.3/stdlib/linalg/?highlight=lanczos
        Ψ=zeros(ComplexF64, HRank)
        # I don't understand why this copying is necessary, it is a type conversion thing
        #
        #evals = eigs(H, nev=1, which=:SR,tol=1e-13,v0=getΨ0_trial(c[:t],V0,boundary,basis, HRank, CycleSize, InvCycles_Id))[1]
        #println(evals)
        #exit(0)

        Ψ = eigs(H, nev=1, which=:SR,tol=1e-13,v0=getΨ0_trial(c[:t],V,boundary,basis, HRank, CycleSize, InvCycles_Id))[2][1: HRank].*ones(ComplexF64, HRank)
        #Ψ = eigs(H, nev=1, which=:SR,tol=1e-13,v0=getΨ0_trial(c[:t],V0,boundary,basis, HRank,))[2][1: HRank].*ones(ComplexF64, HRank)
        #println("size(complex) = ", Base.summarysize(Ψ[1]))
        H= Nothing
        Ψ.= Ψ./sqrt(dot(Ψ,Ψ))
 
    else  
        #---------------------------------------load the states---------------------------------------
        Cycles= Nothing
        InvCycles_order= Nothing
        ######file_header1=Array{Int64}(4)
        ######file_header2=Array{Float64}(4)
        file_header1 =zeros(Int64,3)
        file_header2 =zeros(Float64,2)
        Ψf=open(Ψt_output, "r")
            read!(Ψf,file_header1)
            read!(Ψf,file_header2)
            M_f=file_header1[1]
            N_f=file_header1[2] 
            basis_num_f=file_header1[3] 
            V_f=file_header2[1] 
            Vp_f=file_header2[2]   
            if (M_f!=M) || (N_f!=N)  ||(abs(V_f- V)> 1.0E-12)||(abs(Vp_f- Vp)> 1.0E-12)  
                println("the file of states is not compatible with the input parameters" )
                println("M=",M," N=",N," V=",V," Vp=",Vp)
                println("M_f=",M_f," N_f=",N_f," V_f=",V_f," Vp_f=",Vp_f)
                exit(1)
            end 
            Ψ=zeros(ComplexF64, NumOfCycles)  
            read!(Ψf, Ψ)
        close(Ψf)
        HRank = basis_num_f
        HqRank = basis_num_f
    end
    if ~c[:save_states] 
        #---------------------------------------calculate the entanglement---------------------------------------
        AmatrixStructure =PE_StructureMatrix(basis, Asize, InvCycles_Id)

        for j=1: HRank
            Ψ[j]=Ψ[j]/sqrt(CycleSize[j])
        end

        if c[:spatial]
            s_spatial = spatial_entropy(basis, ℓsize, Ψ, InvCycles_Id)
            write(f_spat, @sprintf "%24.12E%24.12E\n" s_spatial[1] s_spatial[2])
            flush(f_spat)
        end

        # measure the pair correlation function
        if c[:g2]
            g2 = pair_correlation(basis,Ψ, InvCycles_Id) 
            for x=1:M
                write(f_g2, @sprintf "%24.12E" g2[x])
            end
            write(f_g2,"\n")
            flush(f_g2)
        end

        if c[:obdm] && Asize == 1
            s_particle,obdm[:] = particle_entropy_Ts(basis, Asize, Ψ,c[:obdm], AmatrixStructure)

        else
            s_particle = particle_entropy_Ts(basis, Asize, Ψ,c[:obdm], AmatrixStructure)

        end

        write(f_part, @sprintf "%24.12E%24.12E\n" s_particle[1] s_particle[2]) 
        
        
        # close output files
        close(f_part)
        if c[:spatial]
        close(f_spat)
        end
        if c[:g2]
            close(f_g2)
        end

        #_______________________________________________________________________________
        # output the time dependent OBDM to disk
        if c[:obdm] && Asize == 1
            obdm_name = @sprintf "obdm_%02d_%02d_%+5.3f_%+5.3f.dat" M N V Vp
                obdm_f = open(obdm_name, "w")
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
    else 
        #---------------------------------------save the states---------------------------------------

        file_header= FileHeader(M, N, HRank, V, Vp)
        Ψf=open(Ψt_output, "w")
            write(Ψf, file_header.M, file_header.N, file_header.basis_num, file_header.V, file_header.Vp,Ψ)
            flush(Ψf)
        close(Ψf)
    end

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

