# Entanglement entropy of the t-V-V` Model at half-filling after a quantum quench 

push!(LOAD_PATH, joinpath(dirname(@__FILE__), "src"))

using tVVpDiagonalize
using ArgParse
using IntFermionicbasis 
using Printf 

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
            Ψ_output = @sprintf "psioft_%02d_%02d_%+5.3f_%+5.3f.dat" M N V Vp
        else
            Ψ_output=c[:states_file]
        end
    end

    # open and prepare files for output
    if ~c[:save_states]
        f_part = open(output, "w")
        write(f_part, "# M=$(M), N=$(N), V=$(V), Vp=$(Vp), $(boundary)\n")
        write(f_part,@sprintf "#%24s#%24s%24s#%24s#%24s#%24s#%24s#%24s#%24s#%24s#%24s\n"  "S₀₋₅(n=$(Asize))" "S₁(n=$(Asize))" "S₂(n=$(Asize))" "S₃(n=$(Asize))" "S₄(n=$(Asize))" "S₅(n=$(Asize))" "S₆(n=$(Asize))" "S₇(n=$(Asize))" "S₈(n=$(Asize))" "S₉(n=$(Asize))" "S₁₀(n=$(Asize))")
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
    # setting up the basis
    basis = Fermionsbasis(M, N) 
    # the one body density matrix
    obdm=zeros(Float64,M,1)
    # Exploit symmetries of the hamiltonian to perform a bloack diagonalization
    Cycles, CycleSize, NumOfCycles, InvCycles_Id, InvCycles_order =Symmetry_Cycles_q0R1PH1(basis) 


    if ~c[:load_states]  
        #---------------------------------------find the ground state using Lanczos diagonalization---------------------------------------
        Ψ, HRank = ground_state(basis,Cycles, CycleSize, NumOfCycles, InvCycles_Id, InvCycles_order, c[:t],V,Vp,boundary)  
    else  
        Cycles= Nothing
        InvCycles_order= Nothing
        #---------------------------------------load the ground state---------------------------------------
        Ψ, HRank = load_ground_state(Ψ_output,M,N,V,Vp,NumOfCycles)
    end

    if ~c[:save_states] 
        #---------------------------------------calculate the entanglement---------------------------------------
        AmatrixStructure =PE_StructureMatrix(basis, Asize, InvCycles_Id)

        Ψ_coeff = Ψ
        for j=1: HRank
            Ψ_coeff[j]=Ψ_coeff[j]/sqrt(CycleSize[j])
        end

        if c[:spatial]
            s_spatial = spatial_entropy(basis, ℓsize, Ψ_coeff, InvCycles_Id)
            write(f_spat, @sprintf "%24.12E%24.12E\n" s_spatial[1] s_spatial[2])
            flush(f_spat)
        end

        # measure the pair correlation function
        if c[:g2]
            g2 = pair_correlation(basis,Ψ_coeff, InvCycles_Id) 
            for x=1:M
                write(f_g2, @sprintf "%24.12E" g2[x])
            end
            write(f_g2,"\n")
            flush(f_g2)
        end

        if c[:obdm] && Asize == 1
            s_particle,obdm[:] = particle_entropy_Ts(basis, Asize, Ψ_coeff,c[:obdm], AmatrixStructure)

        else
            s_particle = particle_entropy_Ts(basis, Asize, Ψ_coeff,c[:obdm], AmatrixStructure)

        end

        write(f_part, @sprintf "%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E\n" s_particle[11] s_particle[1] s_particle[2] s_particle[3] s_particle[4] s_particle[5] s_particle[6] s_particle[7] s_particle[8] s_particle[9] s_particle[10]) 
        
        
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
            save_obdm(obdm,M,N,V,Vp)
        end
    else 
        #---------------------------------------save the ground state---------------------------------------
        save_ground_state(Ψ_output,Ψ,M,N,HRank,V,Vp)
    end

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

