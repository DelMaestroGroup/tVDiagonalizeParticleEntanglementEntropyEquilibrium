# Entanglement entropy of the t-V-V` Model at half-filling - dependence on interaction strength V

push!(LOAD_PATH, joinpath(dirname(@__FILE__), "src"))

using tVVpDiagonalize
using ArgParse
using IntFermionicbasis 
using Printf 
using ProgressBars
using Dates

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
            metavar = "FOLDER"
            help = "path to output folder"
        "--tmp"
            metavar = "FOLDER"
            help = "folder for hamiltonian storage location"
        "--out-obdm" 
            metavar = "FOLDER"
            help = "folder for obdm storage location (if --obdm provided)"
        "--g2"
            help = "output the pair correlation function ⟨ρ_iρ_0⟩"
            action = :store_true
        "--obdm"
            help = "output the spatial dependence of the OBDM"
            action = :store_true
        "--spatial"
            help = "output the spatial entanglement entropy for ℓ = M/2"
            action = :store_true
        "--skip_hoffdiag_saving"
            help = "do not save offdiagonal terms of Hamiltonian for V=0" 
            action = :store_true
        "--skip_hoffdiag_loading"
            help = "do not load offdiagonal terms of Hamiltonian from V=0 case" 
            action = :store_true
        "--no-flush"
            help = "do not flush write buffer to output files in after computation for each V" 
            action = :store_true
        "--no-recompute-structure-matrix"
            help = "compute structure matrix once and store it in memory (time efficient but memory intensive)" 
            action = :store_true 

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
        "--V_start"
            metavar = "V_start"
            help = "start V"
            arg_type = Float64
            default = -2.0
        "--V_end"
            metavar = "V_end"
            help = "end V"
            arg_type = Float64
            default = 2.0
        "--V_step"
            metavar = "V_step"
            help = "step in V"
            arg_type = Float64
            default = 0.1 
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

"""Use 'write' to write string to IOstream (e.g. write to a file) and flush IOstream if toflush is true."""
function write_flush(stream::IO,str::String,toflush::Bool=true)
    write(stream, str)
    if toflush
        flush(stream)
    end
end

"""Runs entanglement calculation for a range of interaction strenths V based on the input parameters.

Usage: 
        julia tVVp_Vdependence_q0R1PH1_IntF.jl 
                    [--out FOLDER] [--tmp FOLDER]
                    [--out-obdm FOLDER] [--g2] [--obdm]
                    [--spatial] [--skip-hoffdiag-saving]
                    [--skip-hoffdiag-saving]
                    [--no-flush] [--pbc] [--obc]
                    [--V-start V_start] [--V-end V_end]
                    [--V-step V_step] [--Vp Vp] [--t t] --ee ℓ M N

positional arguments:
                    M                     number of sites (type: Int64)
                    N                     number of particles (type: Int64)
                  
optional arguments:
                    --out FOLDER          path to output folder
                    --tmp FOLDER          folder for hamiltonian storage location
                    --out-obdm FOLDER     folder for obdm storage location (if --obdm
                                          provided)
                    --g2                  output the pair correlation function ⟨ρ_iρ_0⟩
                    --obdm                output the spatial dependence of the OBDM
                    --spatial             output the spatial entanglement entropy for ℓ
                                          = M/2
                    --no-recompute-structure-matrix 
                                          compute structure matrix once and store it in 
                                          memory (time efficient but memory intensive)
                    --skip-hoffdiag-saving
                                          do not save offdiagonal terms of Hamiltonian
                                          for V=0 (if already saved or should never be 
                                          saved in combination with --skip-hoffdiag-loading)
                    --skip-hoffdiag-loading
                                          do not load offdiagonal terms of Hamiltonian
                                          from V=0
                    --no-flush            do not flush write buffer to output files in
                                          after computation for each V
                    -h, --help            show this help message and exit
                  
boundary conditions:
                    --pbc                 periodic boundary conditions (default)
                    --obc                 open boundary conditions
                  
tV parameters:
                    --V-start V_start     start V (type: Float64, default: -2.0)
                    --V-end V_end         end V (type: Float64, default: 2.0)
                    --V-step V_step       step in V (type: Float64, default: 0.1)
                    --Vp Vp               final Vp (type: Float64, default: 0.0)
                    --t t                 t value (type: Float64, default: 1.0)
                  
entanglement entropy:
                    --ee ℓ                compute all EEs with partition size ℓ (type:
                                          Int64)
"""
function main()
 # _____________1_Parameter_Setup________________
    c=parse_commandline()
    # Number of sites
    M = c[:M]
    # Number of particles
    N = c[:N]
    # Hopping
    t = c[:t]
    # Boundary conditions
    boundary = c[:boundary]
    # Size of region A
    Asize = c[:ee]
    ℓsize = div(M, 2)
    # Interaction paramers V, V' 
    V_array = c[:V_start]:c[:V_step]:c[:V_end]
    Vp = c[:Vp]
    # check parameters
    if M!=2N
        println("Not at half-filling: the number of sites =", M," and the number of particles =",N )
        exit(1)
    end

 # _____________2_Output_Setup___________________
    if c[:out] === nothing
        out_folder = "./"
    else
        out_folder = c[:out]
    end 
    if c[:tmp] === nothing
        tmp_folder = "./"
    else
        tmp_folder = c[:tmp]
    end 
    if c[:out_obdm] === nothing
        out_folder_obdm = out_folder
    else
        out_folder_obdm = c[:out_obdm]
    end
    calculation_label = @sprintf "M%02d_N%02d_t%+5.3f_Vp%+5.3f_Vsta%+5.3f_Vend%+5.3f_Vstp%+5.3f" M N t Vp c[:V_start] c[:V_end] c[:V_step]
    # 2.1. output of particle entanglement (pe_01)
        # function to convert data to string
        out_str_pe_01 = (V,data)->@sprintf "%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E\n" V data...
        # open file
        path_pe_01 = joinpath(out_folder,@sprintf "particle_entanglement_n%02d_%s.dat" Asize calculation_label)
        file_pe_01 = open(path_pe_01,"w")
        # write initial header
        write(file_pe_01, "# M=$(M), N=$(N), Vp=$(Vp), t=$(t), n=$(Asize), Vstart=$(c[:V_start]), Vstop=$(c[:V_end]), Vstep=$(c[:V_step]), $(boundary)\n")
        write(file_pe_01, "# start time $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))\n")
        write_flush(file_pe_01,@sprintf "#%24s#%24s#%24s%24s#%24s#%24s#%24s#%24s#%24s#%24s#%24s#%24s\n" "V" "S₁(n=$(Asize))" "S₂(n=$(Asize))" "S₃(n=$(Asize))" "S₄(n=$(Asize))" "S₅(n=$(Asize))" "S₆(n=$(Asize))" "S₇(n=$(Asize))" "S₈(n=$(Asize))" "S₉(n=$(Asize))" "S₁₀(n=$(Asize))" "S₀₋₅(n=$(Asize))")
 
    # 2.2. output of spatial entanglement (se_02)
    if c[:spatial]   
        # function to convert data to string
        out_str_se_02 = (V,data)->@sprintf "%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E\n" V data...
        # open file
        path_se_02 = joinpath(out_folder,@sprintf "spatial_entanglement_l%02d_%s.dat" ℓsize calculation_label)
        file_se_02 = open(path_se_02,"w")
        # write initial header
        write(file_se_02, "# M=$(M), N=$(N), Vp=$(Vp), t=$(t), l=$(ℓsize), Vstart=$(c[:V_start]), Vstop=$(c[:V_end]), Vstep=$(c[:V_step]), $(boundary)\n")
        write(file_se_02, "# start time $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))\n")
        write_flush(file_se_02,@sprintf "#%24s#%24s#%24s%24s#%24s#%24s#%24s#%24s#%24s#%24s#%24s#%24s\n" "V" "S₁(n=$(ℓsize))" "S₂(n=$(ℓsize))" "S₃(n=$(ℓsize))" "S₄(n=$(ℓsize))" "S₅(n=$(ℓsize))" "S₆(n=$(ℓsize))" "S₇(n=$(ℓsize))" "S₈(n=$(ℓsize))" "S₉(n=$(ℓsize))" "S₁₀(n=$(ℓsize))" "S₀₋₅(n=$(ℓsize))")      
    end

    # 2.3. output of pair correlations (pcf_03)
    if c[:g2]
        # function to convert data to string
        function out_str_pcf_03(V::Float64,data::Vector{Float64})::String
            str = @sprintf "%24.12E" V
            for x=1:M
                str = @sprintf "%s%24.12E" str data[x] 
            end
            str = @sprintf "%s\n" str
            return str
        end
        # open file
        path_pcf_03 = joinpath(out_folder,@sprintf "pair_correlations_g2_n%02d_%s.dat" Asize calculation_label)
        file_pcf_03 = open(path_pcf_03,"w")
        # write initial header
        write(file_pcf_03, "# M=$(M), N=$(N), Vp=$(Vp), t=$(t), n=$(Asize), Vstart=$(c[:V_start]), Vstop=$(c[:V_end]), Vstep=$(c[:V_step]), $(boundary)\n")
        write(file_pcf_03, "# start time $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))\n")
        write(file_pcf_03,@sprintf "%24s" "V" )
        for x=0:M-1
            write(file_pcf_03,@sprintf "%24d" x )
        end 
        write_flush(file_pcf_03,"\n" ) 
    end

 # _____________3_Calculation______________________ 
    # compute symmetry cycles 
    Cycle_leaders, Cycle_sizes, NumOfCycles = symmetry_cycles_q0R1PH1(M, N)
    # if requested, compute structure matrix only once here 
    if c[:no_recompute_structure_matrix]
        AmatrixStructure =PE_StructureMatrix(M, N, Cycle_leaders, Asize) 
    end
    # save off-diagonal terms once for V=0=V'
    if ~c[:skip_hoffdiag_saving] 
        ground_state(M,N,Cycle_leaders, Cycle_sizes, NumOfCycles, t, 0.0, 0.0, boundary, false,true,tmp_folder)
    end
    for V in ProgressBar(V_array)
        # compute ground state Ψ and devide by cycle size below to get Ψ_coeff (only use single variable Ψ_coeff here)
        Ψ_coeff, HRank = ground_state(M, N, Cycle_leaders, Cycle_sizes, NumOfCycles, t, V, Vp, boundary,~c[:skip_hoffdiag_loading],false,tmp_folder)
        # coefficents of basis states appearing in each cycle (renaming just for clarity) 
        for j=1: HRank
            Ψ_coeff[j]=Ψ_coeff[j]/sqrt(Cycle_sizes[j])
        end

        # 3.1 particle entanglement
        if c[:no_recompute_structure_matrix]
            if c[:obdm] && Asize==1 
                s_particle, obdm = particle_entropy_Ts(M, N, Asize, Ψ_coeff, true, AmatrixStructure) 
                # save obdm to file (one file for each V)
                save_obdm(obdm,M,N,V,Vp,out_folder_obdm)
            else
                s_particle = particle_entropy_Ts(M, N, Asize, Ψ_coeff, false, AmatrixStructure) 
            end
        else
            if c[:obdm] && Asize==1 
                s_particle, obdm = particle_entropy_Ts_and_structureMatrix(M, N, Asize, Ψ_coeff, Cycle_leaders, true) 
                # save obdm to file (one file for each V)
                save_obdm(obdm,M,N,V,Vp,out_folder_obdm)
            else
                s_particle = particle_entropy_Ts_and_structureMatrix(M, N, Asize, Ψ_coeff, Cycle_leaders, false)  
            end
        end
            # save to file
            write_flush(file_pe_01, out_str_pe_01(V,s_particle), ~c[:no_flush]) 

        # 3.2 calculate spatial entanglement 
        if c[:spatial]
            s_spatial = spatial_entropy(M, N, ℓsize, Ψ_coeff, Cycle_leaders) 
            # save to file
            write_flush(file_se_02, out_str_se_02(V,s_spatial), ~c[:no_flush])
        end

        # 3.3 pair correlation function
        if c[:g2]
            g2 = pair_correlation(M, N, Ψ_coeff, Cycle_leaders)  
            # save to file
            write_flush(file_pcf_03, out_str_pcf_03(V,g2), ~c[:no_flush])
        end
        Ψ_coeff = nothing  
    end

 # ________4_Output_Finalization___________________
    # 4.1. output of particle entanglement (pe_01)
        write_flush(file_pe_01,"\n\n Calculation finished at  $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        close(file_pe_01)    
    # 4.2. output of spatial entanglement (se_02)
    if c[:spatial]
        write_flush(file_se_02,"\n\n Calculation finished at  $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        close(file_se_02) 
    end 
    # 4.3. output of pair correlations (pcf_03)
    if c[:g2]
        write_flush(file_pcf_03,"\n\n Calculation finished at  $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        close(file_pcf_03)  
    end

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

