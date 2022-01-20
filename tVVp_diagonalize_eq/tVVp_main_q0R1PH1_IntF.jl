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
    time_num:: Int64
    basis_num:: Int64
    V0::Float64
    V::Float64
    Vp0::Float64
    Vp::Float64
    time_min::Float64
    time_max::Float64
end

# ------------------------------------------------------------------------------
function getΨ0_trial(t::Float64, V0::Float64, boundary::BdryCond, basis::AbstractFermionsbasis, Rank::Int64, CycleSize::Vector{Int64},InvCycles_Id::Vector{Int64})


    if -1.95 < V0/t < 1.95
        Ψ0_trial = ones(Float64, Rank)

    else
        Ψ0_trial = 0.01*ones(Float64, Rank)

        if boundary==OBC
            num_link = basis.K-1
        elseif boundary==PBC
            num_link = basis.K
        end

        for i=1: Rank
            bra=basis.vectors[Cycles[i,1]]
            cc=0
            for j=1:num_link j_next = j % basis.K + 1
		cc+=CheckSite(bra,j)*CheckSite(bra,j_next)
            end

            if (cc== basis.N-1) && (V0/t < -1.95)
                Ψ0_trial[i]=1.0
            elseif (cc==0) && (V0/t > 1.95)
                Ψ0_trial[i]=1.0
            end
        end
    end

    for j=1: Rank
       Ψ0_trial[j]= Ψ0_trial[j]*sqrt(CycleSize[j])
    end

    Ψ0_trial.=Ψ0_trial./sqrt(dot(Ψ0_trial,Ψ0_trial))

    return Ψ0_trial
end

# ------------------------------------------------------------------------------
function pair_correlation(basis::AbstractFermionsbasis, d::Vector{ComplexF64}, InvCycles_Id::Vector{Int64})
""" We exploit translational symmetry to speed things up. """

    # setup the needed vectors
    g2 = zeros(Float64, basis.K)

    # measure the pair correlation function (density-density) exploiting
    # translational symmetry
    for i=1:basis.K
        for b=1:basis.D
            config = basis.vectors[b]
            weight = d[b]
            n0 = CheckSite(config,1)
            ni = CheckSite(config,i)
            g2[i] += n0*ni*abs(weight)^2
        end
    end

    return g2
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
    "--Allqs"
        help = "must be used, if the initial state is not an eigenstate of the one site translation operator with a unity eigenvalue (q=0)"
        action = :store_true
    "--Allps"
        help = "must be used, if the initial state is not eigenstate of the particle-hole-exchange operator with a unity eigenvalue (p=1)"
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
    "--V0"
        metavar = "V0"
        help = "initial V"
        arg_type = Float64
        default = 0.0
    "--V"
        metavar = "V"
        help = "final V"
        arg_type = Float64
        default = 1.0
    "--Vp0"
        metavar = "Vp0"
        help = "initial Vp"
        arg_type = Float64
        default = 0.0
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
add_arg_group(s, "time parameters")
@add_arg_table s begin
    "--time-min"
        metavar = "time"
        help = "minimum time"
        arg_type = Float64
        default = 0.0
    "--time-max"
        metavar = "time"
        help = "maximum time"
        arg_type = Float64
        default = 5.0
    "--time-step"
        metavar = "time"
        help = "time step"
        arg_type = Float64
        default = 0.1
    "--time-num"
        metavar = "N"
        help = "number of time"
        arg_type = Int
    "--time-log"
        help = "use logarithmic scale for time"
        action = :store_true
    "--ftime-min"
        metavar = "time"
        help = "file minimum time"
        arg_type = Float64
    "--ftime-max"
        metavar = "time"
        help = "file maximum time"
        arg_type = Float64

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

# Initial V
V0 = c[:V0]
Vp0 = c[:Vp0]

# Final V
V = c[:V]
Vp = c[:Vp]

# Initial time
time_min = c[:time_min]

if c[:time_log] && c[:time_num] === nothing
    println("--time-log must be used with --time-num")
    exit(1)
end

if c[:time_step] === nothing
    if c[:time_num] === nothing
        time_range = c[:time_min]:0.5:c[:time_max]
        time_num=length(time_range) 
        Δt = 0.5
    else
        if c[:time_log]
            time_range = logspace(c[:time_min], c[:time_max], c[:time_num])
        else
            time_range = linspace(c[:time_min], c[:time_max], c[:time_num])
            if length(time_range) > 1
                Δt = time_range[2]-time_range[1]
            else
                Δt = time_range[1]
            end
        end
        time_num=c[:time_num]
    end
else
    if c[:time_num] === nothing
        time_range = c[:time_min]:c[:time_step]:c[:time_max]
        time_num=length(time_range) 
        Δt = c[:time_step]
    else
        println("--time-step and --time-num may not both be supplied")
        exit(1)
    end
end

# Output file
if c[:out] === nothing
     output = @sprintf "partEE_%02d_%02d_%+5.3f_%+5.3f_%+5.3f_%+5.3f_%6.4f_%06.3f_%06.3f_%1d.dat" M N V0 Vp0 V Vp Δt time_range[1] time_range[end] Asize
else
      output = c[:out]
 end

# output file if we are measuring the spatial entanglement entropy
 if c[:spatial]
     spat_output = @sprintf "spatEE_%02d_%02d_%+5.3f_%+5.3f_%+5.3f_%+5.3f_%6.4f_%06.3f_%06.3f_%1d.dat" M N V0 Vp0 V Vp Δt time_range[1] time_range[end] Asize
 end

# output file if we are measuring the pair correlation function 
if c[:g2]
     g2_output = @sprintf "g2_%02d_%02d_%+5.3f_%+5.3f_%+5.3f_%+5.3f_%6.4f_%06.3f_%06.3f_%1d.dat" M N V0 Vp0 V Vp Δt time_range[1] time_range[end] Asize
end

# state output file
if c[:save_states] || c[:load_states]
    if c[:states_file] === nothing
        if (c[:load_states])
            if c[:ftime_max] === nothing
                 ftime_max= time_range[end]
            else
                 ftime_max= c[:ftime_max]
            end 
            if c[:ftime_min] === nothing
                 ftime_min= time_range[1]
            else
                 ftime_min= c[:ftime_min]
            end
        elseif (c[:save_states])
            ftime_max= time_range[end]
            ftime_min= time_range[1]
        end
        Ψt_output = @sprintf "psioft_%02d_%02d_%+5.3f_%+5.3f_%+5.3f_%+5.3f_%6.4f_%06.3f_%06.3f.dat" M N V0 Vp0 V Vp Δt ftime_min ftime_max
    else
        Ψt_output=c[:states_file]
    end
end
basis = Fermionsbasis(M, N)

q0 = c[:Allqs] ? 1 : 0
p1 = c[:Allps] ? 1 : -1

#_______________________________________________________________________________
ll=length(basis)
μ=zeros(Float64, M)
exp_q=zeros(ComplexF64, basis.K)
#  Initial wave function in terms of the spatial basis
Ψn=ComplexF64
# the one body density matrix
obdm=zeros(Float64,M,length(time_range))

#_______________________________________________________________________________

# open and prepare files for output
if ~c[:save_states]
     f_part = open(output, "w")
     write(f_part, "# M=$(M), N=$(N), V0=$(V0), Vp0=$(Vp0), V=$(V), Vp=$(Vp), $(boundary)\n")
     write(f_part,@sprintf "#%11s%24s%24s\n" "time (tJ)" "S₁(n=$(Asize))" "S₂(n=$(Asize))")
     if c[:spatial]
        ℓsize = div(M, 2)
        f_spat = open(spat_output, "w")
        write(f_spat, "# M=$(M), N=$(N), V0=$(V0), Vp0=$(Vp0), V=$(V), Vp=$(Vp), $(boundary)\n")
        write(f_spat,@sprintf "#%11s%24s%24s\n" "time (tJ)" "S₁(ℓ=$(ℓsize))" "S₂(ℓ=$(ℓsize))")
     end
    if c[:g2]
        f_g2 = open(g2_output, "w")
        write(f_g2, "# M=$(M), N=$(N), V0=$(V0), Vp0=$(Vp0), V=$(V), Vp=$(Vp), $(boundary)\n")
        write(f_g2,@sprintf "#%11s" "time (tJ)")
        for x=0:M-1
            write(f_g2,@sprintf "%24d" x )
        end
        write(f_g2,"\n")
    end

end
#_______________________________________________________________________________

# Exploit symmetries of the hamiltonian to perform a bloack diagonalization
Cycles, CycleSize, NumOfCycles, InvCycles_Id, InvCycles_order =Symmetry_Cycles_q0R1PH1(basis)
      #println("size(Cycles) = ",Base.summarysize(Cycles)/1024^3," gb")

if ~c[:load_states]  
   #---------------------------------------find the states---------------------------------------
   # Create the Hamiltonian 
   #H = sparse_hamiltonian(basis,c[:t],V0,Vp0 ,boundary=boundary) 
   H, HRank = sparse_Block_Diagonal_Hamiltonian_q0R1PH1(basis, Cycles, CycleSize, NumOfCycles, InvCycles_Id, InvCycles_order, c[:t],V0,Vp0) 
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

   Ψ = eigs(H, nev=1, which=:SR,tol=1e-13,v0=getΨ0_trial(c[:t],V0,boundary,basis, HRank, CycleSize, InvCycles_Id))[2][1: HRank].*ones(ComplexF64, HRank)
   #Ψ = eigs(H, nev=1, which=:SR,tol=1e-13,v0=getΨ0_trial(c[:t],V0,boundary,basis, HRank,))[2][1: HRank].*ones(ComplexF64, HRank)
   #println("size(complex) = ", Base.summarysize(Ψ[1]))
   H= Nothing
   Ψ.= Ψ./sqrt(dot(Ψ,Ψ))

   if abs(time_range[1])> 1.0E-12 || length(time_range) >1
      EigenEnergie=ComplexF64
      #  wave function in terms of the Symmetry basis at time t
      Ψt=zeros(ComplexF64, NumOfCycles, time_num)
      TimeEvolutionFactor=zeros(ComplexF64, time_num)

      for q =0: (basis.K-1)*q0
         for P=1:-2:-1*p1 
            #Create the Hamiltonian

            #H0 = full_hamiltonian(basis, c[:t], V, Vp,boundary=boundary)
            Hq,HqRank = Block_Diagonal_Hamiltonian_q0R1PH1(basis, Cycles, CycleSize, NumOfCycles, InvCycles_Id, InvCycles_order, c[:t],V,Vp) 
            EigenEnergies,VV = eigen(Symmetric(Hq))
            #EigenEnergies,VV = eigen(Hq)
            println("size(Hq) = ",Base.summarysize(Hq)/1024^3," gb")
            Hq= Nothing
            for i_HqRank =1: HqRank
               EigenEnergie= EigenEnergies[i_HqRank]
               Ψn =dot(VV[:, i_HqRank],Ψ)
               TimeEvolutionFactor.=exp.(-(0.0+1.0im)* time_range*EigenEnergie)*Ψn
               Ψt+=kron(VV[:, i_HqRank],transpose(TimeEvolutionFactor))
            end
            VV= Nothing
            EigenEnergies= Nothing
            Cycles= Nothing
            InvCycles_order= Nothing
         end
      end   
      print(" Block_Diagonal_Hamiltonian finished\n ")
   else
      Ψt=zeros(ComplexF64, NumOfCycles, time_num)
      Ψt[:,1]=Ψ[:]
   end
else  
   #---------------------------------------load the states---------------------------------------
   Cycles= Nothing
   InvCycles_order= Nothing
   ######file_header1=Array{Int64}(4)
   ######file_header2=Array{Float64}(4)
   file_header1 =zeros(Int64,4)
   file_header2 =zeros(Float64,6)
   Ψtf=open(Ψt_output, "r")
      read!(Ψtf,file_header1)
      read!(Ψtf,file_header2)
      M_f=file_header1[1]
      N_f=file_header1[2]
      time_num_f=file_header1[3]
      basis_num_f=file_header1[4]
      V0_f=file_header2[1]
      V_f=file_header2[2]
      Vp0_f=file_header2[3]
      Vp_f=file_header2[4]
      time_min_f=file_header2[5]
      time_max_f=file_header2[6]
      time_range_f = LinRange(time_min_f, time_max_f, time_num_f)
      println(length(time_range_f))
      if length(time_range_f) > 1
          Δt_f = time_range_f[2]-time_range_f[1]
      else
          Δt_f = time_range_f[1]
      end
      println(Δt_f)
      if (M_f!=M) || (N_f!=N) || (abs(Δt_f-Δt)> 1.0E-12)||(abs(V0_f- V0)> 1.0E-12) ||(abs(V_f- V)> 1.0E-12)||(abs(Vp0_f- Vp0)> 1.0E-12)||(abs(Vp_f- Vp)> 1.0E-12) ||((time_range_f[1]- time_range[1])> 1.0E-12) ||(time_range_f[end]- time_range[end]< -1.0E-12)
         println("the file of states is not compatible with the input parameters" )
         println("M=",M," N=",N," V0=",V0," V=",V," Vp0=",Vp0," Vp=",Vp," Δt=",Δt," time_max=", time_range[end]," time_min=",time_min)
         println("M_f=",M_f," N_f=",N_f," V0_f=",V0_f," V_f=",V_f," Vp0_f=",Vp0_f," Vp_f=",Vp_f," Δt_f=",Δt_f," time_max_f=", time_max_f," time_min_f=",time_min_f)
         exit(1)
      end
      Ψt = zeros(ComplexF64, NumOfCycles, time_num)
      Ψ=zeros(ComplexF64, NumOfCycles)
      it=0
      for (it_f, time_f) in enumerate(time_range_f)
        read!(Ψtf, Ψ)
        if (abs(time_f- time_range[it+1])< 1.0E-12)
            it+=1
            Ψt[:,it]=Ψ[:]
            if it== time_num
                break
            end
        end
      end
   close(Ψtf)
   HRank = basis_num_f
   HqRank = basis_num_f
end
if ~c[:save_states] 
   #---------------------------------------calculate the entanglement---------------------------------------
   AmatrixStructure =PE_StructureMatrix(basis, Asize, InvCycles_Id)

   it = 1

   for time in time_range[1:end] 

      Ψ.=Ψt[:,it]./sqrt(dot(Ψt[:,it],Ψt[:,it]))
      for j=1: HRank
         Ψ[j]=Ψ[j]/sqrt(CycleSize[j])
      end

       if c[:spatial]
          s_spatial = spatial_entropy(basis, ℓsize, Ψ, InvCycles_Id)
          write(f_spat, @sprintf "%12.6f%24.12E%24.12E\n" time s_spatial[1] s_spatial[2])
          flush(f_spat)
       end

      # measure the pair correlation function
      if c[:g2]
          g2 = pair_correlation(basis,Ψ, InvCycles_Id)
          write(f_g2, @sprintf "%12.6f" time)
          for x=1:M
              write(f_g2, @sprintf "%24.12E" g2[x])
          end
          write(f_g2,"\n")
          flush(f_g2)
      end

       if c[:obdm] && Asize == 1
          s_particle,obdm[:,it] = particle_entropy_Ts(basis, Asize, Ψ,c[:obdm], AmatrixStructure)

       else
         s_particle = particle_entropy_Ts(basis, Asize, Ψ,c[:obdm], AmatrixStructure)

       end

       write(f_part, @sprintf "%12.6f%24.12E%24.12E\n" time s_particle[1] s_particle[2])
       flush(f_part)
      it += 1
   end
    close(f_part)

    if c[:spatial]
       close(f_spat)
    end

   # close the pair correlation function file
   if c[:g2]
       close(f_g2)
   end

#_______________________________________________________________________________
   # output the time dependent OBDM to disk
   if c[:obdm] && Asize == 1
       obdm_name = @sprintf "obdm_%02d_%02d_%+5.3f_%+5.3f_%+5.3f_%+5.3f_%6.4f_%06.3f_%06.3f.dat" M N V0 Vp0 V Vp Δt time_range[1] time_range[end]
         obdm_f = open(obdm_name, "w")
           write(obdm_f, @sprintf "#%11s" "|i-j|")
           for time in time_range
              write(obdm_f, @sprintf "%16.6f" time)
           end
           write(obdm_f, "\n")
           flush(obdm_f)

           for i = 1:M
               write(obdm_f, @sprintf "%16d" (i-Int(M/2)))
               for (time_index, time) in enumerate(time_range)
                   write(obdm_f, @sprintf "%16.6E" obdm[i,time_index])
               end
               write(obdm_f, "\n")
               flush(obdm_f)
           end
        close(obdm_f)
   end
else 
   #---------------------------------------save the states---------------------------------------

   file_header= FileHeader(M, N, time_num, HRank, V0, V, Vp0, Vp, c[:time_min], c[:time_max])
   Ψtf=open(Ψt_output, "w")
      write(Ψtf, file_header.M, file_header.N, file_header.time_num, file_header.basis_num, file_header.V0, file_header.V, file_header.Vp0, file_header.Vp ,file_header.time_min, file_header.time_max,Ψt)
      flush(Ψtf)
   close(Ψtf)
end

end

main()

