using LinearAlgebra
using Arpack
using Printf

"""Obtain a starting point for the Lanczos algorithm based on the phase of the system."""
function getΨ0_trial(t::Float64, V0::Float64, boundary::BdryCond, Rank::Int64, Cycle_sizes::Vector{Int64}, NumberOfCycles::Int64)


    if -1.95 < V0 / t < 1.95
        Ψ0_trial = ones(Float64, Rank)

    else
        Ψ0_trial = 0.01 * ones(Float64, Rank)

        if V0 / t < 0
            # attractive: cycle with largest entry wins
            Ψ0_trial[1] = 1.0 #?
        else
            # repulsive: smallest cycle wins
            i_smallest = argmin(Cycle_sizes[1:NumberOfCycles])
            Ψ0_trial[i_smallest] = 1.0 #?
        end
    end

    for j = 1:Rank
        Ψ0_trial[j] = Ψ0_trial[j] * sqrt(Cycle_sizes[j])
    end

    Ψ0_trial .= Ψ0_trial ./ sqrt(dot(Ψ0_trial, Ψ0_trial))

    return Ψ0_trial
end

"""Setup Hamiltonian q=0, P=1, R=1 block in sparse format and compute the ground state wavefunction in the symmetry basis.
also returns the rank of the Hamiltonian, i.e. length of Ψ.
"""
function ground_state(L::Int, N::Int, Cycle_leaders::Vector{Int64}, Cycle_sizes::Vector{Int64}, NumOfCycles::Int64,
    t::Float64, V::Float64, Vp::Float64, boundary::BdryCond,
    load_offdiag::Bool=false, storage_path::String="./")

    H, HRank = sparse_Block_Diagonal_Hamiltonian_q0R1PH1(L, N, Cycle_leaders, Cycle_sizes, NumOfCycles, t, V, Vp, load_offdiag, storage_path)
    #print(" sparse_hamiltonian finish\n ")
    Ψ = zeros(ComplexF64, HRank)
    Ψ .= eigs(H, nev=1, which=:SR, tol=1e-13, v0=getΨ0_trial(t, V, boundary, HRank, Cycle_sizes, NumOfCycles))[2][1:HRank] .* ones(ComplexF64, HRank)
    H = Nothing
    Ψ .= Ψ ./ sqrt(dot(Ψ, Ψ))

    return Ψ, HRank
end

struct FileHeader
    M::Int64
    N::Int64
    basis_num::Int64
    V::Float64
    Vp::Float64
end

function load_ground_state(Ψ_output::String, M::Int64, N::Int64, V::Float64, Vp::Float64, NumOfCycles::Int64)

    ######file_header1=Array{Int64}(4)
    ######file_header2=Array{Float64}(4)
    file_header1 = zeros(Int64, 3)
    file_header2 = zeros(Float64, 2)
    Ψf = open(Ψ_output, "r")
    read!(Ψf, file_header1)
    read!(Ψf, file_header2)
    M_f = file_header1[1]
    N_f = file_header1[2]
    basis_num_f = file_header1[3]
    V_f = file_header2[1]
    Vp_f = file_header2[2]
    if (M_f != M) || (N_f != N) || (abs(V_f - V) > 1.0E-12) || (abs(Vp_f - Vp) > 1.0E-12)
        println("the file of states is not compatible with the input parameters")
        println("M=", M, " N=", N, " V=", V, " Vp=", Vp)
        println("M_f=", M_f, " N_f=", N_f, " V_f=", V_f, " Vp_f=", Vp_f)
        exit(1)
    end
    Ψ = zeros(ComplexF64, NumOfCycles)
    read!(Ψf, Ψ)
    close(Ψf)
    HRank = basis_num_f

    return Ψ, HRank
end

function save_ground_state(Ψ_output::String, Ψ::Vector{ComplexF64}, M::Int64, N::Int64, HRank::Int64, V::Float64, Vp::Float64)
    file_header = FileHeader(M, N, HRank, V, Vp)
    Ψf = open(Ψ_output, "w")
    write(Ψf, file_header.M, file_header.N, file_header.basis_num, file_header.V, file_header.Vp, Ψ)
    flush(Ψf)
    close(Ψf)
end