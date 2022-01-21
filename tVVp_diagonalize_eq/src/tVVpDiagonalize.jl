module tVVpDiagonalize

using IntFermionicbasis
using SparseArrays
using LinearAlgebra: svdvals!,Symmetric, svdvals

export
    BdryCond,
    OBC,
    PBC,
    sparse_hamiltonian,
    spatial_entropy,
    particle_entropy_Ts,
    PE_StructureMatrix,
    Translational_Symmetry_Cycles,
    Symmetry_Cycles_q0R1PH1,
    Block_Diagonal_Hamiltonian_q0R1PH1,
    sparse_Block_Diagonal_Hamiltonian_q0R1PH1,
#   full_hamiltonian, 
    pair_correlation,
    ground_state,
    load_ground_state,
    save_ground_state,
    save_obdm
"""
Boundary conditions.
"""
@enum BdryCond PBC OBC
@doc "Periodic boundary conditions." PBC
@doc "Open boundary conditions." OBC

include("sparse_hamiltonian.jl")
include("spatial_entropy.jl")
include("particle_entropy_Ts.jl")
include("PE_StructureMatrix.jl")
include("Translational_Symmetry_Cycles.jl")
include("Symmetry_Cycles_q0R1PH1.jl")
include("Block_Diagonal_Hamiltonian_q0R1PH1.jl")
include("sparse_Block_Diagonal_Hamiltonian_q0R1PH1.jl")
#include("full_hamiltonian.jl") 
include("pair_correlations.jl")
include("ground_state.jl")
end
