using LinearAlgebra

"""Obtain a starting point for the Lanczos algorithm based on the phase of the system."""
function getΨ0_trial(t::Float64, V0::Float64, boundary::BdryCond, basis::AbstractFermionsbasis, Rank::Int64, CycleSize::Vector{Int64},InvCycles_Id::Vector{Int64})


    if -1.95 < V0/t < 1.95
        Ψ0_trial = ones(Float64, Rank)

    else
        Ψ0_trial = 0.01*ones(Float64, Rank)
        
        if V0/t < 0 
            # attractive: cycle with largest entry wins
            Ψ0_trial[1] = 1.0 #?
        else
            # repulsive: smallest cycle wins
            i_smallest = argmin(CycleSize[CycleSize.>0]) 
            Ψ0_trial[i_smallest] = 1.0 #?
        end 
    end

    for j=1: Rank
       Ψ0_trial[j]= Ψ0_trial[j]*sqrt(CycleSize[j])
    end

    Ψ0_trial.=Ψ0_trial./sqrt(dot(Ψ0_trial,Ψ0_trial))

    return Ψ0_trial
end