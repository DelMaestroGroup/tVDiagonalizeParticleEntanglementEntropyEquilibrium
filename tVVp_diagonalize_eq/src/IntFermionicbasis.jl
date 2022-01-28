module IntFermionicbasis

export
    AbstractFermionsbasis,
    Fermionsbasis,
    num_vectors,
    serial_num,
    sub_serial_num,
    CheckSite,
    OccupySite,
    EmptySite,
    FlipKet,
    ReverseKet,
    CircshiftKet,
    SubKet,
    BasisKetGenerator

include("IntFermionsbasis_functions.jl")

end
