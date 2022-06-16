abstract type AbstractFermionsbasis end

"""
Basis of Fermionic occupation vectors.
"""
struct Fermionsbasis <: AbstractFermionsbasis
    "Number of sites."
    K::Int
    "Number of fermions."
    N::Int

    "Number of basis vectors."
    D::Int

    "Occupation vectors (1 by D)."
    vectors::Vector{Int}
end
num_vectors(N::Int, K::Int) = binomial(K, N)
num_vectors(basis::Fermionsbasis, N::Int, K::Int) = num_vectors(N, K)

"""
    Fermionsbasis(K::Int, N::Int)
Create a basis for `K` sites and `N` fermions.
"""
function Fermionsbasis(K::Int, N::Int)
    K >= 1 || throw(DomainError(K, "At least 1 site is required."))
    N >= 0 || throw(DomainError(N, "At least 0 particles are required."))
    N <= K || throw(DomainError(N, "fermions do not fit on the sites."))

    # Basis size.
    D = num_vectors(N, K)
    v = 2^N - 1
    vectors = Vector{Int}(undef, D)

    vectors[1] = v
    for i in 2:D
        if CheckSite(v, 1) == 1
            j = findfirstEmpty(v)
            v = EmptySite(v, j - 1)
            v = OccupySite(v, j)
        else
            j = findfirstOccupide(v)
            k = findFrsEmpAfterFrsOcc(v)
            v = EmptySite(v, k - j)
            v = OccupySite(v, k)
            for l in 1:(k-j-1)
                v = OccupySite(v, l)
            end
            for l in (k-j+1):(k-1)
                v = EmptySite(v, l)
            end
        end
        vectors[i] = v
    end
    Fermionsbasis(K, N, D, vectors)
end

"""
Generator for basis states on the fly.
"""
struct BasisKetGenerator
    L::Int
    N::Int
    D::Int
end
BasisKetGenerator(L::Int, N::Int) = BasisKetGenerator(L, N, num_vectors(N, L))

struct StateTracker
    i::Int64
    last_ket::Int64
end

Base.iterate(x::BasisKetGenerator) = (2^x.N - 1, StateTracker(1, 2^x.N - 1))

function Base.iterate(x::BasisKetGenerator, state::StateTracker)
    if state.i == x.D
        return nothing
    end
    v = state.last_ket
    if CheckSite(v, 1) == 1
        j = findfirstEmpty(v)
        v = EmptySite(v, j - 1)
        v = OccupySite(v, j)
    else
        j = findfirstOccupide(v)
        k = findFrsEmpAfterFrsOcc(v)
        v = EmptySite(v, k - j)
        v = OccupySite(v, k)
        for l in 1:(k-j-1)
            v = OccupySite(v, l)
        end
        for l in (k-j+1):(k-1)
            v = EmptySite(v, l)
        end
    end
    return (v, StateTracker(state.i + 1, v))
end


"""
    serial_num(K::Int, N::Int, M::Int, v::AbstractVector{Int})
Compute the serial number of occupation vector `v` in a basis with `K` sites
and `N` fermions.
"""
function serial_num(K::Int, N::Int, v::Int)
    I = 1
    for mu in 1:K
        s = 0
        for nu in (mu+1):K
            s += CheckSite(v, Int64(nu))
        end
        for i in 0:(CheckSite(v, Int64(mu))-1)
            I += num_vectors(N - s - i, mu - 1)
        end
    end

    I
end


function serial_num(basis::Fermionsbasis, v::Int)
    return searchsortedfirst(basis.vectors, v)
end


"""
    sub_serial_num(basis::AbstractFermionsbasis, v::AbstractVector{Int})
Compute the serial number of the reduced occupation vector `v`, which has only
a subset of the sites present in `basis`.
"""
function sub_serial_num(basis::Fermionsbasis, v::Int)
    K = basis.K
    N = count_ones(v)

    # Only one way to have no sites.
    K >= 1 || return 1
    # Only one way to have no particles.
    N >= 1 || return 1

    # Count the zero-particle case.
    I = 1

    for n in 1:(N-1)
        I += num_vectors(n, K)
    end

    I + serial_num(basis, v) #serial_num(K, N, v)
end

Base.getindex(basis::AbstractFermionsbasis, i::Int) = @view basis.vectors[i]
Base.view(basis::AbstractFermionsbasis, i::Int) = reverse(bits(basis.vectors[i]))
Base.eltype(::Type{AbstractFermionsbasis}) = Int
Base.length(basis::AbstractFermionsbasis) = basis.D

function Base.in(v::Int, basis::Fermionsbasis)
    v <= 2^basis.K && count_ones(v) == basis.N
end


# find first bit=0 from the extreme right bit
function findfirstEmpty(v::Int)
    return trailing_ones(v) + 1
end

# find first bit=1 from the extreme right bit
function findfirstOccupide(v::Int)
    return trailing_zeros(v) + 1
end

# find first bit=0 after the first bit=1 from the extreme right bit
function findFrsEmpAfterFrsOcc(v::Int)
    return findfirstOccupide(v) + trailing_ones(v >>> (findfirstOccupide(v))) + 1
end

# set a bit to 1  
function OccupySite(v::Int, bit::Int64)
    return ((1 << (bit - 1)) | v)::Int64
end

# set a bit to 0 
function EmptySite(v::Int, bit::Int64)
    return convert(Int64, (~(1 << (bit - 1))) & Int64(v))
end

# set a bit to bitvalue 0 or 1  
function SetSite(v::Int, bit::Int64, bitvalue::Int64)
    return convert(Int64, (~((1 - bitvalue) << (bit - 1))) & Int64(((bitvalue << (bit - 1)) | v)))

end

# read a bit 
function CheckSite(v::Int, bit::Int64)
    return Int64((v >> (bit - 1)) & 1)
end

# Flip k bits (0 <--> 1) (101100 --> 010011)
function FlipKet(v::Int, k::Int64)
    return ((~v << (64 - k)) >>> (64 - k))::Int64
end
FlipKet(v::Int, basis::Fermionsbasis) = FlipKet(v, basis.K)

# Reverse k bits from right to lift (101100 --> 001101)
function ReverseKet(v::Int, k::Int64)
    return (Mybswap(v) >>> (64 - k))::Int64
end
ReverseKet(v::Int, basis::Fermionsbasis) = ReverseKet(v, basis.K)


function Mybswap(i::Int64)
    i = (i >>> 32) | (i << 32)
    i = ((i & -281470681808896) >>> 16) | ((i & 281470681808895) << 16)
    i = ((i & -71777214294589696) >>> 8) | ((i & 71777214294589695) << 8)
    i = ((i & -1085102592571150096) >>> 4) | ((i & 1085102592571150095) << 4)
    i = ((i & -3689348814741910324) >>> 2) | ((i & 3689348814741910323) << 2)
    return (((i & -6148914691236517206) >>> 1) | ((i & 6148914691236517205) << 1))
end

function Mybswap(i::UInt64)
    i = (i >>> 32) | (i << 32)
    i = ((i & 0xffff0000ffff0000) >>> 16) | ((i & 0x0000ffff0000ffff) << 16)
    i = ((i & 0xff00ff00ff00ff00) >>> 8) | ((i & 0x00ff00ff00ff00ff) << 8)
    i = ((i & 0xf0f0f0f0f0f0f0f0) >>> 4) | ((i & 0x0f0f0f0f0f0f0f0f) << 4)
    i = ((i & 0xcccccccccccccccc) >>> 2) | ((i & 0x3333333333333333) << 2)
    return (((i & 0xaaaaaaaaaaaaaaaa) >>> 1) | ((i & 0x5555555555555555) << 1))
end


# circle shifts k bits to the left by 1 step (101101 --> 110110)
function CircshiftKet(v::Int, k::Int64)
    return (v << (64 - k + 1) >>> (64 - k)) | (v >>> (k - 1))::Int64
end
CircshiftKet(v::Int, basis::Fermionsbasis) = CircshiftKet(v, basis.K)

#creates a subset u out of v using the indices list A
function SubKet(v::Int, A::Array{Int,1})
    u = Int64(0)
    for i = 1:length(A)
        u = SetSite(u, i, CheckSite(v, A[i]))
    end
    return convert(Int64, u)
end
SubKet(v::Int, A::Array{Int,1}, basis::Fermionsbasis) = SubKet(v, A)


