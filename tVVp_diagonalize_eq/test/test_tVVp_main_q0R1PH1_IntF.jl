""" Definitions of the unit tests run test by running "runtests.jl" """
# provide additional print output during tests
const VERBOSE = false

# import modules used in tests  
push!(LOAD_PATH, joinpath(dirname(@__FILE__), "..", "src"))
using Test 

using tVVpDiagonalize  
using ArgParse
using IntFermionicbasis 
using Printf
using LinearAlgebra

test_save_location = "./test/"

# test helper functions
function test_approx(got, truth, atol, verbose=true)
    passed = isapprox(got , truth , atol=atol)
    if verbose
        if passed
            print("\nPASS: ")
        else
            print("\nFAIL: ")
        end
        print("got: ", got,", truth: ", truth, '\n')
    end
    @test passed
end

# define test cases
@time @testset "Test Entropy and Pair Correlations" begin
    
    """
    Runs test cases for spatial entanglement entropy and particle entanglement entropy
    comparing with S1_spatial_truth, S2_spatial_truth, S1_particle_truth, S2_particle_truth
    allowing for a tolerance tol.
    """
    function tests_entropies(t::Float64,V::Float64,Vp::Float64,M::Int64,N::Int64,n::Int64,l::Int64,S1_spatial_truth::Float64, S2_spatial_truth::Float64, S1_particle_truth::Float64, S2_particle_truth::Float64,pair_correlations_truth::Vector{Float64},boundary::BdryCond=PBC, tol::Float64=0.0, load_offdiag::Bool=false, save_offdiag::Bool=false; verbose=VERBOSE)  
        # init 
        obdm=zeros(Float64,M) 
        # compute cycles 
        Cycle_leaders, Cycle_sizes, NumOfCycles = symmetry_cycles_q0R1PH1(M, N)
        # compute ground state 
        if save_offdiag
            save_offdiag_sparse_Block_Diagonal_Hamiltonian_q0R1PH1(M, N ,Cycle_leaders , Cycle_sizes, NumOfCycles, t,test_save_location)
        end
        Ψ_coeff, HRank = ground_state(M, N, Cycle_leaders, Cycle_sizes, NumOfCycles, t, V, Vp, boundary,load_offdiag,test_save_location) 
        # coefficents of basis states appearing in each cycle (renaming just for clarity) 
        for j=1: HRank
            Ψ_coeff[j]=Ψ_coeff[j]/sqrt(Cycle_sizes[j])
        end
        # calculate spatial entanglement
        ℓsize = div(M, 2)
        s_spatial = spatial_entropy(M, N, ℓsize, Ψ_coeff, Cycle_leaders) 
        # pair correlation function
        g2 = pair_correlation(M,N,Ψ_coeff,Cycle_leaders) 
        test_approx(g2, pair_correlations_truth, tol, verbose)
        # structure matrix
        if n==1
            s_particle, obdm = particle_entropy_Ts_and_structureMatrix(M, N, n, Ψ_coeff, Cycle_leaders, true) 
        else
            s_particle = particle_entropy_Ts_and_structureMatrix(M, N, n, Ψ_coeff, Cycle_leaders, false) 
        end

        # test against production code result 
        if verbose
            print('\n')
        end
        test_approx(s_spatial[1], S1_spatial_truth, tol, verbose)
        test_approx(s_spatial[2], S2_spatial_truth, tol, verbose)
        test_approx(s_particle[1], S1_particle_truth, tol, verbose)
        test_approx(s_particle[2], S2_particle_truth, tol, verbose)   
        if verbose
            print('\n')
        end
    end 
    
    testcases = Vector()

    # test case 1:
    params = (
            t=0.8,
            V=1.3,
            Vp=0.4,
            M=12,
            N=6,
            n=3,
            l=6,
            S1_spatial_truth= 1.180537978291E+00, 
            S2_spatial_truth= 8.585299146596E-01, 
            S1_particle_truth=1.261431311869E-01, 
            S2_particle_truth=5.137504874197E-02, 
            pair_correlations_truth= [5.000000000000E-01,1.247136397872E-01,2.719361733294E-01,2.275256708928E-01,2.612951246071E-01,2.345422692789E-01,2.599742442092E-01,2.345422692789E-01,2.612951246071E-01,2.275256708928E-01,2.719361733294E-01,1.247136397872E-01],
            boundary=PBC,
            tol=0.0)
    push!(testcases,params)
    # test case 2 
    params = (
            t=2.3,
            V=0.3,
            Vp=1.3,
            M=14,
            N=7,
            n=3,
            l=7,
            S1_spatial_truth=1.238733545053E+00, 
            S2_spatial_truth=9.893307250029E-01, 
            S1_particle_truth=5.064243992800E-02, 
            S2_particle_truth=1.482187045819E-02, 
            pair_correlations_truth= [5.000000000000E-01,1.571939107499E-01,2.310417986060E-01,2.453536553746E-01,2.493565776083E-01,2.448472790070E-01,2.491579823775E-01,2.460975925536E-01,2.491579823775E-01,2.448472790070E-01,2.493565776083E-01,2.453536553746E-01,2.310417986060E-01,1.571939107499E-01],
            boundary=PBC,
            tol=0.0)
    push!(testcases,params)
    # test case 3
    params = (
            t=0.1,
            V=0.003,
            Vp=0.0,
            M=14,
            N=7,
            n=1,
            l=7,
            S1_spatial_truth=1.226149079417E+00, 
            S2_spatial_truth=9.781005910398E-01, 
            S1_particle_truth=1.688303810217E-04, 
            S2_particle_truth=2.862259822134E-05,
            pair_correlations_truth= [5.000000000000E-01,1.462109470708E-01,2.507619036892E-01,2.366314147620E-01,2.502528490923E-01,2.435806164301E-01,2.501697681511E-01,2.447850016088E-01,2.501697681511E-01,2.435806164301E-01,2.502528490923E-01,2.366314147620E-01,2.507619036892E-01,1.462109470708E-01],
            boundary=PBC, 
            tol=0.0)
    push!(testcases,params)
    # test case 4
    params = (
            t=1.0,
            V=100.0,
            Vp=0.0,
            M=14,
            N=7,
            n=5,
            l=7,
            S1_spatial_truth=6.951898117047E-01, 
            S2_spatial_truth=6.935472606275E-01, 
            S1_particle_truth=6.931479839384E-01, 
            S2_particle_truth=6.931465943021E-01,
            pair_correlations_truth= [5.000000000000E-01,9.997000789381E-05,4.998000699802E-01,1.999500076619E-04,4.998000200165E-01,1.999799967261E-04,4.998000199821E-01,1.999800178869E-04,4.998000199821E-01,1.999799967261E-04,4.998000200165E-01,1.999500076619E-04,4.998000699802E-01,9.997000789381E-05],
            boundary=PBC, 
            tol=0.0)
    push!(testcases,params)
    # test case 5
    params = (
            t=0.1,
            V=0.003,
            Vp=0.0,
            M=14,
            N=7,
            n=1,
            l=7,
            S1_spatial_truth=1.226149079417E+00, 
            S2_spatial_truth=9.781005910398E-01, 
            S1_particle_truth=1.688303810217E-04, 
            S2_particle_truth=2.862259822134E-05,
            pair_correlations_truth= [5.000000000000E-01,1.462109470708E-01,2.507619036892E-01,2.366314147620E-01,2.502528490923E-01,2.435806164301E-01,2.501697681511E-01,2.447850016088E-01,2.501697681511E-01,2.435806164301E-01,2.502528490923E-01,2.366314147620E-01,2.507619036892E-01,1.462109470708E-01],
            boundary=PBC, 
            tol=0.0,
            load_offdiag=false,
            save_offdiag=true)
    push!(testcases,params)
    # test case 6
    params = (
            t=0.1,
            V=0.003,
            Vp=0.0,
            M=14,
            N=7,
            n=1,
            l=7,
            S1_spatial_truth=1.226149079417E+00, 
            S2_spatial_truth=9.781005910398E-01, 
            S1_particle_truth=1.688303810217E-04, 
            S2_particle_truth=2.862259822134E-05,
            pair_correlations_truth= [5.000000000000E-01,1.462109470708E-01,2.507619036892E-01,2.366314147620E-01,2.502528490923E-01,2.435806164301E-01,2.501697681511E-01,2.447850016088E-01,2.501697681511E-01,2.435806164301E-01,2.502528490923E-01,2.366314147620E-01,2.507619036892E-01,1.462109470708E-01],
            boundary=PBC, 
            tol=0.0,
            load_offdiag=true,
            save_offdiag=false)
    push!(testcases,params)
    # test case 7
    params = (
            t=1.0,
            V=100.0,
            Vp=0.0,
            M=14,
            N=7,
            n=5,
            l=7,
            S1_spatial_truth=6.951898117047E-01, 
            S2_spatial_truth=6.935472606275E-01, 
            S1_particle_truth=6.931479839384E-01, 
            S2_particle_truth=6.931465943021E-01,
            pair_correlations_truth= [5.000000000000E-01,9.997000789381E-05,4.998000699802E-01,1.999500076619E-04,4.998000200165E-01,1.999799967261E-04,4.998000199821E-01,1.999800178869E-04,4.998000199821E-01,1.999799967261E-04,4.998000200165E-01,1.999500076619E-04,4.998000699802E-01,9.997000789381E-05],
            boundary=PBC, 
            tol=0.0,
            load_offdiag=true,
            save_offdiag=false)
    push!(testcases,params)
    # test case 8
    params = (
            t=2.3,
            V=0.3,
            Vp=1.3,
            M=14,
            N=7,
            n=3,
            l=7,
            S1_spatial_truth=1.238733545053E+00, 
            S2_spatial_truth=9.893307250029E-01, 
            S1_particle_truth=5.064243992800E-02, 
            S2_particle_truth=1.482187045819E-02, 
            pair_correlations_truth= [5.000000000000E-01,1.571939107499E-01,2.310417986060E-01,2.453536553746E-01,2.493565776083E-01,2.448472790070E-01,2.491579823775E-01,2.460975925536E-01,2.491579823775E-01,2.448472790070E-01,2.493565776083E-01,2.453536553746E-01,2.310417986060E-01,1.571939107499E-01],
            boundary=PBC,
            tol=0.0)
    push!(testcases,params)

    # run tests 
    @testset "Test case $i" for i in 1:length(testcases)
            tests_entropies(testcases[i]...) 
    end
end

# define test cases
@time @testset "Test Serial Number" begin
    
    """
    Runs test cases for the serial_num function.
    """
    basis = Fermionsbasis(12,6)
    @test serial_num(12,6,63)==1
    @test serial_num(basis,63)==1
    @test serial_num(12,6,95)==2
    @test serial_num(basis,95)==2
    @test serial_num(12,6,119)==4
    @test serial_num(basis,119)==4
    basis = Fermionsbasis(30,15)
    @test serial_num(30,15,1073680384)==155117517
    @test serial_num(basis,1073680384)==155117517
    @test serial_num(30,15,65471)==10
    @test serial_num(basis,65471)==10
end