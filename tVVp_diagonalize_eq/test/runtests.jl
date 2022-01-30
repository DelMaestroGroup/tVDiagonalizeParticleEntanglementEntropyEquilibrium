using Test

using Pkg

Pkg.add("Arpack")
Pkg.add("ArgParse")
Pkg.add("ProgressBars")
Pkg.add("BenchmarkTools")

@testset "All tests" begin
    
    include("test_tVVp_main_q0R1PH1_IntF.jl")

end
