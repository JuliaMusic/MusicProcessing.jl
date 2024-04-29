using TestItemRunner

@testitem "audio.jl" begin include("audio.jl") end

@run_package_tests verbose=true
