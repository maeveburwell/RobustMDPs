

@testset "Simple L1" begin
    # Test case 1
    z = [0.5, 0.2, 0.9, 0.1]
    p̄ = [0.25, 0.25, 0.25, 0.25]
    ξ = 0.5

    # Run the function
    popt, obj = worstcase_l1(z, p̄, ξ)

    # Display the results
    #println("Test case 1")
    #println("Optimal p: ", popt)
    #println("Objective value: ", obj)

    @test obj ≤ z' * p̄d
    @test obj ≥ minimum(z)
end

