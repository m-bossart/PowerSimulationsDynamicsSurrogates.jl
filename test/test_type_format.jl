@testset "Test SurrogatePerturbations format" begin
    for T in subtypes(SurrogatePerturbation)
        @test T() !== nothing
        @test in(:type, fieldnames(T))
    end
end

@testset "Test SurrogateOperatingPoint format" begin
    for T in subtypes(SurrogateOperatingPoint)
        @test T() !== nothing
        @test in(:type, fieldnames(T))
    end
end

@testset "Test SurrogateDataset format" begin
    for T in subtypes(SurrogateDataset)
        @test T() !== nothing
        @test in(:type, fieldnames(T))
    end
end
