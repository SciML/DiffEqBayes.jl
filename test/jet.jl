using DiffEqBayes
using JET
using Test

@testset "JET static analysis" begin
    @testset "Package-level analysis" begin
        result = JET.report_package("DiffEqBayes"; target_modules = (DiffEqBayes,))
        @test length(JET.get_reports(result)) == 0
    end
end
