using Pkg , Printf
Pkg.add("JuMP")
Pkg.add("GLPK")
#Pkg.add("GLPKMathProgInterface")
Pkg.add("Plots")
Pkg.add("Distributions")
Pkg.add("SpecialFunctions")
using GLPK
using JuMP, GLPKMathProgInterface, Plots, Distributions, SpecialFunctions
include("../src/OT_sampler.jl")
plot_flag = 1
savedata_flag = 0
a = @time OT_sampler(plot_flag,savedata_flag)
# write your own tests here
# @test 1 == 1
