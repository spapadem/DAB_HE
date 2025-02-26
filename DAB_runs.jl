include("basis.jl")
include("basisx.jl")
include("bubble_sort.jl")
include("ChkNodeEx.jl")
include("ex_sol.jl")
include("femerror2.jl")
include("f_ex.jl")
include("Fvec.jl")
include("GLpw.jl")
include("ID.jl")
include("intNodes.jl")
include("intWeights.jl")
include("Mesh2Drect.jl")
include("norm_calc.jl")
include("I1.jl")
include("I2.jl")
include("I3.jl")
include("I4.jl")
include("I5.jl")
include("I6.jl")
include("green.jl")
include("HE2D_mult.jl")
using SparseArrays
using LinearAlgebra
using Printf



for r = 1 : 5
	for P = 1 : 5
		println();
		println("r = ",r,", P = ",P)
		println();
		HE2D(r,P)
	end
end
