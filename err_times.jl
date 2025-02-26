using SparseArrays
using LinearAlgebra
using Printf

include("basis.jl")
include("basisx.jl")
include("bubble_sort.jl")
include("ChkNodeEx.jl")
include("ex_sol.jl")
include("femerror.jl")
include("f_ex.jl")
include("Fvec.jl")
include("GLpw.jl")
include("ID.jl")
include("intNodes.jl")
include("intWeights.jl")
include("KMtrx.jl")
include("Mesh2Drect.jl")
include("MMtrx.jl")


include("W2DFE_Dbc_RK4.jl")
#include("W2DFE_Dbc_RK4_LC.jl")

rmax   =   5;
hxinit = 0.5;
hyinit = 0.5;
T      = 1.0;


ndisc  =   5;

println("Change Test")
hx = hxinit./(2 .^(Array(1:ndisc) .- 1));
hy = hyinit./(2 .^(Array(1:ndisc) .- 1));

osc_f = 5;
re = zeros(rmax,ndisc);
for r = 1 : rmax

	for h  = 1 : ndisc

#		dt = (hx[h] ./2 .^(h-1))/(100*r);
		dt = hx[h]/(100*r);
		@printf("r = %5d, h = %7.6f \n",r,hx[h])
#		[re(r,h),~] = W2DFE_Dbc(r,hx(h),hy(h),T,osc_f);
		re[r,h] = W2DFE_Dbc_RK4(r,hx[h],hy[h],dt,T,osc_f)
#		re[r,h] = W2DFE_Dbc_RK4_LC(r,hx[h],hy[h],dt,T,osc_f)
		println()
		println(re[r,h])
		println()
	end
end


filename = "errors_wrt_h_r_n.dat";

f = open(filename,"w")
for r = 1 : rmax
    for h = 1 : ndisc
        print(f,re[r,h]," ")
    end
    println(f,"")
end
close(f)





#filename = "errors_wrt_h_r_n";

