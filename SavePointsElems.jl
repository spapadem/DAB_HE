using SparseArrays
using LinearAlgebra
using Printf

#function  W2DFE_Dbc_RK4(r,hx,hy,dt,T,osc_f)

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
include("I3.jl")
include("I4.jl")
include("I5.jl")
include("I6.jl")
#close all


#################################################################
#						Problem Setup							#
#################################################################
# Create points for mesh generation
hx = 0.05;
hy =  hx;
#Waveguide characteristics.
b  = 3;  # Width.
xI = 5;  # Location of interface.
xE = 6; # Location of ABC.
nL = (xE-xI)/hx;
#xF = 3;  # Far boundary.
c  = sqrt(3);  # Wave speed.
s  = 0;
P  = 5;  # Order of ABC.
#Problem setup
r = 4;
osc_f =1;
x  = Array(0 : hx : xE);    
y  = Array(0 : hy : b );
Nx = length(x);
Ny = length(y);
N  = Nx*Ny;
# Time discretization setup
#dt = min(hx,hy)/(100*r);          % Time step
T = 30;
dt = hx/(100*r);
t = Array(0 : dt :T);      # Time vector

SaveSol = 1; # Option if we wish to save the solution in each time step to a file.

# Gaussian pulse as an initial condition. Can use as many sources as we want,
# if we supply their locations in u_0 (see below).
radius = 1.0;   				# Radius of the Gaussian pulse.
a = log(10^6)/(radius*radius);  # Amplitude.
xstar =  1.5;						# Location in x-direction.
ystar =  1.5;						# Location in y-direction.
xstar2 = 8;						# Location in x-direction for 2nd source.
ystar2 = 7.5;					# Location in y-direction for 2nd source

u0(x,y,xstar,ystar)=sqrt(2*a/pi)*exp.(-a*((x .-xstar).^2+(y .-ystar).^2)) .*
					((((x.-xstar).^2+(y.-ystar).^2)).<=radius);



# A more complex initial condition.
#u0(x,y) = ((x .- 8.5).^2 .- 1).^2 .*sin.(2*pi*y/b).*(x .>= 7.5).*(x .<= 9.5);

# A Gaussian in x and sin function in y
#u0(x,y,xstar,ystar) = sqrt(2*a/pi)*exp.(-a*((x .-xstar).^2)).*sin.(2*pi*y/b);

#################################################################
#				Mesh creation & DOF (ID) array					#
#################################################################
# Create a rectangular mesh, based on the discretization lengths provided above.
println("Creating initial mesh.")
p_in,e_in,el_in = Mesh2Drect(x,y);

# Create the additional nodes for shape functions with order r>1.
println("Creating additional DOFs.")
p,el = intNodes(el_in,p_in,r,hx,hy,x,y);
Nel = size(el,1); # Number of elements.



# Converting element indices to Ints, might be possibly fixed in next Julia version.
el = floor.(Int,el);

# Create the ID array for the degrees of freedom.
DOF = ID(p);
DOF = floor.(Int,DOF);
findind = zeros(size(DOF));
# Find the correspondence between the DOF indices and node indices.
findind .= DOF;
for i = 1 : length(findind)
	tempind = [findind[1:i-1]; NaN; findind[i+1:end]];
	ind3 = findind[i] .== tempind;
	if maximum(ind3) != false
		ind4 = findall(x -> x!=0,ind3);
		findind[ind4[end]] = 0;
	end
end

findind2 = findall(x -> x != 0, vec(findind));

#Open Nodes
ON = p[findind2,:];
#NdofL = length(vec(findall(x -> x - xI > -1e6*hx, vec(ON[:,1]))));
NdofL = 0;
for i = 1 : size(ON,1)
	global NdofL
	if ON[i,1] >= xI
		NdofL = NdofL + 1;
	end
end
Ndof = Int(maximum(DOF));



#################################################################
#				Saving DOF and solution at DOF					#
#################################################################
#
if SaveSol == 1
f = open("points_r$r.dat","w")
pf = p[findind2,:];
for i = 1 : size(pf,1)
    for j = 1 : 3
        print(f,pf[i,j]," ")
    end
    println(f,"")
end

close(f)

#function end
#end
