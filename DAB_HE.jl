using SparseArrays
using LinearAlgebra
using Printf

include("basis.jl")
include("basisx.jl")
include("bubble_sort.jl")
include("ChkNodeEx.jl")
include("ex_sol.jl")
include("f_ex.jl")
include("femerror.jl")
include("Fvec.jl")
include("GLpw.jl")
include("green.jl")
include("I1.jl")
include("I2.jl")
include("I3.jl")
include("I4.jl")
include("I5.jl")
include("I6.jl")
include("ID.jl")
include("intNodes.jl")
include("intWeights.jl")
include("Mesh2Drect.jl")



#################################################################
#						Problem Setup							#
#################################################################
f = 24; # Frequency
f_0 = 75; # Reference frequency
c_0 = 1500; # Speed of sound
s  = 0; # Damping parameter
lambda_0 = c_0/f_0; # Wavelength at reference frequency
lambda   = c_0/f; # Wavelength at frequency f
kn = 2*pi*f/c_0; # Wave number

# Create points for mesh generation
hx = lambda_0/5; # Discretization length in x-direction
hy = hx; # Discretization length in y-direction

#Waveguide characteristics.
b  = 10*lambda_0;  # Width.
xI = 15*lambda_0;  # Location of interface.
xE = 17*lambda_0;  # Location of ABC.
nL = (xE-xI)/hx; # Number of discretization points in the ABC
P  = 10;  		   # Order of ABC.

#Problem setup
r = 3; # Order of the basis functions
x  = Array(12*lambda_0 : hx : xE); # x-coordinates of the mesh points
y  = Array(0  : hy : b ); # y-coordinates of the mesh points
Nx = length(x); # Number of mesh points in x-direction
Ny = length(y); # Number of mesh points in y-direction
N  = Nx*Ny; # Total number of mesh points


SaveSol = 1; # Option if we wish to save the solution to a file.

#################################################################
#				Mesh creation & DOF (ID) array					#
#################################################################
# Create a rectangular mesh, based on the discretization lengths provided above.
println("Creating initial mesh.")
p_in,e_in,el_in = Mesh2Drect(x,y);

# Create the additional nodes for shape functions with order r>1.
println("Creating additional nodes.")
p,el = intNodes(el_in,p_in,r,hx,hy,x,y);
Nel  = size(el,1); # Number of elements.
	
println("Creating DOFs.")
# Converting element indices to Ints, might be possibly fixed in next Julia version.
el = floor.(Int,el);

# Create the ID array for the degrees of freedom.
DOF = ID(p);
DOF = floor.(Int,DOF);
findind = zeros(size(DOF));

# Find the correspondence between the DOF indices and node indices.
findind .= DOF;
for i = 1 : length(findind)
	ind3 = findind[i] .== findind;
	if (maximum(ind3) != false && ind3 != ind3)
		ind4 = findall(x -> x!=0,ind3);
		findind[ind4[end]] = 0;
	end
end

# Find the indices of the DOFs that are open.
findind2 = findall(x -> x != 0, vec(findind));
Ndof = Int(maximum(DOF)); # Number of DOFs.

#################################################################
#						Matrix Assembly							#
#################################################################
# Create the reference element and its integration points.
pref = [ -1 -1 1;
			1 -1 1;
			1  1 1;
			-1  1 1];
elref = [1 2 3 4];

pref,elref = intNodes(elref,pref,r,2,2,[-1 1],[-1 1]);
elref = elref[1,:];



println("Creating element matrices.")
# Create Gauss-Lobatto weights for the reference element points.
gj = intWeights(r,2,2);
# Element Stiffness and Mass matrices for rectangular elements.
I1e = I1(r,hx,hy,gj,pref,elref);
I2e = I2(r,hx,hy,gj,pref,elref);
I3e = I3(r,hy,gj,pref,elref);
I4e = I4(r,hy,gj,pref,elref);
I5e = I5(r,hx,hy,gj,pref,elref);
I6e = I6(r,hx,hy,gj,pref,elref);
# Assemble the stiffness and mass matrices for all the D.O.F.`s.
println("Assembling matrices.")
M  = spzeros(Complex{Float64},Ndof,Ndof);
MJ = spzeros(Complex{Float64},Ndof,Ndof);
K  = spzeros(Complex{Float64},Ndof,Ndof);
KJ = spzeros(Complex{Float64},Ndof,Ndof);
G  = spzeros(Complex{Float64},Ndof,Ndof);
H  = spzeros(Complex{Float64},Ndof,Ndof);
F  = spzeros(Complex{Float64},Ndof,1);
#

for l = 1 : Nel
	for i = 1 : (r+1)^2
		m = Int(DOF[Int(el[l,i])]);
		if m!=0
			for j = 1 : (r+1)^2
				n = Int(DOF[Int(el[l,j])]);
				if n!=0
					K[m,n] = K[m,n] + I2e[i,j] - kn^2*I1e[i,j];
					M[m,n] = M[m,n] + I1e[i,j];

					if ((xI - p[el[l,1],1]) < 1e-6*hx)
						KJ[m,n] = KJ[m,n] + I2e[i,j] - kn^2*I1e[i,j];
						MJ[m,n] = MJ[m,n] + I1e[i,j];
					end
					if (abs(p[el[l,1],1]-xI)<1e-6*hx)
							KJ[m,n] = KJ[m,n] - 1im*kn*I4e[i,j];
							G[m,n] =  G[m,n] + 1im*kn*I4e[i,j] - I6e[i,j];
					end

					if (abs(p[el[l,2],1]-xE)<1e-6*hx)
							K[m,n] =  K[m,n] - 1im*kn*I3e[i,j];
							KJ[m,n] = KJ[m,n] - 1im*kn*I3e[i,j];
							H[m,n] =  H[m,n] + 1im*kn*I3e[i,j] + I5e[i,j];
					end
				end
			end
		end
	end
end

A = 1im;
ky  = pi/b;
kx  = sqrt(Complex(kn*kn-ky*ky));

#   Assembly of right hand side.
for l = 1 : Nel
	fel = zeros(Complex,(r+1)^2,1);
	Fel = zeros(Complex,(r+1)^2,1);
	for m = 1 : (r+1)^2
		ind1 = Int(DOF[Int(el[l,m])]);
		if ind1==0
			# Green's function for the Helmholtz equation on the left side of the domain.
			fel[m] = green(xstar,ystar,p[el[l,m],1],p[el[l,m],2],kn,b,floor(2*b/lambda));
		end
	end

	Fel = -(I2e-kn^2*I1e-1im*kn*I3e)*fel;
	for k = 1 : (r+1)^2
		ind1 = Int(DOF[el[l,k]]);
		if ind1 != 0
			F[ind1] = F[ind1] + Fel[k];
		end
	end
end

if P == 0
	Mbig = spzeros(Complex{Float64},Ndof*(P+2),Ndof*(P+2));
	Kbig = spzeros(Complex{Float64},Ndof*(P+2),Ndof*(P+2));
	Fbig =   zeros(Complex{Float64},Ndof*(P+2),1);
else
	Mbig = spzeros(Complex{Float64},Ndof*(P+1),Ndof*(P+1));
	Kbig = spzeros(Complex{Float64},Ndof*(P+1),Ndof*(P+1));
	Fbig =   zeros(Complex{Float64},Ndof*(P+1),1);
end

Kbig[1:Ndof,1:Ndof] = K;
Mbig[1:Ndof,1:Ndof] = M;
Kbig[1:Ndof,Ndof+1:2*Ndof] = H;
Fbig[1:Ndof] = F;
	
for i = 1 : P - 1
	Mbig[Ndof*i+1 : Ndof*(i+1),Ndof*(i  )+1 : Ndof*(i+1)] = MJ;

	Kbig[Ndof*i+1 : Ndof*(i+1),Ndof*(i-1)+1 : Ndof*(i  )] = G;
	Kbig[Ndof*i+1 : Ndof*(i+1),Ndof*(i  )+1 : Ndof*(i+1)] = KJ;
	Kbig[Ndof*i+1 : Ndof*(i+1),Ndof*(i+1)+1 : Ndof*(i+2)] = H;
end
	
if(P>=1)
	Mbig[Ndof*P+1 : Ndof*(P+1),Ndof*(P  )+1 : Ndof*(P+1)] = MJ;

	Kbig[Ndof*P+1 : Ndof*(P+1),Ndof*(P-1)+1 : Ndof*(P  )] = G;
	Kbig[Ndof*P+1 : Ndof*(P+1),Ndof*(P  )+1 : Ndof*(P+1)] = KJ;
end


println("Resizing.")
reindx = findall(x -> x != 0, diag(Mbig)); # Choosing the indices for which the mass matrix is non-zero on its diagonal.
Kbig = Kbig[reindx,reindx]; # Resizing the stiffness matrix.
Fbig = Fbig[reindx]; # Resizing the right hand side.

println("Solving Linear Problem.")
u = Kbig\Fbig;


# Saving solution begins.
if SaveSol == 1
	println("Saving solution.")
	fp = open("points_r$r-h$hx.dat","w")
	pf = p[findind2,:];
	for i = 1 : size(pf,1)
		if p[findind2[i],1] <= xI
			for j = 1 : 3
				print(fp,pf[i,j]," ")
			end

			println(fp,"")
		end
	end

	close(fp)

	g = green(xstar,ystar,p[findind2,1],p[findind2,2],kn,b,floor(2*b/lambda));
	
	fr = open("sol_f$f-r$r-P$P-nL$nL-h$hx-real.dat","w");
	fi = open("sol_f$f-r$r-P$P-nL$nL-h$hx-imag.dat","w");

	fre = open("exact_f$f-r$r-h$hx-real.dat","w");
	fie = open("exact_f$f-r$r-h$hx-imag.dat","w");

	for tt = 1 : size(pf,1)
		if p[findind2[tt],1] <= xI
			print(fr ,real(u[tt])," ")
			print(fi ,imag(u[tt])," ")
			print(fre,real(g[tt])," ")
			print(fie,imag(g[tt])," ")
		end
	end

	close(fr)
	close(fi)
	close(fre)
	close(fie)

end
# Saving solution end.