using SparseArrays
using LinearAlgebra
using Printf
using NonlinearEigenproblems
using Arpack
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
hy = hx;
#Waveguide characteristics.
b  = 3;  # Width.
xI = 3;  # Location of interface.
nL = 2;
xE = 4;#round(xI + nL*hx,digits=5); # Location of ABC.

c  = sqrt(3);  # Wave speed.
s  = 0;
P  = 3;  # Order of ABC.
#Problem setup
r  = 1;
osc_f =1;
x  = Array(0 : hx : xE);    
y  = Array(0 : hy : b );
Nx = length(x);
Ny = length(y);
N  = Nx*Ny;
# Time discretization setup
#dt = min(hx,hy)/(100*r);          % Time step
T = 500;
dt = hx/(100*r);
t = Array(0 : dt :T);      # Time vector

SaveSol = 0; # Option if we wish to save the solution in each time step to a file.

# Gaussian pulse as an initial condition. Can use as many sources as we want,
# if we supply their locations in u_0 (see below).
radius = 0.5;   				# Radius of the Gaussian pulse.
a = log(10^6)/(radius*radius);  # Amplitude.
xstar =  1.5;						# Location in x-direction.
ystar =  1.5;						# Location in y-direction.
xstar2 = 8;						# Location in x-direction for 2nd source.
ystar2 = 7.5;					# Location in y-direction for 2nd source

u0(x,y,xstar,ystar)=hx*sqrt(2*a/pi)*exp.(-a*((x .-xstar).^2+(y .-ystar).^2));



# A more complex initial condition.
#u0(x,y) = ((x .- 8.5).^2 .- 1).^2 .*sin.(2*pi*y/b).*(x .>= 7.5).*(x .<= 9.5);

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

Ndof = Int(maximum(DOF));

# Initial condition.
# Initial condition for udot.
upI = zeros(size(p[findind2,1],1),1);
#u = ex_sol(osc_f,0,c,p[findind2,1],p[findind2,2]);
uI  = u0(p[findind2,1],p[findind2,2],xstar,ystar);
#uI  = u0(p[findind2,1],p[findind2,2]);



#################################################################
#						Matrix Assembly							#
#################################################################

# Create Gauss-Lobatto weights for the reference element points.
gj = intWeights(r,2,2);
# Element Stiffness and Mass matrices for rectangular elements.
I1e = MMtrx(r,hx,hy,gj);
I2e = KMtrx(r,hx,hy,gj);
I3e = I3(r,hy,gj);
I4e = I4(r,hy,gj);
I5e = I5(r,hy,gj);
I6e = I6(r,hy,gj);
# Assemble the stiffness and mass matrices for all the D.O.F.`s.
println("Assembling matrices.")
M  = spzeros(Ndof,Ndof);
MJ = spzeros(Ndof,Ndof);
C0 = spzeros(Ndof,Ndof);
K  = spzeros(Ndof,Ndof);
KJ = spzeros(Ndof,Ndof);
G  = spzeros(Ndof,Ndof);
H  = spzeros(Ndof,Ndof);
#
CJ = spzeros(Ndof,Ndof);
A  = spzeros(Ndof,Ndof);
B  = spzeros(Ndof,Ndof);
#
CP = spzeros(Ndof,Ndof);

for l = 1 : Nel
	for i = 1 : (r+1)^2
		m = Int(DOF[Int(el[l,i])]);
		if m!=0
			for j = 1 : (r+1)^2
				n = Int(DOF[Int(el[l,j])]);
				if n!=0
					M[m,n]  = M[m,n] +     I1e[i,j];
					K[m,n]  = K[m,n] + c^2*I2e[i,j] + s*I1e[i,j];

					if ((xI - p[el[l,1],1]) < 1e-6*hx)
						MJ[m,n] = MJ[m,n] +     I1e[i,j];
						KJ[m,n] = KJ[m,n] + c^2*I2e[i,j] + s*I1e[i,j];
					end
					if (abs(p[el[l,1],1]-xI)<1e-6*hx)
						 A[m,n] =  A[m,n] - c  *I4e[i,j];
						 B[m,n] =  B[m,n] - c^2*I6e[i,j];
						CP[m,n] = CP[m,n] + c  *I4e[i,j];
						CJ[m,n] = CJ[m,n] + c  *I4e[i,j];
					end

					if (abs(p[el[l,2],1]-xE)<1e-6*hx)
						 G[m,n] =  G[m,n] - c  *I3e[i,j];
						 H[m,n] =  H[m,n] + c^2*I5e[i,j];
						C0[m,n] = C0[m,n] + c  *I3e[i,j];
						CJ[m,n] = CJ[m,n] + c  *I3e[i,j];
					end
				end
			end
		end
	end
end 

if P == 0
	Mbig = spzeros(Ndof*(P+2),Ndof*(P+2));
	Kbig = spzeros(Ndof*(P+2),Ndof*(P+2));
	Cbig = spzeros(Ndof*(P+2),Ndof*(P+2));
else
	Mbig = spzeros(Ndof*(P+1),Ndof*(P+1));
	Kbig = spzeros(Ndof*(P+1),Ndof*(P+1));
	Cbig = spzeros(Ndof*(P+1),Ndof*(P+1));
end



Kbig[1:Ndof,1:Ndof] = K;
Kbig[1:Ndof,Ndof+1:2*Ndof] = H;
#
Cbig[1:Ndof,1:Ndof] = C0;
Cbig[1:Ndof,Ndof+1:2*Ndof] = G;
#
Mbig[1:Ndof,1:Ndof] = M;

for i = 1 : P - 1 
	Kbig[Ndof*(i  )+1 : Ndof*(i+1),Ndof*(i-1)+1 : Ndof*(i  )] = B;
	Kbig[Ndof*(i  )+1 : Ndof*(i+1),Ndof*(i  )+1 : Ndof*(i+1)] = KJ;
	Kbig[Ndof*(i  )+1 : Ndof*(i+1),Ndof*(i+1)+1 : Ndof*(i+2)] = H;
	

	Cbig[Ndof*(i  )+1 : Ndof*(i+1),Ndof*(i-1)+1 : Ndof*(i  )] = A;
	Cbig[Ndof*(i  )+1 : Ndof*(i+1),Ndof*(i  )+1 : Ndof*(i+1)] = CJ;
	Cbig[Ndof*(i  )+1 : Ndof*(i+1),Ndof*(i+1)+1 : Ndof*(i+2)] = G;
	
	Mbig[Ndof*(i  )+1 : Ndof*(i+1),Ndof*(i  )+1 : Ndof*(i+1)] = MJ;
		
end



if(P>=1)
	Kbig[Ndof*(P  )+1 : Ndof*(P+1),Ndof*(P-1)+1 : Ndof*(P  )] = B;
	Kbig[Ndof*(P  )+1 : Ndof*(P+1),Ndof*(P  )+1 : Ndof*(P+1)] = KJ;

	Cbig[Ndof*(P  )+1 : Ndof*(P+1),Ndof*(P-1)+1 : Ndof*(P  )] = A;
	Cbig[Ndof*(P  )+1 : Ndof*(P+1),Ndof*(P  )+1 : Ndof*(P+1)] = CJ;

	Mbig[Ndof*(P  )+1 : Ndof*(P+1),Ndof*(P  )+1 : Ndof*(P+1)] = MJ;

end


i = 1 ;
while i <= size(Mbig,1)
	global Mbig,Kbig,Cbig,i
	if Mbig[i,i] == 0;
	Mbig = [Mbig[1:i-1,1:i-1] Mbig[1:i-1,i+1:end] ; Mbig[i+1:end,1:i-1] Mbig[i+1:end,i+1:end]];
	Kbig = [Kbig[1:i-1,1:i-1] Kbig[1:i-1,i+1:end] ; Kbig[i+1:end,1:i-1] Kbig[i+1:end,i+1:end]];
	Cbig = [Cbig[1:i-1,1:i-1] Cbig[1:i-1,i+1:end] ; Cbig[i+1:end,1:i-1] Cbig[i+1:end,i+1:end]];
	else
		i = i + 1;
	end
end




println("Solving QEP")
#skb = size(Kbig,1);
#A = spzeros(2*skb,2*skb);
#B = spzeros(2*skb,2*skb);
#A[1:skb,skb+1:2:skb] = sparse(Matrix(I,skb,skb));
#A[skb+1:2*skb,1:skb] = Cbig;
#A[skb+1:2*skb,skb+1:2*skb] = Kbig;
#A[skb+1:2*skb,1:skb] = Kbig;

#B[1:skb,1:skb] = sparse(Matrix(I,skb,skb));
#B[skb+1:2*skb,skb+1:2*skb] = Mbig;

#lambda,vevs = eigs(A,B);
#println(lambda)
#
#
qep = PEP([Kbig,Cbig,Mbig]);
#println("QEP created")
#unit_square = float([1+1im, 1-1im, -1-1im,-1+1im]);
#lambda,eigvecs = nleigs(qep,unit_square)


E,A = companion(qep);
#println("Companion Linearization done")
#lambda,vevs = eigs(A,E, nev=10);
#println(maximum(real(lambda)))
lambda, = eigs(A,E,nev=size(Kbig,1));
println(maximum(real(lambda)))
println("Eigs computed")
#println(lambda)
#lambda,vecs = eigs(A,B,v0 = ones(size(A,1)));
#println(maximum(real(lambda)))
#println(maximum(imag(lambda)))

