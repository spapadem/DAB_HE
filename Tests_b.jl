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
hx = 0.1;
hy = hx;
#Waveguide characteristics.
b  = 3;  # Width.
xI = 3;  # Location of interface.
nL = 10;
xE = round(xI + nL*hx,digits=5); # Location of ABC.

c  = sqrt(3);  # Wave speed.
s  = 0;
P  = 4;  # Order of ABC.
#Problem setup
r = 2;
osc_f =1;
x  = Array(0 : hx : xE);    
y  = Array(0 : hy : b );
Nx = length(x);
Ny = length(y);
N  = Nx*Ny;
# Time discretization setup
#dt = min(hx,hy)/(100*r);          % Time step
T = 50;
dt = hx/(100*r);
t = Array(0 : dt :T);      # Time vector
@printf("| r: %3.2f | P: %3.3g \n",r,P);

SaveSol = 0; # Option if we wish to save the solution in each time step to a file.

# Gaussian pulse as an initial condition. Can use as many sources as we want,
# if we supply their locations in u_0 (see below).
radius = 0.5;   				# Radius of the Gaussian pulse.
a = log(10^6)/(radius*radius);  # Amplitude.
xstar =  1.5;						# Location in x-direction.
ystar =  1.5;						# Location in y-direction.
xstar2 = 8;						# Location in x-direction for 2nd source.
ystar2 = 7.5;					# Location in y-direction for 2nd source

u0(x,y,xstar,ystar)=hx/r*sqrt(2*a/pi)*exp.(-a*((x .-xstar).^2+(y .-ystar).^2));



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
	Cbig[Ndof*(P  )+1 : Ndof*(P+1),Ndof*(P  )+1 : Ndof*(P+1)] = CP;

	Mbig[Ndof*(P  )+1 : Ndof*(P+1),Ndof*(P  )+1 : Ndof*(P+1)] = MJ;

end


# End of assemble.
# Calculating the inverse of M. It is lumped due to spectral elements, 
# so we simply invert the diagonal elements.
Minv = spzeros(size(Mbig,1),size(Mbig,2));
for i = 1 : size(Minv,1)
	if Mbig[i,i] != 0
		Minv[i,i] = 1/Mbig[i,i];
#	else
#		Minv[i,i] = 1;
	end
end

rel   = zeros(size(t));
uel   = zeros((r+1)^2,1);
uddel = zeros((r+1)^2,1);
fel1  = zeros((r+1)^2,1);
fel2  = zeros((r+1)^2,1);
fel4  = zeros((r+1)^2,1);

#################################################################
#						Time-Stepping							#
#################################################################

# Time-stepping with 4th-order Runge-Kutta for the problem \ddot{u} = g(t,\dot{u},u).
println("Starting time-stepping.")
ufull = zeros(size(p,1),1);

u = spzeros(size(Mbig,2));
up = spzeros(size(Mbig,2));

u[1:Ndof]  = uI;
up[1:Ndof] = upI;

if SaveSol == 1
	filename = "solt_r$r-P$P.dat";
	ft = open(filename,"w")
end

for i = 1 :length(t)
	global u
	global up
#=
#	Variables for the values of F(t) for the intermediate steps for the RK4 method.
#	We do not create F3 since the times t2 & t3 are the same.
	F1 = zeros(Int(Ndof),1);
	F2 = zeros(Int(Ndof),1);
	F4 = zeros(Int(Ndof),1);
#	Assembly of right hand side.
	for l = 1 : Nel
		for m = 1 : (r+1)^2
			ind1 = Int(DOF[Int(el[l,m])]);
			if ind1!=0
				fel1[m]   = 0; #f_ex(osc_f,t[i]     ,c,p[el[l,m],1],p[el[l,m],2]);
				fel2[m]   = 0;#f_ex(osc_f,t[i]+dt/2,c,p[el[l,m],1],p[el[l,m],2]);
				fel4[m]   = 0;#f_ex(osc_f,t[i]+dt  ,c,p[el[l,m],1],p[el[l,m],2]);
			end
		end

		Fel1 =  Mel*fel1;
		Fel2 =  Mel*fel2;
		Fel4 =  Mel*fel4;

		for k = 1 : (r+1)^2
			ind1 = Int(DOF[el[l,k]]);
			if ind1 != 0
				F1[ind1] = F1[ind1] + Fel1[k];
				F2[ind1] = F2[ind1] + Fel2[k];
				F4[ind1] = F4[ind1] + Fel4[k];
			end
		end
	end


=#
#################################################################
#			Printing Solution to file, if  specified			#
#################################################################
if SaveSol == 1 && mod(i,100) == 0
#	ufull[findall(vec(DOF))] = u;
	for tt = 1 : size(uI,1)
		if p[findind2[tt],1] <= xI
	    print(ft,u[tt]," ")
		end
	end
	println(ft,"")
end
	g(ddot,d) = Minv*(-Cbig*ddot-Kbig*d); 
#	Intermediate time steps.
    y1 = u ;
    x1 = up;
    y2 = u  + .5*dt*x1;
    x2 = up + .5*dt*g(x1,y1);
    y3 = u  + .5*dt*x2;
    x3 = up + .5*dt*g(x2,y2);
    y4 = u  +    dt*x3;
    x4 = up +    dt*g(x3,y3);
    up = up + dt/6*(g(x1,y1) + 2*g(x2,y2) + 2*g(x3,y3) + g(x4,y4));
    u  = u  + dt/6*(x1       + 2*x2       + 2*x3       + x4);

#	errt,nu,nuh = femerror(p,el,ufull,hx,hy,r,osc_f,t[i],c);
#	rel[i] = errt/nuh;
	if mod(10*t[i],1) == 0
#		@printf("| Time: %3.2f |\n",t[i]);
		@printf("| Time: %3.2f |Norm Sol. : %3.3g \n",t[i],norm(u));
#		@printf("| Time: %3.2f | Norm Sol.: %3.3g | Norm. App.: %3.3g | Rel. err. %3.3e\n",t[i],nu,nuh,rel[i]);
#		println("Max PhiP: ",maximum(abs.(u[end-Ndof:end])))
#		println("Max Phi0: ",maximum(abs.(u[1:Ndof])))
	end


end

println("Time-stepping finished.")

@printf("| Time: %3.2f |Norm Sol. : %3.3g \n",t[end],norm(u));

if SaveSol == 1
	close(ft)
end


#################################################################
#			Calculating the error at the final time				#
#################################################################

#ufull[findall(vec(DOF))] = u;
#errt,nu,nuh = femerror(p,el,ufull,hx,hy,r,osc_f,t[end]+dt,c);
#rel_err = errt/nu;


#################################################################
#				Saving DOF and solution at DOF					#
#################################################################
f = open("points.dat","w")
pf = p[findind2,:];
for i = 1 : size(pf,1)
    for j = 1 : 3
        print(f,pf[i,j]," ")
    end
    println(f,"")
end

close(f)


f = open("sol.dat","w")


for tt = 1 : size(uI,1)
        print(f,u[tt]," ")
    end
close(f)



rel_err = 0;
return  rel_err

#function end
#end
