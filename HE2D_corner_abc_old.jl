using SparseArrays
using LinearAlgebra
using Printf

#function  W2DFE_Dbc_RK4(r,hx,hy)

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
include("I7.jl")
include("I8.jl")
include("I9.jl")
include("I10.jl")
include("green_term.jl")

#################################################################
#						Problem Setup							#
#################################################################
f = 24;
f_0 = 75;
c_0 = 1500;
s  = 0;
lambda_0 = c_0/f_0;
lambda   = c_0/f;
kn = 2*pi*f/c_0;
# Create points for mesh generation
hx = lambda_0/10;
hy = hx;
#Waveguide characteristics.
xIN = 3*lambda_0;  # Width.
xN  = 5*lambda_0;  # Width.
xIE = 3*lambda_0;  # Location of interface.
xE  = 5*lambda_0;  # Location of ABC.
nLE = (xE-xIE)/hx;
nLN = (xN-xIN)/hy;
P = 7;
PE  = P;  		   # Order of ABC in right layer.
PN  = P;  		   # Order of ABC in top layer.
#Problem setup
r = 1;
x  = Array(0 : hx : xE);    
y  = Array(0 : hy : xN);
Nx = length(x);
Ny = length(y);
N  = Nx*Ny;
# Time discretization setup
#dt = min(hx,hy)/(100*r);          % Time step
SaveSol = 1; # Option if we wish to save the solution in each time step to a file.

# Gaussian pulse as an initial condition. Can use as many sources as we want,
# if we supply their locations in u_0 (see below).
radius = 5;   					# Radius of the Gaussian pulse.
a = log(10^6)/(radius*radius);  # Amplitude.
xstar =  80;					# Location in x-direction.
ystar =  52;					# Location in y-direction.
xstar2 = 8;						# Location in x-direction for 2nd source.
ystar2 = 7.5;					# Location in y-direction for 2nd source

u0(x,y,xstar,ystar)= sqrt(a/pi)*exp.(-a*((x .- xstar).^2 + (y .- ystar).^2)) .*
                   					((  ((x .- xstar).^2 + (y .- ystar).^2)) .<= radius); 

xstar = -20;
ystar =  20;

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


findind2 = findall(x -> x != 0, vec(findind));

#Open Nodes
#ON = p[findind2,:];
#NdofL = length(vec(findall(x -> x - xI > -1e6*hx, vec(ON[:,1]))));
#NdofL = 0;
#for i = 1 : size(ON,1)
#	global NdofL
#	if ON[i,1] >= xI
#		NdofL = NdofL + 1;
#	end
#end
Ndof = Int(maximum(DOF));

#################################################################
#						Matrix Assembly							#
#################################################################

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
I1e  =  I1(r,hx,hy,gj,pref,elref);
I2e  =  I2(r,hx,hy,gj,pref,elref);
I3e  =  I3(r,hy,gj,pref,elref);
I4e  =  I4(r,hy,gj,pref,elref);
I5e  =  I5(r,hx,hy,gj,pref,elref);
I6e  =  I6(r,hx,hy,gj,pref,elref);
I7e  =  I7(r,hy,gj,pref,elref);
I8e  =  I8(r,hy,gj,pref,elref);
I9e  =  I9(r,hx,hy,gj,pref,elref);
I10e = I10(r,hx,hy,gj,pref,elref);
# Assemble the stiffness and mass matrices for all the D.O.F.`s.
println("Assembling matrices.")
# Main domain.
W  = spzeros(Complex{Float64},Ndof,Ndof);
X  = spzeros(Complex{Float64},Ndof,Ndof);
Y  = spzeros(Complex{Float64},Ndof,Ndof);
# Right Layer.
E  = spzeros(Complex{Float64},Ndof,Ndof);
D  = spzeros(Complex{Float64},Ndof,Ndof);
G  = spzeros(Complex{Float64},Ndof,Ndof);
# North Layer.
N  = spzeros(Complex{Float64},Ndof,Ndof);
M  = spzeros(Complex{Float64},Ndof,Ndof);
O  = spzeros(Complex{Float64},Ndof,Ndof);
# Corner.
C  = spzeros(Complex{Float64},Ndof,Ndof);
CD = spzeros(Complex{Float64},Ndof,Ndof);
CG = spzeros(Complex{Float64},Ndof,Ndof);
CM = spzeros(Complex{Float64},Ndof,Ndof);
CO = spzeros(Complex{Float64},Ndof,Ndof);

IMW = spzeros(Complex{Float64},Ndof,Ndof);
IME = spzeros(Complex{Float64},Ndof,Ndof);
IMN = spzeros(Complex{Float64},Ndof,Ndof);
IMC = spzeros(Complex{Float64},Ndof,Ndof);

# Right-hand side vector.
F  = spzeros(Complex{Float64},Ndof,1);
#

for l = 1 : Nel
	for i = 1 : (r+1)^2
		m = Int(DOF[Int(el[l,i])]);
		if m!=0
			for j = 1 : (r+1)^2
				n = Int(DOF[Int(el[l,j])]);
				if n!=0
					W[m,n] = W[m,n] + I2e[i,j] - kn*kn*I1e[i,j];
					IMW[m,n] = IMW[m,n] + I1e[i,j];
					#-----------------Right layer. ---------------------#
					# Inside the domain.								#
					if ((xIE - p[el[l,1],1]) < 1e-6*hx)					#
						E[m,n] = E[m,n] + I2e[i,j] - kn*kn*I1e[i,j];	#
						IME[m,n] = IME[m,n] + I1e[i,j];					#
					end													#
					# Right interface boundary.							#
					if (abs(p[el[l,1],1]-xIE)<1e-6*hx)					#
						D[m,n] = D[m,n] + 1im*kn*I4e[i,j] - I6e[i,j];	#
						E[m,n] = E[m,n] - 1im*kn*I4e[i,j];				#
					end													#
					# Right boundary.									#
					if (abs(p[el[l,2],1]-xE)<1e-6*hx)					#
						W[m,n] =  W[m,n] - 1im*kn*I3e[i,j];				#
						E[m,n] =  E[m,n] - 1im*kn*I3e[i,j];				#
						X[m,n] =  X[m,n] + 1im*kn*I3e[i,j] + I5e[i,j];	#
						G[m,n] =  G[m,n] + 1im*kn*I3e[i,j] + I5e[i,j];	#
					end													#
					#---------------------------------------------------#

					#--------------------Top layer. --------------------#
					# Inside the domain.
					if ((xIN - p[el[l,1],2]) < 1e-6*hy)					#
						N[m,n] = N[m,n] + I2e[i,j] - kn*kn*I1e[i,j];	#
						IMN[m,n] = IMN[m,n] + I1e[i,j];					#
					end													#
					# Top interface boundary.							#
					if (abs(p[el[l,1],2]-xIN)<1e-6*hy)					#
						M[m,n] = M[m,n] + 1im*kn*I8e[i,j] - I10e[i,j];	#
						N[m,n] = N[m,n] - 1im*kn*I8e[i,j];				#
					end													#
					# Top boundary.										#
					if (abs(p[el[l,3],2]-xN)<1e-6*hy)					#
						W[m,n] =  W[m,n] - 1im*kn*I7e[i,j];				#
						N[m,n] =  N[m,n] - 1im*kn*I7e[i,j];				#
						Y[m,n] =  Y[m,n] + 1im*kn*I7e[i,j] + I9e[i,j];	#
						O[m,n] =  O[m,n] + 1im*kn*I7e[i,j] + I9e[i,j];	#
					end													#
					#---------------------------------------------------#
					
					#-----------------Corner layer.---------------------#
					# Inside the domain.								#
					if (((xIE - p[el[l,1],1]) < 1e-6*hx) &&				#
                        ((xIN - p[el[l,1],2]) < 1e-6*hy))			    #
						  C[m,n] =   C[m,n] + I2e[i,j] - kn*kn*I1e[i,j];#
						IMC[m,n] = IMC[m,n] + I1e[i,j];			        #
					end													#
					# Right interface boundary.							#
					if ((abs(p[el[l,1],1] - xIE) < 1e-6*hx) &&		    #
                        ((xIN - p[el[l,1],2])    < 1e-6*hy))			#
						 C[m,n] =  C[m,n] - 1im*kn*I4e[i,j];			#
						CD[m,n] = CD[m,n] + 1im*kn*I4e[i,j] - I6e[i,j];	#
					end													#
					# Right boundary.									#
					if ((abs(p[el[l,2],1] - xE) < 1e-6*hx) &&			#
                        ((xIN - p[el[l,1],2])   < 1e-6*hy))				#
						 N[m,n] =  N[m,n] - 1im*kn*I3e[i,j];	        #
						 C[m,n] =  C[m,n] - 1im*kn*I3e[i,j];			#
						CG[m,n] = CG[m,n] + 1im*kn*I3e[i,j] + I5e[i,j];	#
					end													#
					# Top interface boundary.							#
					if ((abs(p[el[l,1],2]-xIN) < 1e-6*hy) &&			#
                        ((xIE - p[el[l,1],1])   < 1e-6*hx))				#
						 C[m,n]  = C[m,n] - 1im*kn*I8e[i,j];			#
						CM[m,n] = CM[m,n] + 1im*kn*I8e[i,j] - I10e[i,j]	#
					end													#
					# Top boundary.										#
					if ((abs(p[el[l,3],2]-xN) < 1e-6*hy) &&				#
                        ((xIE - p[el[l,1],1])  < 1e-6*hx))				#
   				 	     E[m,n] =  E[m,n] - 1im*kn*I7e[i,j];			#
						 C[m,n] =  C[m,n] - 1im*kn*I7e[i,j];			#
						CO[m,n] = CO[m,n] + 1im*kn*I7e[i,j] + I9e[i,j];	#
					end													#
					#---------------------------------------------------#
				end
			end
		end
	end
end 

#ky  = pi/b;
#kx  = sqrt(Complex(kn*kn-ky*ky));
#   Assembly of right hand side.
for l = 1 : Nel
fel = zeros(Complex,(r+1)^2,1);
Fel = zeros(Complex,(r+1)^2,1);
	for m = 1 : (r+1)^2
		ind1 = Int(DOF[Int(el[l,m])]);
		if ind1==0
#			Gaussian			
#			fel[m] = u0(p[el[l,m],1],p[el[l,m],2],xstar,ystar);
#			Dirac
#			fel[m] = r*Float64((p[el[l,m],1] .== xstar) .* (p[el[l,m],2] .== ystar));
#			Non homogeneous Dirichlet
			dist = sqrt((p[el[l,m],1]-xstar)^2+(p[el[l,m],2]-ystar)^2);
            fel[m] = exp.(1im*kn*dist)/(4*pi*dist);
#            fel[m] = u0(p[el[l,m],1],p[el[l,m],2],xstar,ystar);
#			fel[m] = cos.(ky*p[el[l,m],2]);
		end
	end

#	Fel = -(I2e-kn^2*I1e-1im*kn*I3e-1im*kn*I7e)*fel;
   Fel = -(I2e-kn^2*I1e)*fel;
#	Fel = I1e*fel;
	for k = 1 : (r+1)^2
		ind1 = Int(DOF[el[l,k]]);
		if ind1 != 0
			F[ind1] = F[ind1] + Fel[k];
		end
	end
end

#findx = findall(x->x!=0, (p[findind2,1] .== xstar .* p[findind2,2] .== ystar));

#F[findall(x->x!=0, F)[1][1]] = 1


# Assembling the big matrix.
println("Assembling big matrix.")
Mbig = spzeros(Complex{Float64},Ndof*(PE+PN+1+PE*PN),Ndof*(PE+PN+1+PE*PN));
Kbig = spzeros(Complex{Float64},Ndof*(PE+PN+1+PE*PN),Ndof*(PE+PN+1+PE*PN));
Fbig =   zeros(Complex{Float64},Ndof*(PE+PN+1+PE*PN),1);


# Full domain.
Kbig[1:Ndof,              1 : Ndof       ] = W;
Kbig[1:Ndof,Ndof        + 1 : Ndof*2     ] = X;
Kbig[1:Ndof,Ndof*(PE+1) + 1 : Ndof*(PE+2)] = Y;
Mbig[1:Ndof,1:Ndof] = IMW;
Fbig[1:Ndof] = F;
#

cc = PE+PN+1;
# Right layer
for i = 1 : PE  
	Mbig[Ndof*i+1 : Ndof*(i+1),Ndof*(i  )+1 : Ndof*(i+1)] = IME;

	Kbig[Ndof*i+1 : Ndof*(i+1),Ndof*(i-1)+1 : Ndof*(i  )] = D;
	Kbig[Ndof*i+1 : Ndof*(i+1),Ndof*(i  )+1 : Ndof*(i+1)] = E;
    Kbig[Ndof*i+1 : Ndof*(i+1),Ndof*(cc+(i-1)*PN+(0  ))+1 : Ndof*(cc+(i-1)*PN+1)] = CO;
	if i < PE
		Kbig[Ndof*i+1 : Ndof*(i+1),Ndof*(i+1)+1 : Ndof*(i+2)] = G;
	end
end

# Top layer.
for j = 1 : PN 
	Mbig[Ndof*(PE+j)+1 : Ndof*(PE+j+1),Ndof*(PE+j  )+1 : Ndof*(PE+j+1)] = IMN;

	if j == 1
		Kbig[Ndof*(PE+j)+1 : Ndof*(PE+j+1),1 : Ndof] = M;
	else
		Kbig[Ndof*(PE+j)+1 : Ndof*(PE+j+1),Ndof*(PE+j-1)+1 : Ndof*(PE+j  )] = M;
	end

	Kbig[Ndof*(PE+j)+1 : Ndof*(PE+j+1),Ndof*(PE+j  )+1 : Ndof*(PE+j+1)] = N;
    Kbig[Ndof*(PE+j)+1 : Ndof*(PE+j+1),Ndof*(cc+(0  )*PN+(j-1))+1 : Ndof*(cc+(0  )*PN+j  )] = CG;

	if j < PN
		Kbig[Ndof*(PE+j)+1 : Ndof*(PE+j+1),Ndof*(PE+j+1)+1 : Ndof*(PE+j+2)] = O;
	end
end


for i = 1 : PE
	for j = 1 : PN 

		Mbig[Ndof*(cc+(i-1)*PN+(j-1))+1 : Ndof*(cc+(i-1)*PN+j),
             Ndof*(cc+(i-1)*PN+(j-1))+1 : Ndof*(cc+(i-1)*PN+j  )]     = IMC;



		Kbig[Ndof*(cc+(i-1)*PN+(j-1))+1 : Ndof*(cc+(i-1)*PN+j),
			 Ndof*(cc+(i-1)*PN+(j-1))+1 : Ndof*(cc+(i-1)*PN+j  )] 	  = C;
		if i == 1
			Kbig[Ndof*(cc+(i-1)*PN+(j-1))+1 : Ndof*(cc+(i-1)*PN+j),
				 Ndof*(PE+j  )+1 			: Ndof*(PE+j+1)] 		  = CD;
		else
			Kbig[Ndof*(cc+(i-1)*PN+(j-1))+1 : Ndof*(cc+(i-1)*PN+j),
				 Ndof*(cc+(i-2)*PN+(j-1))+1 : Ndof*(cc+(i-2)*PN+j  )] = CD;
		end
		if i < PE
			Kbig[Ndof*(cc+(i-1)*PN+(j-1))+1 : Ndof*(cc+(i-1)*PN+j),
				 Ndof*(cc+(i  )*PN+(j-1))+1 : Ndof*(cc+(i  )*PN+j  )] = CG;
		end
		if j == 1
			Kbig[Ndof*(cc+(i-1)*PN+(j-1))+1 : Ndof*(cc+(i-1)*PN+j),
				 Ndof*(i)+1                 : Ndof*(i+1)] 			  = CM;
		else
			Kbig[Ndof*(cc+(i-1)*PN+(j-1))+1 : Ndof*(cc+(i-1)*PN+j),
				 Ndof*(cc+(i-1)*PN+(j-2))+1 : Ndof*(cc+(i-1)*PN+j-1)] = CM;
		end
		if j < PN
			Kbig[Ndof*(cc+(i-1)*PN+(j-1))+1 : Ndof*(cc+(i-1)*PN+j),
				 Ndof*(cc+(i-1)*PN+(j  ))+1 : Ndof*(cc+(i-1)*PN+j+1)] = CO;
		end
	end
end


#Kbig[Ndof*(PE+PN+PE*PN+1)+1:Ndof*(PE+PN+PE*PN+2),Ndof*(PE+PN+PE*PN)+1:Ndof*(PE+PN+PE*PN+1)] = IMC;

println("Resizing.")
reindx = findall(x -> x != 0, diag(Mbig));
Kbig = Kbig[reindx,reindx];
Fbig = Fbig[reindx];

println("Solving Linear Problem.")
u = Kbig\Fbig;


# Saving solution begins.
if SaveSol == 1
	println("Saving solution.")
	fp = open("points_r$r-h$hx.dat","w")
	pf = p[findind2,:];
	for i = 1 : size(pf,1)
		if p[findind2[i],1] <= xIE && p[findind2[i],2] <= xIN
    		for j = 1 : 3
		        print(fp,pf[i,j]," ")
		    end
	
		    println(fp,"")
		end
	end

	close(fp)

   dist = sqrt.((p[findind2,1] .-xstar).^2 .+(p[findind2,2] .-ystar).^2);
   g= exp.(1im*kn*dist)./(4*pi*dist);


	fr = open("sol_f$f-r$r-PE$PE-PN$PN-nLE$nLE-nLN$nLN-h$hx-real.dat","w");
	fi = open("sol_f$f-r$r-PE$PE-PN$PN-nLE$nLE-nLN$nLN-h$hx-imag.dat","w");

	fre = open("exact_f$f-r$r-h$hx-real.dat","w");
	fie = open("exact_f$f-r$r-h$hx-imag.dat","w");

	for tt = 1 : size(pf,1)
		if p[findind2[tt],1] <= xIE && p[findind2[tt],2] <= xIN
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

return
#function end
