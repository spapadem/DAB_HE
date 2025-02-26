function I2(r,hx,hy,gj,pref,elref)
K = zeros((r+1)^2,(r+1)^2);
glp,glw = GLpw(r);


i = 1;
itomn = zeros((r+1)^2,2);
#=
pref = [ -1 -1 1;
    	  1 -1 1;
		  1  1 1;
	     -1  1 1];
elref = [1 2 3 4];

pref,elref = intNodes(elref,pref,r,2,2,[-1 1],[-1 1]);
elref = elref[1,:];

gj = intWeights(r,2,2);
=#

for i = 1 : size(K,1)
	itomnx = findall(x -> x .< 1e-10, vec(abs.(pref[Int(elref[i]),1] .- glp)));
	itomny = findall(y -> y .< 1e-10, vec(abs.(pref[Int(elref[i]),2] .- glp)));
    for j = 1 : size(K,2)
		jtomnx = findall(x -> x .< 1e-10, vec(abs.(pref[Int(elref[j]),1] .- glp)));
		jtomny = findall(y -> y .< 1e-10, vec(abs.(pref[Int(elref[j]),2] .- glp)));


		for p = 1 : (r+1)^2
			xp = pref[p,1];
			yp = pref[p,2];
			
			p2x = findall(x -> x.<1e-10, vec(abs.(xp .- glp)));
			q2y = findall(x -> x.<1e-10, vec(abs.(yp .- glp)));


			K[i,j] = K[i,j] + hy/hx*gj[p]*basisx(itomnx[1],p2x[1],r)*basisx(jtomnx[1],p2x[1],r)*
										(itomny[1] == q2y[1])*(jtomny[1] == q2y[1]) +
							  hx/hy*gj[p]*basisx(itomny[1],q2y[1],r)*basisx(jtomny[1],q2y[1],r)*
                                        (itomnx[1] == p2x[1])*(jtomnx[1] == p2x[1]);
		end




    end

end

return K

end
