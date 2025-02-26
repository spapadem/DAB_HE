function I5(r,hx,hy,gj,pref,elref)
I = zeros((r+1)^2,(r+1)^2);
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
=#
yR = zeros(r+1,1);
ind = 1;
for i = 1 : size(pref,1);
    if abs(pref[i,1] - pref[2,1]) < 1e-6
        yR[ind] = pref[i,2];
        ind = ind + 1;
    end
end

for i = 1 : size(I,1)
	itomnx = findall(x -> x .< 1e-10, vec(abs.(pref[Int(elref[i]),1] .- glp)));
	itomny = findall(y -> y .< 1e-10, vec(abs.(pref[Int(elref[i]),2] .- glp)));
    for j = 1 : size(I,2)
		jtomnx = findall(x -> x .< 1e-10, vec(abs.(pref[Int(elref[j]),1] .- glp)));
		jtomny = findall(y -> y .< 1e-10, vec(abs.(pref[Int(elref[j]),2] .- glp)));

		for p = 1 : r+1
			yp = yR[p];
			
			q2y = findall(x -> x.<1e-10, vec(abs.(yp .- glp)));

			I[i,j] = I[i,j] + hy/hx*glw[q2y[1]]*basisx(jtomnx[1],r+1,r)*(jtomny[1] == q2y[1])*
												 basis(itomnx[1],r+1,r)*(itomny[1] == q2y[1]);
		end
    end
end

return I

end
