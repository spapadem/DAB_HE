function  norm_calc(nodes,el,uh,uph,hx,hy,r,c)


itomn = zeros((r+1)^2,2);

pref = [ -1 -1 1;
          1 -1 1;
          1  1 1;
         -1  1 1];
elref = [1 2 3 4];
pref,elref = intNodes(elref,pref,r,2,2,[-1 1],[-1 1]);
elref = elref[1,:];

glp,glw = GLpw(r);

gj = intWeights(r,2,2);
xj = zeros((r+1)^2,2);


err = zeros(size(el,1),1);
nu  = 0; 
nug = 0;
for l = 1 : size(el,1)
	 pt = nodes[el[l,:],1:2];
	 for j = 1 : (r+1)^2
		xp = pt[j,1];
		yp = pt[j,2];
		indx = el[l,j];
	    jtomnx = findall(x -> x .< 1e-10, vec(abs.(pref[Int(elref[j]),1] .- glp)));
    	jtomny = findall(y -> y .< 1e-10, vec(abs.(pref[Int(elref[j]),2] .- glp)));
		guhx = 0;
		guhy = 0;
		for k = 1 : (r+1)^2
			
		    ktomnx = findall(x -> x .< 1e-10, vec(abs.(pref[Int(elref[k]),1] .- glp)));
	    	ktomny = findall(y -> y .< 1e-10, vec(abs.(pref[Int(elref[k]),2] .- glp)));
            guhx = guhx + uh[el[l,k]]*basisx(ktomnx[1],jtomnx[1],r)*(ktomny[1]==jtomny[1]);
            guhy = guhy + uh[el[l,k]]*basisx(ktomny[1],jtomny[1],r)*(ktomnx[1]==jtomnx[1]);	
		end
    	nu  = nu  + .25*hx*hy*gj[j]*(uph[indx]).^2 
		nug = nug + gj[j]*((guhx)^2+(guhy)^2);
	end
end
nu   = sqrt(nu);
nug  = sqrt(nug)
return [nu,nug]

end

