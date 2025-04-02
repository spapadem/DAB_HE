function  femerror(nodes,el,uh,g,hx,hy,r)


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
j = 1;


gi = zeros(Complex,size(nodes,1),1);

for i = 1 : size(nodes,1)
	pe = zeros(Complex,size(el,1),1);
	
	for l = 1 : size(el,1)
		 pt = nodes[el[l,:],1:2];
		 for j = 1 : (r+1)^2
				xp = pt[j,1];
				yp = pt[j,2];
				indx = el[l,j];
	            pe[l] = pe[l] + .25*hx*hy*gj[j]*(g[indx]*uh[i]);
#           	 err[l]  = err[l]  + .25*hx*hy*gj[j]*(ex_sol(n,t,c,xp,yp)).^2;
#       	     nuh[l] = nuh[l] + .25*hx*hy*gj[j]*(uh[indx]).^2;
    	    end
	end
	gi[i] = sum(pe);
end
#errt = sqrt(sum(err));
#nu   = sqrt(sum(nu));
#nuh  = sqrt(sum(nuh));
#pe = sqrt(sum(pe));

return gi

end

