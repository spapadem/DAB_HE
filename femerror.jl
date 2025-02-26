function  femerror(nodes,el,uh,hx,hy,r,n,t,c)


itomn = zeros((r+1)^2,2);

pref = [ -1 -1 0;
          1 -1 0;
          1  1 0;
         -1  1 0];
elref = [1 2 3 4];
    pref,elref = intNodes(elref,pref,r,2,2,[-1 1],[-1 1]);
    elref = elref[1,:];


glp,glw = GLpw(r);

gj = intWeights(r,2,2);
xj = zeros((r+1)^2,2);
j = 1;



err = zeros(size(el,1),1);
nu  = zeros(size(el,1),1);
nuh = zeros(size(el,1),1);

for l = 1 : size(el,1)
	 pt = nodes[el[l,:],1:2];
	 for j = 1 : (r+1)^2
			xp = pt[j,1];
			yp = pt[j,2];
			indx = el[l,j];
            err[l] = err[l] + .25*hx*hy*gj[j]*(ex_sol(n,t,c,xp,yp) - uh[indx]).^2;
            nu[l]  = nu[l]  + .25*hx*hy*gj[j]*(ex_sol(n,t,c,xp,yp)).^2;
            nuh[l] = nuh[l] + .25*hx*hy*gj[j]*(uh[indx]).^2;
        end
end

errt = sqrt(sum(err));
nu   = sqrt(sum(nu));
nuh  = sqrt(sum(nuh));


return [errt,nu,nuh]

end

