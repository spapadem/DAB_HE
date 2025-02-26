function  Fvec(r,hx,hy,t,c,l,nodes,el,osc_f,gj)

#F = zeros((r+1)^2,1);
glp,glw = GLpw(r);


xl = nodes[el[l,1],1];
yb = nodes[el[l,1],2];

glpx = xl .+ .5*hx*(glp .+ 1);
glpy = yb .+ .5*hy*(glp .+ 1);


F = zeros((r+1)^2,1);
for i = 1 : (r+1)^2
	itomnx = findall(x -> x .< 1e-10, vec(abs.(nodes[el[l,i],1] .- glpx)));
    itomny = findall(y -> y .< 1e-10, vec(abs.(nodes[el[l,i],2] .- glpy)));

	for p = 1 : (r+1)^2
		xp = nodes[el[l,p],1];
        yp = nodes[el[l,p],2];

        p2x = findall(x -> x.<1e-10, vec(abs.(xp .- glpx)));
        q2y = findall(x -> x.<1e-10, vec(abs.(yp .- glpy)));


		F[i] = F[i] + .5*hx*hy*gj[p]*f_ex(osc_f,t,c,xp,yp)*
				(itomny[1] == q2y[1])*(itomnx[1] == p2x[1]);
		end

	end


return F
end



