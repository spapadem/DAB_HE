function intWeights(r,hx,hy)

#include("GLpw.jl")
#include("ChkNodeEx.jl")

p = [ -1 -1 0;
       1 -1 0;
       1  1 0;
      -1  1 0];
el = [1 2 3 4];

intW = zeros((r+1)^2,1);
Nel = 1;
N = size(p,1);
glp,glw = GLpw(r);
pr = zeros(Nel*(r+1)^2,3);
pr[1,:] .= NaN;

ind = 1;
elr = zeros(1,(r+1)^2);
intW[1] = glw[1]*glw[1];
intW[2] = glw[1]*glw[end];
intW[3] = glw[end]*glw[1];
intW[4] = glw[end]*glw[end];
for l = 1 : Nel
	for indel = 1 : 4
		p1 = p[el[l,indel],:];
		AlEx,ind1 = ChkNodeEx(p1,pr,ind,hx,hy);
		if AlEx == 0
			pr[ind,:] = p1;
			elr[l,indel] = ind;
			ind = ind + 1;
		else
			elr[l,indel] = ind1;
		end	
	end
	indel = 5;
	for k = 1 : r+1
		
		for j = 1 : r+1
			pt = p[el[l,1],1:2];
			px = pt[1] + .5*hx*(glp[j]+1);
		    py = pt[2] + .5*hy*(glp[k]+1);
			itomnx = findall(vec(abs.(px .- glp) .< 1e-10));
		    itomny = findall(vec(abs.(py .- glp) .< 1e-10));
			AlEx,ind1 = ChkNodeEx([px py],pr,ind,hx,hy);
			if AlEx == 0
                pr[ind,1:2] = [px py];
       	        elr[l,indel] = ind;
				intW[indel] = glw[itomnx[1]]*glw[itomny[1]];
				indel = indel + 1;
				ind = ind + 1;
            else
       	        if maximum(ind1 .== elr[l,:]) == false  		
					elr[l,indel] = ind1;
					indel = indel + 1;
				end
			end
		end
		for i = 1 : r+1
   	        pt = p[el[l,2],1:2];
           	px = pt[1] - .5*hx*(glp[k]+1);
            py = .5*hy*(glp[i]+1) + pt[2];
			itomnx = findall(vec(abs.(px .- glp) .< 1e-10));
		    itomny = findall(vec(abs.(py .- glp) .< 1e-10));
        	AlEx,ind1 = ChkNodeEx([px py],pr,ind,hx,hy);
            if AlEx == 0
				pr[ind,1:2] = [px py];
        	    elr[l,indel] = ind;
				intW[indel] = glw[itomnx[1]]*glw[itomny[1]];
                indel = indel + 1;
                ind = ind + 1;
			else
				if maximum(ind1 .== elr[l,:]) == false
        	    	elr[l,indel] = ind1;
					indel = indel + 1;
				end
            end
		end
	    
	    for j = 1 : r+1
        	    pt = p[el[l,3],1:2];
                px = -.5*hx*(glp[j]+1) + pt[1];
	            py = pt[2] - .5*hy*(glp[k]+1);
				itomnx = findall(vec(abs.(px .- glp) .< 1e-10));
			    itomny = findall(vec(abs.(py .- glp) .< 1e-10));
        	    AlEx,ind1 = ChkNodeEx([px py],pr,ind,hx,hy);
                if AlEx == 0
                    pr[ind,1:2] = [px py];
	                elr[l,indel] = ind;
					intW[indel] = glw[itomnx[1]]*glw[itomny[1]];
        	        indel = indel + 1;
                	ind = ind + 1;
	            else	
					if maximum(ind1 .== elr[l,:]) == false
        	          	elr[l,indel] = ind1;
						indel = indel + 1;
					end
	            end	
        	end
			for i = 1 : r+1
				pt = p[el[l,4],1:2];
            	px = pt[1] + .5*hx*(glp[k]+1);
	        	py = -.5*hy*(glp[i]+1) + pt[2];
				itomnx = findall(vec(abs.(px .- glp) .< 1e-10));
			    itomny = findall(vec(abs.(py .- glp) .< 1e-10));
        		AlEx,ind1 = ChkNodeEx([px py],pr,ind,hx,hy);
            	if AlEx == 0
	            	pr[ind,1:2] = [px py];
        	        elr[l,indel] = ind;
					intW[indel] = glw[itomnx[1]]*glw[itomny[1]];
                	indel = indel + 1;
                    ind = ind + 1;
	            else
					if maximum(ind1 .== elr[l,:]) == false
        	          	elr[l,indel] = ind1;
						indel = indel + 1;
					end
	            end
			end

	end

end

intW = intW[1:(r+1)^2];

return intW

end
