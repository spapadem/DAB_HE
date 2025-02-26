function intNodes(el,p,r,hx,hy,x,y)


#include("GLpw.jl")
#include("ChkNodeEx.jl")
xsc = NaN;
ysc = NaN;
b = NaN;


Nel = size(el,1);
N = size(p,1);
glp, = GLpw(r);
pr = zeros(Nel*(r+1)^2,3);
#pr[1,:] .= NaN;

ind = 1;
elr = zeros(size(el,1),(r+1)^2);
for l = 1 : Nel
	for indel = 1 : 4
		p1 = p[Int(el[l,indel]),:];
			pr[ind,:] = p1;
			elr[l,indel] = ind;
			ind = ind + 1;
	end
	indel = 5;
if r > 1
	for k = 1 : r+1
		for j = k : r+1-k
			pt = p[Int(el[l,1]),1:2];
			px = .5*hx*(glp[j]+1) .+ pt[1];
		    py = pt[2] +.5*hy*(glp[k]+1);
				pr[ind,1:2] = [px py];
                elr[l,indel] = ind;
				indel = indel + 1;
				if(px == x[1])
		                    plab = 2;
        		elseif(px == x[end])
                		    plab = 4;
		        elseif(py == y[1])
        		            plab = 1;
                elseif(py == y[end])
	                	    plab = 3;
				elseif(px == xsc && py <= (ysc +b/2) && py >=(ysc-b/2));
							plab = 6;
	        	else
		        	plab = 5;
				end
	                pr[ind,3] = plab;
        	        ind = ind + 1;
		end
		for i = k : r+1-k
        	        pt = p[Int(el[l,2]),1:2];
                	px = pt[1] - .5*hx*(glp[k]+1);
	                py = .5*hy*(glp[i]+1) + pt[2];
	                        pr[ind,1:2] = [px py];
        	                elr[l,indel] = ind;
                	        indel = indel + 1;
                        	if(px == x[1])
	                            plab = 2;
        	                elseif(px == x[end])
                	            plab = 4;
                        	elseif(py == y[1])
	                            plab = 1;
        	                elseif(py == y[end])
                	            plab = 3;
							elseif(px == xsc && py <= (ysc +b/2) && py >=(ysc-b/2));
							plab = 6;
                        	else
	                            plab = 5;
        	                end
                	        pr[ind,3] = plab;
                        	ind = ind + 1;
        	end
	        for j = k : r+1-k
        	        pt = p[Int(el[l,3]),1:2];
                	px = -.5*hx*(glp[j]+1) + pt[1];
	                py = pt[2] - .5*hy*(glp[k]+1);
                        	pr[ind,1:2] = [px py];
	                        elr[l,indel] = ind;
        	                indel = indel + 1;
                	        if(px == x[1])
                        	    plab = 2;
	                        elseif(px == x[end])
        	                    plab = 4;
                	        elseif(py == y[1])
                        	    plab = 1;
	                        elseif(py == y[end])
        	                    plab = 3;
							elseif(px == xsc && py <= (ysc +b/2) && py >=(ysc-b/2));
								plab = 6;
                	        else
                        	    plab = 5;
	                        end
        	                pr[ind,3] = plab;
                	        ind = ind + 1;
        	end
		for i = k : r+1-k
        	        pt = p[Int(el[l,4]),1:2];
                	px = pt[1] + .5*hx*(glp[k]+1);
	                py = -.5*hy*(glp[i]+1) + pt[2];
	                        pr[ind,1:2] = [px py];
        	                elr[l,indel] = ind;
                	        indel = indel + 1;
                        	if(px == x[1])
	                            plab = 2;
        	                elseif(px == x[end])
                	            plab = 4;
                        	elseif(py == y[1])
	                            plab = 1;
        	                elseif(py == y[end])
								plab = 3;
							elseif(px == xsc && py <= (ysc +b/2) && py >=(ysc-b/2));
								plab = 6;
                        	else
	                            plab = 5;
        	                end
                	        pr[ind,3] = plab;
                        	ind = ind + 1;
        	end

	end
end

    if mod(r,2)==0
        pr[ind,1:2] = [(pr[Int(elr[l,1]),1]+pr[Int(elr[l,2]),1])/2, (pr[Int(elr[l,1]),2]+pr[Int(elr[l,3]),2])/2];
        pr[ind,3]   = 5;
        elr[l,indel] = ind;
        ind = ind + 1;
        indel = indel + 1;
    end

end
pr = pr[1:maximum(findall(x->x!=0,pr[:,3])),:];

for i = 1 : size(pr,1)
	for j = i : size(pr,1)
		if i != j
			if pr[i,1] == pr[j,1] && pr[i,2] == pr[j,2]
#				prem += 1
				pr[j,3] = 0;
				elind = findall(x->x==j,elr)
#				println(elind)
				if elind != []
					elr[elind[1][1],elind[1][2]] = i
#					elr[elind[1][1],elind[1][2]+1:end] = elr[elind[1][1],elind[1][2]+1:end] .- 1;
				end
			end
		end
	end
end

pr = pr[findall(x->x!=0,pr[:,3]),:];


prem = 0;

for i = 1 : size(pr,1)
	j = i;
	elind = []
	while elind == []
		elind = findall(x->x==j,elr[:])
		j = j + 1;
	end
	println(elind[1])
	elr[elind] .= i;
end

#elr = elr .- (elr[end]-size(pr,1)).*(elr .> size(pr,1));
return [pr,elr]

end
