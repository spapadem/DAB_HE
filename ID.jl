function ID(p)

# Labels.
# Bottom  : 1
# Left    : 2
# Top     : 3
# Right   : 4
# Interior: 5

# Dirichlet b.c.
D = [1 2 3];
# Neummann b.c.
#N = [1 0 1 1];
# Periodic b.c.
P = [0 0];
DOFS = zeros(size(p,1),1);
ind = 1;

for i = 1 : size(p,1)
	if (findall(x->x.!=0, vec(p[i,3] .== D)) .!= 0) != []
		DOFS[i] = 0;
	else
		if P == [1 3]
			if p[i,3] == 3
				inds1 = findall(x -> x == 1, p[:,3]);
				xind = findall(x -> x .< 1e-10, vec(abs.(p[i,1] .- p[inds1,1])));
				if xind != []
#					display(xind)
#					display(inds1[xind])
#					display(DOFS[inds1[xind]])
#					display(DOFS[inds1[xind]][1])
					DOFS[i] = DOFS[inds1[xind][1]];
					continue;
				end
			end
		end
		DOFS[i] = ind;
		ind = ind + 1;
	end
	
end

return DOFS

end

