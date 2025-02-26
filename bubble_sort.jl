function bubble_sort(p,tol)
# Sort the vector p, comparing the values of p up to a tolerance tol.
# We return the sorted vector and the vector containing the index swaps.

sort_tol = 1e6;
n = length(p);
# Original indices of the unsorted vector p.
inds = Array(1 : n);
psort = zeros(size(p));

for j = 1 : n-1
#	Comparing each number with the next and swapping. We check a floor version of the values,
#	up to the number of digits, defined by the exponent of sort_tol.
	for i = 1 : n-1
		if (floor(sort_tol*p[i])/sort_tol)>(floor(sort_tol*p[i+1])/sort_tol);
			temp      = p[i];
			p[i]      = p[i+1];
			p[i+1]    = temp;	
			temp      = inds[i];
			inds[i]   = inds[i+1];
			inds[i+1] = temp;
		end
	end
end

return [p,inds]

end

