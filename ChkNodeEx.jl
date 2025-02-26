function ChkNodeEx(p,pr,ind,hx,hy)

# Check if a Node already exists in our mesh. If found, also return where it it located.

indx = 0;
px = p[1];
py = p[2];
AlEx = 0;
ind1 = findall(x -> x .< 1e-6*hx, vec(abs.(px .- pr[1:ind-1,1])));
ind2 = findall(y -> y .< 1e-6*hy, vec(abs.(py .- pr[1:ind-1,2])));


flag = 0;
for m = 1 : length(ind1)
	for n = 1 : length(ind2)
		if (ind1[m] == ind2[n])
			AlEx = 1;
			indx = ind1[m];
			flag = 1;
			break;
		end
	end
	if flag == 1;
		break
	end
end

#println(AlEx)
#println(indx)
return [AlEx,indx]

end
