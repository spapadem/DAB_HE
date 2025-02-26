function GLpw(r)
# Create Gauss-Legendre integration points.


r1  = r+1;
glp = - cos.(pi*(0:r)/r);
P=zeros(r1,r1);
xold=2;

while maximum(abs.(glp .- xold)) > eps()
	xold    = glp;
	P[:,1] .=   1;
	P[:,2] .= glp;
	for k = 2 : r
		P[:,k+1] .=((2*k-1)*glp.*P[:,k]-(k-1)*P[:,k-1] )/k;
	end
	glp = xold-( glp.*P[:,r1]-P[:,r])./( r1*P[:,r1] );
end

w = 2 ./(r*r1*P[:,r1].^2);

return [glp,w]

end
