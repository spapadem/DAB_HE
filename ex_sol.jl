function  ex_sol(n,t,c,x,y)
#	u = t^2*sin(n*pi*x).*sin(n*pi*y);
	return t^2*sin.(n*pi*x).*sin.(n*pi*y);
#	return sin.(n*pi*x).*sin.(n*pi*y);
#	u = sin(t)*ones(size(x));
end
