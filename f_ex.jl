function f_ex(n,t,c,x,y)

	return 2*sin(n*pi*x).*sin(n*pi*y)*(1+(n^2)*(c^2)*(t^2)*(pi^2));
#	return 2*sin.(n*pi*x).*sin.(n*pi*y)*n^2*c^2*pi^2;
#	u = -sin(t);
end
