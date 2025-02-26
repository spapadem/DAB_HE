close all
clear all

xs = 0.5;
ys = 0.5;
f0 = 75;
r = 0.5;
a = log(10^6)/(r*r);
b = 1;
t_c = 0.05;

hx = 0.01;
hy = hx;

T = 40;
dt = hx/0.1;
x = 0 : hx : b;
y = 0 : hy : b;

h = @(t) sqrt(2)*f0*(1-4*pi*pi*f0*f0*(t-t_c).^2)...
        .*exp(-(sqrt(2)*pi*f0*(t-t_c)).^2);
g = @(x,y) sqrt(2*a/pi)*exp(-a*((x-xs).^2+(y-ys).^2));

uex = @(t,x,y) sin(t)*sin(3*pi*x).*sin(3*pi*y);

[X,Y] = meshgrid(x,y);

t = 0 : dt : T;
for i = 1 : length(t)

%	imagesc(x,y,(uex(t(i),X,Y)))
	surf(X,Y,uex(t(i),X,Y))
	axis([0 b 0 b -1 1])
	shading interp
	colormap jet
%	colorbar
	drawnow

end
