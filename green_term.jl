function  green_term(x1,y1,x2,y2,kn,depth,nmod,R)
# y1 = column array
# y2 = row array
m = size(x1,1);
n = size(x1,2);
g = zeros(Complex{Float64},m,n);
g1 = zeros(m,n);
k = Int64.(Array(1:nmod));
lambda = (k*pi) / (depth);
mu = sqrt.(Complex.(kn*kn .- lambda.*lambda));

for n = k

    Xn(x) =  sin.(lambda[n]*x);
   g = g .+ 1/mu[n]*(exp.(1im*mu[n]*abs.(x1.-x2)).-
      exp.(1im*mu[n]*abs.(x1.-2*R.+x2)) ).*Xn(y1).*Xn(y2);

end

return 1im*(g/depth)

end
