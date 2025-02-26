function basisx(a,j,r)
# Compute the value of the derivative of a Lagrange polynomial L_a'(x_j).
#a: index of shape function
#j: index of point in the set of G-L points in [-1,1]


#include("GLpw.jl")

# Create Gauss-Legendre points.
glp, = GLpw(r);

Na_x = 0;

for p = 1 : r + 1
    if p !=a
        k = 1/(glp[a]-glp[p]);
        for m = 1 : r+1
            if m!=a && m != p
            	k = k * (glp[j] - glp[m])/(glp[a] - glp[m]);
            end
        end
        Na_x = Na_x + k;
    end
end

return Na_x

end
