function basis(a,j,r)
# Compute the value of a Lagrange polynomial L_a(x_j).
#a: index of shape function
#j: index of point in the set of G-L points in [-1,1]

# Create Gauss-Legendre points.
glp, = GLpw(r);

Na = 1;

for p = 1 : r + 1
        if p==a
            continue;
        end
        Na = Na * ((glp[j]-glp[p]))/(glp[a]-glp[p]);
end

return Na;

end
