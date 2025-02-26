function  Mesh2Drect(x,y)

xsc = 6;
ysc = 1.5;
b = 1;


# Length of discretizations.
Nx = length(x)-1;
Ny = length(y)-1;
# Create points
p = zeros((Nx+1)*(Ny+1),3);
ind = 1;
for j = 1 : Ny + 1
    for i = 1 : Nx + 1
        p[ind,1:2] = [x[i] y[j]];
        # Assign labels.
		# Bottom  : 1
		# Left    : 2
		# Top     : 3
		# Right   : 4
		# Interior: 5
#		if(x(i) <= (xsc+b/2) &&  x(i)>=(xsc-b/2) && y(j) <= (ysc +b/2) && y(j) >=(ysc-b/2));
#                            plab = 6;
        if(x[i] == x[1])
            plab = 2;
        elseif(y[j] == y[1])
            plab = 1;
        elseif(y[j] == y[end])
            plab = 3;
        elseif(x[i] == x[end])
            plab = 4;
        else
            plab = 5;
        end
        p[ind,3] = plab;
        ind = ind + 1;  
    end
    
end

#Create edges
e = zeros(4*Nx*Ny,3);
#Horizontal edges first
ind = 1;
for j = 1 : Ny+1
    for i = 1 : Nx
        e[ind,1:2] = [(j-1)*(Nx+1) + i, (j-1)*(Nx+1) + (i + 1)];
        ind = ind + 1;
    end
end

#Vertical edges now.
for i = 1 : Nx +1
    for j = 1 : Ny
        e[ind,1:2] = [(j-1)*(Nx+1)+i, j*(Nx+1)+i];
        ind = ind + 1;
    end
end

# Assign edge labels.
for i = 1 :ind-1
   ed = e[i,1:2];
   p1 = p[Int(ed[1]),1:2];
   p2 = p[Int(ed[2]),1:2];
   if(p1[1] == x[1] && p2[1] == x[1])
       elab = 2;
   elseif (p1[1] == x[end] && p2[1] == x[end])
       elab = 4;
   elseif(p1[2] == y[1] && p2[2] == y[1])
       elab = 1;
   elseif(p1[2] == y[end] && p2[2] == y[end])
       elab = 3;
   else
       elab = 5;
   end
   e[i,3] = elab;
end
e = e[1:ind-1,:];
#Create elements
Nel = Nx*Ny;
el = zeros(Nel,4);

for i = 1 : Nel
    rem = floor((i-1)/(Nx));
    el[i,:] = [i+rem, i+1+rem, i+Nx+2+rem,i+Nx+1+rem];
end
    
    
return [p,e,el]
    
end

    
    
    
    
