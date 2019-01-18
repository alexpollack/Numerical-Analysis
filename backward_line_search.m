h = 1/(N+1);
z = zeros(2*N,1);
for i = 1:2*N
z(i,1) = i*h;
end
z = z - 1;
X = zeros(2*N,1);
for i = 1:2*N
    X(i,1) = i*0.001;
end
X = X - 0.001;
alpha = 1;
rho = 0.5;
c = 10^-4;
[f, g] = geodesic(z,100,1,1);
F = zeros(2*N,1);
[fk, gk] = geodesic(X,100,1,1);
while norm(z-X) > 10^-8 
    [fk, gk] = geodesic(X,100,1,1);
    pk = -gk;
   xk = X - alpha*pk;
   X = xk;
   [fk1, gk1] = geodesic(X+alpha*pk,100,1,1);
   gk1 = -gk1;
   while fk1 > fk + c*alpha*gk'*pk
       alpha = alpha*rho;
       [fk1, gk1] = geodesic(X+alpha*pk,100,1,1);
       gk1 = -gk1;
   end
   
   %[fk, gk] = geodesic(X,100,1,1);
   norm(g)
   norm(gk)
end