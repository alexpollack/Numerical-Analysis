function [x,v] = PCCGM(A,D)
tic
N=sqrt(size(A,1));
x = zeros(size(A,1),1);
b = zeros(size(A,1),1);
for i = 1:size(A,1)
    b(i,1) = 1;
end
r = b;
p = r;
z = inv(D)*r;
n = 1;
while norm(r) > 10^-8 && n < 10000
    AP = A*p;
    alpha = (r'*r)/(p'*AP);
    x = x + alpha*p;
    rn = r - alpha*AP;
    if norm(rn) < 10^-8
        break
    end
    zn = inv(D)*rn;
    beta = (zn'*rn)/(z'*r);
    r = rn;
    z = zn;
    p = z + beta*p;
    v(n,1) = n;A
    v(n,2) = norm(r);
    n = n +1;
    
end
n
norm(r)
toc