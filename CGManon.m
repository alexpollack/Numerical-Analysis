function [x,v] = CGManon(A,B)
tic
N=sqrt(size(A,1));
NA = size(A,1);
x = zeros(NA,1);
b = zeros(NA,1);
for i = 1:NA
    b(i,1) = 1;
end
r = b;
p = r;
n = 1;
while norm(r) > 10^-8 && n < 10000
    Avmult = @(A,B,p) p + A\B*(B'*(A\p));
    AP = Avmult(A,B,p);
    alpha = (r'*r)/(p'*AP);
    if alpha < 0 
        v(n,1) = n;
        v(n,2) = norm(r);
        break
    end
    x = x + alpha*p;
    rn = r - alpha*AP;
    beta = (rn'*rn)/(r'*r);
    r = rn;
    p = r + beta*p;
    v(n,1) = n;
    v(n,2) = norm(r);
    n = n +1;
    
end
toc
n
norm(r)
semilogy(v(:,1),v(:,2))
title('Semilogy for CGM')
xlabel('Iteration n')
ylabel('semlilog(x)')
xm = reshape(x,N,N);
figure
mesh(xm)
title('Mesh of CGManon x')
end