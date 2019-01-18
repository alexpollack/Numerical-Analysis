function [x,v] = CGM(A)
tic
N=sqrt(size(A,1));
x = zeros(size(A,1),1);
b = zeros(size(A,1),1);
for i = 1:size(A,1)
    b(i,1) = 1;
end
r = b;
p = r;
n = 1;
alpha = 1;
while norm(r) > 10^-8 && n < 10000
    AP = A*p;
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
    n = n + 1;
    
end
toc
n
norm(r)
semilogy(v(:,1),v(:,2))
title('Semilogy for CGM')
xlabel('Iteration n')
ylabel('semlilog(x)')
tic 
R = chol(A);
toc
tic
xc = A\b;
toc
xm = reshape(x,N,N);
xcm = reshape(xc,N,N);
figure
mesh(xm)
title('Mesh of CGM x')
figure
mesh(xcm)
title('Mesh of bslashA')
end