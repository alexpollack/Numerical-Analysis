function [Q,R] = mgs(A)
% modified Gram-Schmidt to compute 
%    "reduced" QR factorization of A
% input: A (m by n, with m >= n)
% output: Q (orthogonal columns, m by n)
%         R (upper triangular, n by n)
% properties: A = Q*R (in exact arithmetic)
%
[m,n] = size(A);
% check dimensions are valid
if m < n
   error('invalid dimension for A')
end
% compute QR via Gram-Schmidt
for i = 1:n
    v = A(:,i);   %v_i is the ith column of A for each iteration
    R(i,i) = norm(v);  % could be + or -
    Q(:,i) = v/R(i,i); % normalize
    for j = (i + 1):n
        R(i,j) = Q(:,i)'*A(:,j);
       A(:,j) = A(:,j) - R(i,j)*Q(:,i); % constructed so v is orthogonal
    end
end  % closes outer for loop
end