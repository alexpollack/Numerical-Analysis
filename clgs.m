function [Q,R] = clgs(A)
% classical Gram-Schmidt (unstable) to compute 
%    "reduced" QR factorization of A
% input: A (m by n, with m >= n)
% output: Q (orthogonal columns, m by n)
%         R (upper triangular, n by n)
% properties: A = Q*R (in exact arithmetic)
% reference: Trefethen and Bau, Alg 7.1
[m,n] = size(A);
% check dimensions are valid
if m < n
   error('invalid dimension for A')
end
% compute QR via Gram-Schmidt
for j = 1:n 
   % at jth step through loop, we take jth column of A and
   % add multiples of previous columns of Q to it to construct
   % jth column of Q with property that 
   %   - the first j columns of Q and A span same space
   %   - the first j columns of Q are mutually orthonormal 
   v = A(:,j);   
   for i=1:j-1   % does nothing if i=j
      R(i,j) = Q(:,i)'*A(:,j);
      v = v - R(i,j)*Q(:,i); % constructed so v is orthogonal
         % to all columns of Q generated so far
   end % closes inner for loop
   R(j,j) = norm(v);  % could be + or -
   Q(:,j) = v/R(j,j); % normalize
end  % closes outer for loop
end