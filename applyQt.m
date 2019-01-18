function [v] = applyQt(A,b)
[m,n] = size(A);
%v is storing all of the c and s values from the givens of A
v = givensQR(A);
Q = eye(m);
J = zeros(2);
k = 1;
for j = 1:n
    for i = m:-1:(j+1)
        %this essentially does the multiplication that would have happened
        %were Q updated every iteration within givens,
        G = eye(m);
        %'J' is recreated, and reinserted in a matrix to be multiplied to
        %find Q
        J = [v(1,k) v(2,k); -v(2,k) v(1,k)];
        G([i-1, i], [i-1,i]) = J;
        Q = Q*G';
        k = k+1;
    end
end
%Out put the vector matrix product
v = Q'*b;
end