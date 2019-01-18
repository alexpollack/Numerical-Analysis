function [Q,R] = givensQR(A)
[m, n] = size(A); %gets size of A
%Initializes a matrix for Q since Q will be all of 
%G_i multiplied, gives a space to begin this multiplication
Q = eye(m);
%Q = zeros(2,n);
%R is a modified verion of A each iteration, starting 
%from the initial value of A
R = A; 
k = 1;
for j = 1:n
    for i = m:-1:(j+1)
        %Space for a new G each time we need to eliminate 
        %a value, which here is reset for each new elimination
        G = eye(m);
        %Creates a vector a from the current
        %R to make the givens matrix
        a = [R(i-1, j) R(i, j)];
        %Uses the vector a from R in the function for computing J to  
        %create the newportion of G that will eliminate the current slot.
        J = compJ(a);
        %Q is a matrix that just stores the values of c and s from
        % each J used to create a new givens rotation in each column of Q
        %Q is here is only the saved values of each c and s, which can be
        %used to find the actual Q matrix
        %Q(1,k) = J(1,1);
        %Q(2,k) = J(1,2);
        
        k = k+1;
        %inserts J within G
        G([i-1, i], [i-1,i]) = J;
        %R is computed by multiplication with each new G 
        %to eliminate each entry under a
        R = G*R;  
        Q =Q*G';
        
        
    end
end     
end

