function [T] = trapezoid(a,b,m)
    h = (b - a)/m;  %spacing h
    x = zeros(m+1,1); %initialize nodes
    x(1,1) = a; %first node is bound a
    for i = 1:m 
        x(i+1,1) = a + i*h; %set nodes 
    end
    sum = 0; %initialize summation value 
    for i = 1:(m+1)
       if i == 1 || i == (m+1) %a, b in [a,b] have factor 1/2
           sum = sum + 0.5*exp(-2*x(i,1))/(1+4*x(i,1)); %f part a
           %sum = sum + 0.5*sin(x(i,1)^(1/3)); %f part b,c,d
       else %all other points
           sum = sum + exp(-2*x(i,1))/(1+4*x(i,1)); %f part a
           %sum = sum + sin(x(i,1)^(1/3));  %f part b,c,d
       end
    end
    T = h*sum; %final sum times factor h
end

