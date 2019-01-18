function [f,g,H] = geodesic(z,N,g1,g2)
h = 1/(N+1);
x = z(1:N,1);
y = z(1:N,1);



f = 0;
for i = 1:(N-1)
   s = 1 + g1*exp(-g2*(x(i,1)^2 +y(i,1)^2));
   f = f + h*s*(((x(i+1,1)-x(i,1))/h)^2 + ((y(i+1,1)-y(i,1))/h)^2);
end
g = zeros(2*N,1);
gx = zeros(N,1);
gy = zeros(N,1);
if nargout > 1
    %gx = 0;
    for i = 1:N%(N-1)
        if i == N
            sx = h*(1 + g1*exp(-g2*((x(i,1)+h)^2 +y(i,1)^2)));
            xph = sx*(((x(i,1)-x(i-1,1))/h)^2 + ((y(i,1)-y(i-1,1))/h)^2);
            sx = h*(1 + g1*exp(-g2*((x(i,1)-h)^2 +y(i,1)^2)));
            xmh = sx*(((x(i,1)-x(i-1,1))/h)^2 + ((y(i,1)-y(i-1,1))/h)^2);
            gx(i,1) = (xph - xmh)/(2*h);
        else
        sx = h*(1 + g1*exp(-g2*((x(i,1)+h)^2 +y(i,1)^2)));
        xph = sx*(((x(i+1,1)-x(i,1))/h)^2 + ((y(i+1,1)-y(i,1))/h)^2);
        sx = h*(1 + g1*exp(-g2*((x(i,1)-h)^2 +y(i,1)^2)));
        xmh = sx*(((x(i+1,1)-x(i,1))/h)^2 + ((y(i+1,1)-y(i,1))/h)^2);
        gx(i,1) = (xph - xmh)/(2*h); 
        end
       
    end
    %gy = 0;
    for i = 1:N %(N-1)
        if i == N
        sy = h*(1 + g1*exp(-g2*(x(i,1)^2 +(y(i,1)+h)^2)));
        yph = sy*(((x(i,1)-x(i-1,1))/h)^2 + ((y(i,1)-y(i-1,1))/h)^2);
        sy = h*(1 + g1*exp(-g2*(x(i,1)^2 +(y(i,1)-h)^2)));
        ymh = sy*(((x(i,1)-x(i-1,1))/h)^2 + ((y(i,1)-y(i-1,1))/h)^2);
        gy(i,1) = (yph - ymh)/(2*h);
        else
            sy = h*(1 + g1*exp(-g2*(x(i,1)^2 +(y(i,1)+h)^2)));
        yph = sy*(((x(i+1,1)-x(i,1))/h)^2 + ((y(i+1,1)-y(i,1))/h)^2);
        sy = h*(1 + g1*exp(-g2*(x(i,1)^2 +(y(i,1)-h)^2)));
        ymh = sy*(((x(i+1,1)-x(i,1))/h)^2 + ((y(i+1,1)-y(i,1))/h)^2);
        gy(i,1) = (yph - ymh)/(2*h);
        end
    end
    g(1:N,1) = gx;
    g((N+1):2*N,1) = gy;
end
if nargout > 2
u0 = zeros(N,1); 
for i = 1:N
    if i == N 
        s = 1 + g1*exp(-g2*(x(i,1)^2 +y(i,1)^2));
        u0(i,1) = h*s*(((x(i,1)-x(i-1,1))/h)^2 + ((y(i,1)-y(i-1,1))/h)^2);
    else
        s = 1 + g1*exp(-g2*(x(i,1)^2 +y(i,1)^2));
        u0(i,1) = h*s*(((x(i+1,1)-x(i,1))/h)^2 + ((y(i+1,1)-y(i,1))/h)^2);
    end
end
H = zeros(N,1);
for i = 1:N
    for j = 1:N
        if i == j && i ~= 1 && i ~= N
            H(i,i) = 2*((2*u0(i,1) - 1)^2 - 1);
            H(i,i-1) = -((2*u0(i-1,1) - 1)^2 - 1);
            H(i,i+1) = -((2*u0(i+1,1) - 1)^2 - 1);
        elseif i == 1
            H(i,i) = 2*((2*u0(i,1) - 1)^2 - 1);
            H(i,i+1) = -((2*u0(i+1,1) - 1)^2 - 1);
        elseif i == N
            H(i,i-1) = -((2*u0(i-1,1) - 1)^2 - 1);
            H(i,i) = 2*((2*u0(i,1) - 1)^2 - 1);
        end %F' is the tridiagonal matrix created from u''
    end
end
H = H * (2/h^2);
end


end