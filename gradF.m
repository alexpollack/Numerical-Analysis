function [f,g] = geodesic(z,N,g1,g2)
h = 1/(N+1);
x = z(1:N,1);
y = z((N+1):2*N,1);


f = 0;
for i = 1:(N-1)
   s = 1 + g1*exp(-g2*(x(i,1)^2 +y(i,1)^2));
   f = f + s*(((x(i+1,1)-x(i,1))/h)^2 + ((y(i+1,1)-y(i,1))/h)^2);
end
g = zeros(1,2*N);
gx = zeros(1,N);
gy = zeros(1,N);
if nargout > 1
    gx = 0;
    for i = 1:(N-1)
        sx = (1 + g1*exp(-g2*((x(i,1)+h)^2 +y(i,1)^2)));
        xph = sx*(((x(i+1,1)-x(i,1))/h)^2 + ((y(i+1,1)-y(i,1))/h)^2);
        sx = (1 + g1*exp(-g2*((x(i,1)-h)^2 +y(i,1)^2)));
        xmh = sx*(((x(i+1,1)-x(i,1))/h)^2 + ((y(i+1,1)-y(i,1))/h)^2);
        gx(1,i) = (xph - xmh)/(2*h);
    end
    gy = 0;
    for i = 1:(N-1)
        sy = (1 + g1*exp(-g2*(x(i,1)^2 +(y(i,1)+h)^2)));
        yph = sy*(((x(i+1,1)-x(i,1))/h)^2 + ((y(i+1,1)-y(i,1))/h)^2);
        sy = (1 + g1*exp(-g2*(x(i,1)^2 +(y(i,1)-h)^2)));
        ymh = sy*(((x(i+1,1)-x(i,1))/h)^2 + ((y(i+1,1)-y(i,1))/h)^2);
        gy(1,i) = (yph - ymh)/(2*h);
    end
    g(1,1:N) = gx;
    g(1,(N+1):2*N) = gy;
end

end