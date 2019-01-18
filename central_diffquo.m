format long;
n = 1;
x = 1;
h = 1;
deriv = cos(x);
cdiffquo = 0;
error = 0;
fprintf('Deriv = %f\n', deriv)
%h range form 10^{-1} tp 10^{-20}
fprintf('h\t\t cdiffquo\t error\n')
while n <= 20
    h = h/10; %decrementing h
    cdiffquo = (sin(x+h) - sin(x-h))/(2*h); %Implementation of cent. diff.
    error = abs(deriv - cdiffquo);  %error results
    fprintf('%1.1e \t %1.1e \t %1.1e\n', h, cdiffquo, error);
    n = n + 1;
end


