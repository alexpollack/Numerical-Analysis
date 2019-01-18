format long %digits(6) %Because computed actual values are only 6 places
a = 0; %part a
%a = 10^-3; %part b
%a = 10^-6; %part c
%a = 0; %part d
b = 1;
pass = 0; %did it meet accuracy
actual = 0.2204582194 %part a
%actual = 0.6696579045  %part b
%actual = 0.6697331726  %part c
%actual = 0.6697332001  %part d
 m = 16; %initial m
%while pass ~= 1 %while not at best accuracy
   Tk = zeros(m+1,m+1); %initialize the table for current m
   mm = m; %dummy copy of m for loop
   for i = 1:(m+1) %compute T0
        Tk(i,1) = (4*trapezoid(a,b,2*mm) - trapezoid(a,b,mm))/3;
        mm = 2*mm;
    end
    for k = 1:m %compute Tk's
        for i = 1:(m+1-k)
            Tk(i,k+1) = ((4^k)*Tk(i+1,k) - Tk(i,k))/((4^k)-1);
        end
    end
    %Tk(:,:) = vpa(Tk(:,:));
    for i = 1:(m+1) %checks best aproximation
        for j = 1:(m+1)
            if Tk(i,j) ~= 0 && abs(Tk(i,j) - actual) <= 10^-11
                pass = 1; %if accuracy met
            end
        end
    end
    %if pass ~= 1 %if accuracy not met, try again with another m
     %   m = m*2; 
    %end
%end
Tm = Tk(1,m+1) %best approx found
m %in m rounds
error = abs(actual - Tk(1,m+1)) %how close to actual integral value
