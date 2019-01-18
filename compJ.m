function [J] = compJ(v)
c = v(1)/norm(v);   %computes the 'cos' compontent from first entry of v
s = v(2)/norm(v);   %computes the 'sin' compontent from second entry of v
J = [c s; -s c];    %matrix J to be returned
end