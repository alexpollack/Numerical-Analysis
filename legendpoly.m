function [Q,R] = legendpoly(n,m,method)
    x = linspace(-1,1,m)';
    A = zeros(m,n);
    for i = 1:n
            A(:,i) = x.^(i-1);
    end
    if method == 'Q'
        [Q,R] = qr(A,0);
    elseif method == 'G'
        [Q,R] = clgs(A);
    elseif method == 'M'
        [Q,R] = mgs(A);
    end
    scale = Q(m,:);
    Q = Q*diag(1 ./scale);
    hold on
    xlabel('mth Index of the polynomial')
    ylabel('Value of x^n normalized between [-1, 1]')
    
    if method == 'Q'
        plot(Q)
        title(['First ',num2str(n),' Legendre Polynomials Using QR()'])
    elseif method == 'G'
        plot(Q, '--')
        title(['First ',num2str(n),' Legendre Polynomials Using CLGS()'])
    elseif method == 'M'
        plot(Q, ':')
        title(['First ',num2str(n),' Legendre Polynomials Using MGS()'])
    end
    
end

