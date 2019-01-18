 function [x,S] = sptruss(k, b)
 a = sqrt(2)/2; %alpha
if k == 1
    m = 13; %case as in question 2, one section of ABCDEFG
elseif k > 1
        m = 13 + 8*(k-1);   %more than one section
end
sf = sparse([a 0 0 -1 -a 0 0 0 0 0 0 0 0;   %top section always the same
            a 0 1 0 a 0 0 0 0 0 0 0 0]);
smid = sparse([0 1 0 0 0 -1 0 0 0 0 0 0 0; %middle section k = 1                                          
    0 0 1 0 0 0 0 0 0 0 0 0 0;              
    0 0 0 1 0 0 0 -1 0 0 0 0 0;
    0 0 0 0 0 0 1 0 0 0 0 0 0;
    0 0 0 0 a 1 0 0 -a -1 0 0 0;
    0 0 0 0 a 0 1 0 -a 0 0 0 0;
    0 0 0 0 0 0 0 1 a 0 0 -a 0;
    0 0 0 0 0 0 0 0 a 0 1 -a 0;
    0 0 0 0 0 0 0 0 0 1 0 0 -1;
    0 0 0 0 0 0 0 0 0 0 1 0 0]);
   smk = sparse([0 1 0 0 0 -1 0 0 0 0 0 0 0; %middle section for k>1 
    0 0 1 0 0 0 0 0 0 0 0 0 0;
    0 0 0 1 0 0 0 -1 0 0 0 0 0;
    0 0 0 0 0 0 1 0 0 0 0 0 0;
    0 0 0 0 a 1 0 0 -a -1 0 0 0;
    0 0 0 0 a 0 1 0 -a 0 0 0 0;
    0 0 0 0 0 0 0 1 a 0 0 -1 -1;
    0 0 0 0 0 0 0 0 a 0 1 -a 0;
    0 0 0 0 0 0 0 0 0 1 0 0 -1;
    0 0 0 0 0 0 0 0 0 0 1 0 0]);
    
S = sparse(zeros(m));           %create sparse 0 matrix to input values
S(1:2,1:13) = sf;               %input top 2 rows
S(m,m-1) = a;                   %input last row, always the same values
S(m,m) = 1;
if k == 1
    S(3:12,:) = smid;               %k = 1 original A
else
    w = 0;
    while w < k
        y = w*8+3;              %y and y1 increment the addition of values
        y1=w*8+12;                 %to the matrix as needed k > 1
        S(y:y1,(w*8 +1):(w*8+13)) = smk;        %insert values
        w = w + 1;
    end
end
x = S\b        %solve with the sparse matrix S = A for x
bandw = bandwidth(S) %compute the bandwidth of S i.e. A

 end