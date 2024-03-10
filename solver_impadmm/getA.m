function cAuc = getA(S,gs,b_value,sizex,sizey)
l = size(gs,1); % number of gradient directions   g1,g2,...,gl  (example, l = 21)
A = zeros(15,l);
for i = 1:l
    A(:,i) = v3tov15(gs(i,:));
end
B = A*A';
lb = length(b_value);
b = zeros(sizex,sizey,l);
if lb ==1
    for ir = 1:sizex
        for jc = 1:sizey
            for ii = 1:l
                b(ir,jc,ii) = -log(S(ir,jc,1,ii))/b_value;
            end
        end
    end
else
    for ir = 1:sizex
        for jc = 1:sizey
            for ii = 1:l
                b(ir,jc,ii) = -log(S(ir,jc,1,ii))/b_value(ii);
            end
        end
    end
end
cAuc = zeros(15,sizex,sizey);
for ir = 1:sizex
    for jc = 1:sizey
        bt = b(ir,jc,:);
        cAuc(:,ir,jc) = B\(A*bt(:));
    end
end
end